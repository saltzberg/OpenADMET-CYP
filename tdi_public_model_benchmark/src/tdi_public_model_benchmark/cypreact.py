"""ID-preserving batch adapter for the public CypReact v1.2 R/N output."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem

ENDPOINTS = ["CYP2D6", "CYP3A4"]
PUBLIC_COLUMNS = {"CYP2D6": "2D6", "CYP3A4": "3A4"}
BUNDLE_URL = "https://github.com/Le0nT1/CypReact_old/releases/tag/1.2"
BUNDLE_ZIP_SHA256 = "2778b65dc786125a9f6af489d599f6de6bbe81952f3d484285d0dc73370da48f"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_cypreact_csv(path: Path) -> pd.DataFrame:
    """Parse the historical CSV, retaining its coarse public R/N semantics."""
    raw = pd.read_csv(path, skipinitialspace=True)
    raw.columns = [str(column).strip() for column in raw.columns]
    required = {"Title", *PUBLIC_COLUMNS.values()}
    missing = sorted(required - set(raw.columns))
    if missing:
        raise ValueError(f"CypReact output is missing columns: {missing}")
    if raw["Title"].isna().any() or raw["Title"].duplicated().any():
        raise ValueError("CypReact output titles must be present and unique")

    parsed = pd.DataFrame({"Molecule_Name": raw["Title"].astype(str).str.strip()})
    for endpoint, column in PUBLIC_COLUMNS.items():
        calls = raw[column].astype("string").str.strip().str.upper()
        invalid = calls.notna() & ~calls.isin(["R", "N"])
        if invalid.any():
            examples = calls[invalid].head(5).tolist()
            raise ValueError(f"Unexpected {endpoint} CypReact calls: {examples}")
        parsed[f"{endpoint}_cypreact_call"] = calls
        parsed[f"{endpoint}_cypreact_public_score"] = calls.map({"N": 0.0, "R": 1.0})
    return parsed


def _write_input_sdf(frame: pd.DataFrame, path: Path) -> pd.DataFrame:
    writer = Chem.SDWriter(str(path))
    failures = []
    try:
        for row in frame.itertuples(index=False):
            mol = Chem.MolFromSmiles(row.model_smiles)
            if mol is None:
                failures.append(
                    {"Molecule_Name": row.Molecule_Name, "model_smiles": row.model_smiles,
                     "failure": "rdkit_parse_failure"}
                )
                continue
            mol.SetProp("_Name", row.Molecule_Name)
            mol.SetProp("Molecule_Name", row.Molecule_Name)
            mol.SetProp("input_smiles", row.model_smiles)
            writer.write(mol)
    finally:
        writer.close()
    return pd.DataFrame(failures, columns=["Molecule_Name", "model_smiles", "failure"])


def run_inference(
    input_csv: Path,
    bundle_dir: Path,
    output_csv: Path,
    manifest_path: Path,
    failure_csv: Path,
    work_dir: Path,
    id_column: str = "Molecule_Name",
    smiles_column: str = "SMILES",
    limit: int | None = None,
) -> dict[str, Any]:
    input_csv = input_csv.resolve()
    bundle_dir = bundle_dir.resolve()
    jar_path = bundle_dir / "cypreact.jar"
    zip_path = bundle_dir.parent / "CypReactBundle1.2.zip"
    if not jar_path.is_file():
        raise FileNotFoundError(jar_path)
    if not zip_path.is_file() or sha256_file(zip_path) != BUNDLE_ZIP_SHA256:
        raise ValueError("CypReact v1.2 release ZIP is missing or has the wrong SHA-256")

    source = pd.read_csv(input_csv)
    missing = sorted({id_column, smiles_column} - set(source.columns))
    if missing:
        raise ValueError(f"Input is missing required columns: {missing}")
    if source[id_column].isna().any() or source[id_column].duplicated().any():
        raise ValueError(f"{id_column} must be present and unique")
    if source[smiles_column].isna().any():
        raise ValueError(f"{smiles_column} must be present for every row")
    if limit is not None:
        source = source.head(limit).copy()
    source = source[[id_column, smiles_column]].rename(
        columns={id_column: "Molecule_Name", smiles_column: "model_smiles"}
    )

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    failure_csv.parent.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    input_sdf = work_dir / "cypreact_input.sdf"
    raw_output = work_dir / "cypreact_output.csv"
    log_path = work_dir / "cypreact.log"
    failures = _write_input_sdf(source, input_sdf)
    scoreable = source.loc[~source["Molecule_Name"].isin(failures["Molecule_Name"])].copy()
    if scoreable.empty:
        raise ValueError("No scoreable molecules remain after SDF preparation")

    command = [
        "java", "-jar", str(jar_path), str(bundle_dir) + "/", str(input_sdf),
        str(raw_output), "2D6,3A4",
    ]
    with log_path.open("w", encoding="utf-8") as log_handle:
        result = subprocess.run(
            command,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    if result.returncode != 0:
        raise RuntimeError(f"CypReact exited {result.returncode}; see {log_path}")
    if not raw_output.is_file():
        raise FileNotFoundError(raw_output)

    parsed = parse_cypreact_csv(raw_output)
    expected_ids = set(scoreable["Molecule_Name"])
    returned_ids = set(parsed["Molecule_Name"])
    unexpected = sorted(returned_ids - expected_ids)
    if unexpected:
        raise ValueError(f"CypReact returned unexpected molecule IDs: {unexpected[:10]}")
    missing_ids = sorted(expected_ids - returned_ids)
    if missing_ids:
        failures = pd.concat(
            [
                failures,
                pd.DataFrame(
                    {
                        "Molecule_Name": missing_ids,
                        "model_smiles": source.set_index("Molecule_Name").loc[missing_ids, "model_smiles"].tolist(),
                        "failure": "not_returned_by_cypreact",
                    }
                ),
            ],
            ignore_index=True,
        )

    output = source.merge(parsed, on="Molecule_Name", how="left", validate="one_to_one")
    output.to_csv(output_csv, index=False)
    failures.to_csv(failure_csv, index=False)

    version_result = subprocess.run(
        ["java", "-version"], capture_output=True, text=True, check=False
    )
    asset_paths = [
        jar_path,
        bundle_dir / "supportfiles/CYP2D6/supportfile.csv",
        bundle_dir / "supportfiles/CYP3A4/supportfile.csv",
        *sorted((bundle_dir / "supportfiles/CYP2D6/model").glob("*.model")),
        *sorted((bundle_dir / "supportfiles/CYP3A4/model").glob("*.model")),
    ]
    for path in asset_paths:
        if not path.is_file():
            raise FileNotFoundError(path)

    scored = output[[f"{endpoint}_cypreact_public_score" for endpoint in ENDPOINTS]].notna().all(axis=1)
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "model": "CypReact v1.2 public reactant classifier",
        "release": BUNDLE_URL,
        "release_zip_sha256": BUNDLE_ZIP_SHA256,
        "license": "MIT (current CypReact repository); historical bundle retained as released",
        "score_semantics": "faithful public R/N call encoded N=0 and R=1; not a continuous probability",
        "input_csv": str(input_csv),
        "input_sha256": sha256_file(input_csv),
        "input_rows": int(len(source)),
        "scored_rows": int(scored.sum()),
        "failure_rows": int((~scored).sum()),
        "representation": "original supplied SMILES converted to an ID-labelled RDKit SDF",
        "endpoints": ENDPOINTS,
        "command": command,
        "runtime": {
            "python": platform.python_version(),
            "rdkit": Chem.rdBase.rdkitVersion,
            "java": (version_result.stderr or version_result.stdout).strip(),
        },
        "assets": [
            {"path": str(path), "size": path.stat().st_size, "sha256": sha256_file(path)}
            for path in asset_paths
        ],
        "outputs": {
            "predictions": {"path": str(output_csv), "sha256": sha256_file(output_csv)},
            "failures": {"path": str(failure_csv), "sha256": sha256_file(failure_csv)},
            "input_sdf": {"path": str(input_sdf), "sha256": sha256_file(input_sdf)},
            "raw_output": {"path": str(raw_output), "sha256": sha256_file(raw_output)},
            "log": {"path": str(log_path), "sha256": sha256_file(log_path)},
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--bundle-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--failures", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--id-column", default="Molecule_Name")
    parser.add_argument("--smiles-column", default="SMILES")
    parser.add_argument("--limit", type=int)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    manifest = run_inference(
        input_csv=args.input,
        bundle_dir=args.bundle_dir,
        output_csv=args.output,
        manifest_path=args.manifest,
        failure_csv=args.failures,
        work_dir=args.work_dir,
        id_column=args.id_column,
        smiles_column=args.smiles_column,
        limit=args.limit,
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
