"""Prepare an auditable long-form TDI label registry."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.warning")

ENDPOINT_COLUMNS = {
    "CYP2D6": "CYP2D6_is_TDI",
    "CYP3A4": "CYP3A4_is_TDI",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _parse_label(value: object) -> int | None:
    if pd.isna(value):
        return None
    if isinstance(value, bool):
        return int(value)
    text = str(value).strip().lower()
    if text in {"true", "1", "yes"}:
        return 1
    if text in {"false", "0", "no"}:
        return 0
    raise ValueError(f"Unsupported TDI label value: {value!r}")


def scaffold_for_smiles(smiles: str, molecule_name: str) -> tuple[str, str, bool]:
    """Return canonical SMILES, bootstrap scaffold group, and parse status."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles, f"INVALID::{molecule_name}", False
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    if scaffold:
        return canonical, f"MURCKO::{scaffold}", True
    # Empty Murcko scaffolds are chemically heterogeneous. Use molecular-identity
    # singleton clusters instead of treating every acyclic molecule as one family.
    key = hashlib.sha1(canonical.encode("utf-8")).hexdigest()[:16]
    return canonical, f"ACYCLIC::{key}", True


def prepare_labels(source_csv: Path, output_dir: Path) -> dict[str, object]:
    source_csv = source_csv.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(source_csv)

    required = {"Molecule_Name", "SMILES", *ENDPOINT_COLUMNS.values()}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"TDI source is missing required columns: {missing}")
    if frame["Molecule_Name"].isna().any() or frame["Molecule_Name"].duplicated().any():
        raise ValueError("Molecule_Name must be present and unique in the TDI source")
    if frame["SMILES"].isna().any():
        raise ValueError("SMILES must be present for every TDI source row")

    structure_rows = []
    for row in frame[["Molecule_Name", "SMILES"]].itertuples(index=False):
        canonical, scaffold, valid = scaffold_for_smiles(row.SMILES, row.Molecule_Name)
        structure_rows.append(
            {
                "Molecule_Name": row.Molecule_Name,
                "canonical_smiles": canonical,
                "scaffold_id": scaffold,
                "smiles_valid": valid,
            }
        )
    structures = pd.DataFrame(structure_rows)
    frame = frame.merge(structures, on="Molecule_Name", validate="one_to_one")

    long_parts = []
    for endpoint, column in ENDPOINT_COLUMNS.items():
        labels = frame[column].map(_parse_label)
        mask = labels.notna()
        part = frame.loc[
            mask,
            ["Molecule_Name", "SMILES", "canonical_smiles", "scaffold_id", "smiles_valid"],
        ].copy()
        part["endpoint"] = endpoint
        part["is_tdi"] = labels.loc[mask].astype("int8")
        long_parts.append(part)

    labels_long = pd.concat(long_parts, ignore_index=True)
    labels_long = labels_long.sort_values(["endpoint", "Molecule_Name"]).reset_index(drop=True)
    if labels_long.duplicated(["Molecule_Name", "endpoint"]).any():
        raise ValueError("Duplicate molecule/endpoint labels were generated")

    model_input = (
        labels_long[["Molecule_Name", "SMILES", "canonical_smiles", "scaffold_id", "smiles_valid"]]
        .drop_duplicates("Molecule_Name")
        .sort_values("Molecule_Name")
        .reset_index(drop=True)
    )
    conflicts = labels_long.groupby("Molecule_Name")["SMILES"].nunique()
    if (conflicts > 1).any():
        raise ValueError("A molecule ID maps to conflicting SMILES values")

    labels_path = output_dir / "tdi_labels.csv"
    inputs_path = output_dir / "model_input.csv"
    manifest_path = output_dir / "label_manifest.json"
    labels_long.to_csv(labels_path, index=False)
    model_input.to_csv(inputs_path, index=False)

    endpoint_summary = {}
    for endpoint, group in labels_long.groupby("endpoint", sort=True):
        endpoint_summary[endpoint] = {
            "n_labels": int(len(group)),
            "n_positive": int(group["is_tdi"].sum()),
            "n_negative": int((group["is_tdi"] == 0).sum()),
            "positive_prevalence": float(group["is_tdi"].mean()),
            "n_scaffolds": int(group["scaffold_id"].nunique()),
        }
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_csv": str(source_csv),
        "source_sha256": sha256_file(source_csv),
        "label_endpoints": list(ENDPOINT_COLUMNS),
        "endpoint_summary": endpoint_summary,
        "n_unique_model_inputs": int(len(model_input)),
        "n_invalid_smiles": int((~model_input["smiles_valid"]).sum()),
        "invalid_molecule_names": model_input.loc[
            ~model_input["smiles_valid"], "Molecule_Name"
        ].tolist(),
        "outputs": {
            "tdi_labels.csv": sha256_file(labels_path),
            "model_input.csv": sha256_file(inputs_path),
        },
        "scaffold_policy": {
            "rings": "Bemis-Murcko scaffold without chirality",
            "acyclic": "singleton group by canonical molecular identity",
            "invalid_smiles": "singleton group by Molecule_Name",
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest
