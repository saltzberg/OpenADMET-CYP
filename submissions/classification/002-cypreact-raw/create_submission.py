#!/usr/bin/env python3
"""Reproduce the CypReact v1.2 OpenADMET classification artifact."""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[3]
ADAPTER_SOURCE = REPOSITORY / "tdi_public_model_benchmark" / "src"
sys.path.insert(0, str(ADAPTER_SOURCE))

from tdi_public_model_benchmark.cypreact import run_inference  # noqa: E402

DEFAULT_TEST_CSV = REPOSITORY / "data/cyp-challenge-train-test/cyp-challenge-TEST-BLINDED.csv"
DEFAULT_BUNDLE_DIR = REPOSITORY / "tdi_public_model_benchmark/vendor/CypReact/CypReactBundle"
TEST_SHA256 = "a342f8444a8dcb531ca12f3685293f0bd6c36ae9073f491e44a9bc1cc4b741f9"
TIMESTAMP_RE = re.compile(r"^\d{8}T\d{6}Z$")
CANONICAL_COLUMNS = ["SMILES", "Molecule_Name", "CYP2D6_is_TDI", "CYP3A4_is_TDI"]
CALL_COLUMNS = {
    "CYP2D6_is_TDI": "CYP2D6_cypreact_call",
    "CYP3A4_is_TDI": "CYP3A4_cypreact_call",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test-csv", type=Path, default=DEFAULT_TEST_CSV)
    parser.add_argument("--bundle-dir", type=Path, default=DEFAULT_BUNDLE_DIR)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument(
        "--timestamp",
        default=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
        help="UTC timestamp embedded in the canonical filename (YYYYMMDDTHHMMSSZ)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not TIMESTAMP_RE.fullmatch(args.timestamp):
        raise SystemExit("--timestamp must use YYYYMMDDTHHMMSSZ")
    if not args.test_csv.is_file():
        raise SystemExit(f"Blinded test CSV does not exist: {args.test_csv}")
    digest = sha256_file(args.test_csv)
    if digest != TEST_SHA256:
        raise SystemExit(f"Blinded test-set SHA-256 mismatch: {digest}")

    with args.test_csv.open(newline="", encoding="utf-8") as handle:
        source_rows = list(csv.DictReader(handle))
    if len(source_rows) != 750:
        raise SystemExit(f"Expected 750 blinded compounds, found {len(source_rows)}")
    if not source_rows or not {"SMILES", "Molecule_Name"}.issubset(source_rows[0]):
        raise SystemExit("Blinded test CSV must contain SMILES and Molecule_Name")
    source_ids = [row["Molecule_Name"] for row in source_rows]
    if any(not value for value in source_ids) or len(set(source_ids)) != len(source_ids):
        raise SystemExit("Molecule_Name values must be non-empty and unique")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    predictions_path = args.output_dir / "cypreact_predictions.csv"
    failures_path = args.output_dir / "inference_failures.csv"
    manifest_path = args.output_dir / "inference_manifest.json"
    work_dir = args.output_dir / "inference_work"
    manifest = run_inference(
        input_csv=args.test_csv,
        bundle_dir=args.bundle_dir,
        output_csv=predictions_path,
        manifest_path=manifest_path,
        failure_csv=failures_path,
        work_dir=work_dir,
    )
    if manifest["input_rows"] != 750 or manifest["scored_rows"] != 750 or manifest["failure_rows"] != 0:
        raise SystemExit(
            "CypReact did not return complete blind-set coverage; refusing to finalize the artifact"
        )

    with predictions_path.open(newline="", encoding="utf-8") as handle:
        predictions = list(csv.DictReader(handle))
    returned_ids = [row["Molecule_Name"] for row in predictions]
    if returned_ids != source_ids:
        raise SystemExit("CypReact output IDs or ordering differ from the blinded input")

    artifact = args.output_dir / f"openadmet-cyp_classification-sub_dargason_{args.timestamp}.csv"
    with artifact.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CANONICAL_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for source, prediction in zip(source_rows, predictions, strict=True):
            if prediction["model_smiles"] != source["SMILES"]:
                raise SystemExit(f"SMILES changed for {source['Molecule_Name']}")
            output = {"SMILES": source["SMILES"], "Molecule_Name": source["Molecule_Name"]}
            for artifact_column, call_column in CALL_COLUMNS.items():
                call = prediction[call_column].strip().upper()
                if call not in {"N", "R"}:
                    raise SystemExit(f"Missing or invalid {call_column} for {source['Molecule_Name']}: {call!r}")
                output[artifact_column] = call == "R"
            writer.writerow(output)

    print(f"artifact={artifact}")
    print(f"sha256={sha256_file(artifact)}")
    print("rows=750")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
