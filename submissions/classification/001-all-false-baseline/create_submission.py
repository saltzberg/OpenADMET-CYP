#!/usr/bin/env python3
"""Create the all-False OpenADMET CYP classification submission."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import re
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[3]
DEFAULT_TEST_CSV = (
    REPOSITORY
    / "data"
    / "cyp-challenge-train-test"
    / "cyp-challenge-TEST-BLINDED.csv"
)
TEST_URL = "https://huggingface.co/datasets/openadmet/cyp-challenge-train-test/resolve/main/cyp-challenge-TEST-BLINDED.csv"
TEST_SHA256 = "a342f8444a8dcb531ca12f3685293f0bd6c36ae9073f491e44a9bc1cc4b741f9"
ENDPOINTS = ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]
TIMESTAMP_RE = re.compile(r"^\d{8}T\d{6}Z$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--test-csv", type=Path, default=DEFAULT_TEST_CSV)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument(
        "--timestamp",
        default=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
        help="UTC timestamp used in the filename (YYYYMMDDTHHMMSSZ)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not TIMESTAMP_RE.fullmatch(args.timestamp):
        raise SystemExit("--timestamp must use YYYYMMDDTHHMMSSZ")

    if args.test_csv.is_file():
        test_bytes = args.test_csv.read_bytes()
    elif args.test_csv == DEFAULT_TEST_CSV:
        with urllib.request.urlopen(TEST_URL, timeout=60) as response:
            test_bytes = response.read()
    else:
        raise SystemExit(f"Test CSV does not exist: {args.test_csv}")
    digest = hashlib.sha256(test_bytes).hexdigest()
    if digest != TEST_SHA256:
        raise SystemExit(f"Blinded test-set SHA-256 mismatch: {digest}")
    rows = list(csv.DictReader(io.StringIO(test_bytes.decode("utf-8"))))

    if len(rows) != 750:
        raise SystemExit(f"Expected 750 blinded compounds, found {len(rows)}")
    required = {"SMILES", "Molecule_Name"}
    if not rows or not required.issubset(rows[0]):
        raise SystemExit("Test CSV must contain SMILES and Molecule_Name")
    molecule_names = [row["Molecule_Name"] for row in rows]
    if len(set(molecule_names)) != len(molecule_names):
        raise SystemExit("Molecule_Name values must be unique")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / (
        f"openadmet-cyp_classification-sub_dargason_{args.timestamp}.csv"
    )
    columns = ["SMILES", "Molecule_Name", *ENDPOINTS]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "SMILES": row["SMILES"],
                    "Molecule_Name": row["Molecule_Name"],
                    **{endpoint: False for endpoint in ENDPOINTS},
                }
            )

    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
