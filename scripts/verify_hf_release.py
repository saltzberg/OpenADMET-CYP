#!/usr/bin/env python3
"""Verify the published CYP structure release from its dataset root."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import gemmi
import pandas as pd
import pyarrow.parquet as pq


EXPECTED_METHOD_COUNTS = {
    "experimental": 149,
    "boltz2": 2600,
    "chai1": 2780,
    "esmfold2": 5650,
    "openfold3": 2780,
    "protenix_v1": 2780,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def verify_checksum_manifest(root: Path) -> tuple[int, list[str]]:
    manifest = root / "provenance" / "SHA256SUMS"
    require(manifest.is_file(), "missing_SHA256SUMS")
    failures: list[str] = []
    count = 0
    for line in manifest.read_text().splitlines():
        expected, relative = line.split("  ", 1)
        path = root / relative
        count += 1
        if not path.is_file():
            failures.append(f"missing:{relative}")
        elif sha256(path) != expected:
            failures.append(f"hash:{relative}")
    return count, failures


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path, nargs="?", default=Path("."))
    parser.add_argument("--skip-coordinate-parse", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    parquet_path = root / "data" / "structures.parquet"
    require(parquet_path.is_file(), "missing_structures_parquet")
    frame = pd.read_parquet(parquet_path)
    require(len(frame) == 16739, f"row_count:{len(frame)}")
    require(frame.structure_id.nunique() == 16739, "structure_id_uniqueness")
    require(frame.coordinate_path.nunique() == 16739, "coordinate_path_uniqueness")
    require(frame.method.value_counts().to_dict() == EXPECTED_METHOD_COUNTS, "method_counts")
    require(frame.pdb_id.nunique() == 149, "pdb_id_count")
    require(frame.source_kind.eq("experimental").sum() == 149, "experimental_count")
    require(frame.source_kind.eq("predicted").sum() == 16590, "predicted_count")
    require(frame.alignment_status.eq("PASS").all(), "alignment_status")
    require(frame.alignment_invariant.eq("cyp3a4_1tqn_conserved_scaffold_ca_v1").all(), "alignment_invariant")
    require(frame.qc_protocol.nunique() == 1, "qc_protocol_uniqueness")
    require(frame.qc_status.eq("not_applicable_no_selected_ligand").sum() == 8, "apo_qc_count")
    require(frame.loc[frame.qc_status.ne("not_applicable_no_selected_ligand"), "prolif_status"].eq("ok").all(), "prolif_execution_coverage")
    require(frame.loc[frame.qc_status.ne("not_applicable_no_selected_ligand"), "posebusters_status"].eq("ok").all(), "posebusters_execution_coverage")
    require(not any(str(value).startswith("/") for value in frame.coordinate_path), "absolute_coordinate_path")

    missing: list[str] = []
    coordinate_hash_failures: list[str] = []
    parse_failures: list[str] = []
    for index, row in enumerate(frame[["coordinate_path", "coordinate_sha256"]].itertuples(index=False), start=1):
        path = root / row.coordinate_path
        if not path.is_file():
            missing.append(row.coordinate_path)
            continue
        if sha256(path) != row.coordinate_sha256:
            coordinate_hash_failures.append(row.coordinate_path)
        if not args.skip_coordinate_parse:
            try:
                structure = gemmi.read_structure(str(path))
                if not len(structure) or not len(structure[0]):
                    raise ValueError("empty_structure")
            except Exception as exc:
                parse_failures.append(f"{row.coordinate_path}:{type(exc).__name__}:{exc}")
        if index % 2500 == 0:
            print(json.dumps({"event": "coordinate_verify_progress", "done": index, "total": len(frame)}), flush=True)
    require(not missing, f"missing_coordinates:{missing[:10]}")
    require(not coordinate_hash_failures, f"coordinate_hash_failures:{coordinate_hash_failures[:10]}")
    require(not parse_failures, f"coordinate_parse_failures:{parse_failures[:10]}")

    checksum_count, checksum_failures = verify_checksum_manifest(root)
    require(not checksum_failures, f"checksum_manifest_failures:{checksum_failures[:10]}")
    metadata = pq.read_metadata(parquet_path)
    result = {
        "status": "PASS",
        "rows": len(frame),
        "columns": len(frame.columns),
        "parquet_row_groups": metadata.num_row_groups,
        "coordinate_files": len(frame),
        "coordinate_parse_checked": not args.skip_coordinate_parse,
        "checksum_manifest_files": checksum_count,
        "method_counts": frame.method.value_counts().sort_index().to_dict(),
        "qc_status_counts": frame.qc_status.value_counts(dropna=False).to_dict(),
        "posebusters_overall_pass_counts": {str(key): int(value) for key, value in frame.posebusters_overall_pass.value_counts(dropna=False).items()},
    }
    output = root / "provenance" / "verification.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
