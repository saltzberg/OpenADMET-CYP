#!/usr/bin/env python3
"""Join row-level QC, emit release provenance, and make the HF dataset card."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


REPO = Path(__file__).resolve().parents[1]
DEFAULT_RELEASE = REPO / "release" / "huggingface_cyp_cofold_v1"
DEFAULT_WORK = REPO / "work" / "cyp_hf_release_v1"
DATASET_VERSION = "1.0.0"
APACHE_LICENSE_SOURCE = REPO / "licenses" / "Apache-2.0.txt"
EXPECTED_METHOD_COUNTS = {
    "experimental": 149,
    "boltz2": 2600,
    "chai1": 2780,
    "esmfold2": 5650,
    "openfold3": 2780,
    "protenix_v1": 2780,
}


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_parquet(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".building")
    frame.to_parquet(temporary, index=False, compression="zstd")
    os.replace(temporary, path)


def atomic_json(data: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".building")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True, default=str) + "\n")
    os.replace(temporary, path)


def coverage_tables(frame: pd.DataFrame, release: Path) -> dict[str, Any]:
    experimental = frame[frame.source_kind.eq("experimental")][
        ["pdb_id", "gene", "ligand_selected_ccd", "ligand_selection_status"]
    ].copy()
    experimental["ligand_bearing"] = experimental.ligand_selected_ccd.notna()
    counts = (
        frame[frame.source_kind.eq("predicted")]
        .groupby(["pdb_id", "method"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    pdb_coverage = experimental.merge(counts, on="pdb_id", how="left", validate="one_to_one").fillna(0)
    methods = [method for method in EXPECTED_METHOD_COUNTS if method != "experimental"]
    for method in methods:
        if method not in pdb_coverage:
            pdb_coverage[method] = 0
        pdb_coverage[f"{method}_covered"] = pdb_coverage[method].gt(0)
    pdb_coverage["covered_by_any_prediction_method"] = pdb_coverage[methods].sum(axis=1).gt(0)
    pdb_coverage["covered_by_all_five_methods"] = pdb_coverage[[f"{method}_covered" for method in methods]].all(axis=1)
    pdb_coverage.to_csv(release / "provenance" / "pdb_coverage.csv", index=False)

    rows: list[dict[str, Any]] = []
    ligand_pdbs = set(experimental.loc[experimental.ligand_bearing, "pdb_id"])
    for method in methods:
        subset = frame[frame.method.eq(method)]
        covered = set(subset.pdb_id)
        for gene in sorted(experimental.gene.unique()):
            gene_universe = set(experimental.loc[(experimental.gene == gene) & experimental.ligand_bearing, "pdb_id"])
            rows.append(
                {
                    "method": method,
                    "gene": gene,
                    "structure_rows": int(subset.gene.eq(gene).sum()),
                    "covered_ligand_pdb_ids": len(covered & gene_universe),
                    "ligand_pdb_ids_in_gene": len(gene_universe),
                    "coverage_fraction": len(covered & gene_universe) / len(gene_universe) if gene_universe else None,
                }
            )
        rows.append(
            {
                "method": method,
                "gene": "ALL",
                "structure_rows": len(subset),
                "covered_ligand_pdb_ids": len(covered & ligand_pdbs),
                "ligand_pdb_ids_in_gene": len(ligand_pdbs),
                "coverage_fraction": len(covered & ligand_pdbs) / len(ligand_pdbs),
            }
        )
    method_coverage = pd.DataFrame(rows)
    method_coverage.to_csv(release / "provenance" / "method_coverage.csv", index=False)
    overall = method_coverage[method_coverage.gene.eq("ALL")].set_index("method")
    return {
        method: {
            "structure_rows": int(overall.at[method, "structure_rows"]),
            "covered_ligand_pdb_ids": int(overall.at[method, "covered_ligand_pdb_ids"]),
            "ligand_pdb_ids": int(overall.at[method, "ligand_pdb_ids_in_gene"]),
            "coverage_fraction": float(overall.at[method, "coverage_fraction"]),
        }
        for method in methods
    }


def excluded_attempts(release: Path) -> pd.DataFrame:
    source = REPO / "analysis" / "cyp_cofold_postprocess_20260804" / "failure_ledger.csv"
    failures = pd.read_csv(source, low_memory=False)
    excluded = failures[failures.stage.eq("alignment")].copy()
    if len(excluded) != 40:
        raise ValueError(f"unexpected_alignment_exclusion_count:{len(excluded)}")
    excluded["source_filename"] = excluded.model_path.map(lambda value: Path(str(value)).name)
    excluded = excluded.drop(columns=["model_path"])
    atomic_parquet(excluded, release / "provenance" / "excluded_attempts.parquet")
    return excluded


def field_groups(columns: list[str]) -> dict[str, list[str]]:
    prefixes = {
        "identity_and_files": ["structure_id", "source_kind", "method", "pdb_id", "gene", "accession", "seed", "sample_index", "coordinate_path", "coordinate_sha256"],
        "ligand_and_role": [column for column in columns if column.startswith("ligand_") or column in {"ground_truth_target_ligand_ccd", "ost_target_role_valid"}],
        "alignment": [column for column in columns if column.startswith("alignment_")],
        "cofold_confidence_and_ost": [column for column in columns if column.startswith("ost_") or column in {"pl_iptm", "native_ligand_iptm", "native_iptm"}],
        "prolif": [column for column in columns if column.startswith("prolif_")],
        "posebusters": [column for column in columns if column.startswith("posebusters_")],
        "training_windows": [column for column in columns if "training_" in column or column in {"deposit_date", "initial_release_date"}],
        "experimental_structure": [column for column in columns if column.startswith("heme_") or column in {"resolution", "spacegroup", "unit_cell", "axial_cys_sg_to_1tqn_A", "axial_cys_resnum", "fe_cys_A"}],
        "qc_execution": [column for column in columns if column.startswith("qc_") or column in {"component_extraction_status", "coordinate_heavy_atoms", "protein_heavy_atoms"}],
    }
    return {group: sorted(set(values) & set(columns)) for group, values in prefixes.items()}


def build_readme(summary: dict[str, Any]) -> str:
    method_lines = "\n".join(
        f"| {method} | {data['structure_rows']:,} | {data['covered_ligand_pdb_ids']}/{data['ligand_pdb_ids']} | {data['coverage_fraction']:.1%} |"
        for method, data in summary["method_coverage"].items()
    )
    posebusters_pass = summary["posebusters_overall_pass_counts"].get("True", 0)
    posebusters_fail = summary["posebusters_overall_pass_counts"].get("False", 0)
    return f"""---
pretty_name: OpenADMET CYP PDB cofolding structures v1
license: apache-2.0
task_categories:
- feature-extraction
- table-question-answering
tags:
- structural-biology
- cytochrome-p450
- protein-ligand
- cofolding
- prolif
- posebusters
configs:
- config_name: default
  data_files:
  - split: all
    path: data/structures.parquet
---

# OpenADMET CYP PDB cofolding structures v1

This is a coordinate-complete, one-row-per-structure release for the current UniProt PDB cross-references of human CYP1A2, CYP2C9, CYP2D6, and CYP3A4. It contains **{summary['rows']:,} rows**: {summary['experimental_rows']:,} aligned experimental ground truths and {summary['predicted_rows']:,} valid aligned cofolded structures.

Every coordinate file is gzip-compressed mmCIF and every row points to exactly one relative `coordinate_path`. All structures use the documented `cyp3a4_1tqn_conserved_scaffold_ca_v1` frame: a proper Kabsch rotation over 187 conserved CYP scaffold C-alpha anchors, with at least 180 anchors and every block required. The same rigid transform is applied to every atom; ligand, heme, and pocket atoms are never fitted independently.

## Method coverage

Coverage below is against the 141 ligand-bearing PDB IDs. The other eight experimental structures are apo/heme-only and remain in the dataset with ligand QC marked `not_applicable_no_selected_ligand`.

| Method | Rows | Ligand PDB IDs | Coverage |
|---|---:|---:|---:|
{method_lines}

The union of generated methods covers 139/141 ligand PDB IDs; 4WNU and 4WNV have no valid generated ligand complex. ESMFold2 was run only for CYP3A4. Protenix v2 seed outputs and raw apo cofold attempts are outside this v1 core because they have not passed the same alignment/QC protocol.

## Row-level derived data

`data/structures.parquet` contains identifiers and hashes, ligand chemistry and role fields, alignment metrics, native method confidence, OpenStructure lDDT-PLI where assignment succeeded, training-window labels, ProLIF interaction counts and contact JSON, and the complete PoseBusters `dock` report (all raw report fields are prefixed `posebusters_`). Null lDDT-PLI values with `ost_status` are retained; they are not silently filtered.

ProLIF and PoseBusters use the selected CCD SMILES as authoritative ligand topology, element/graph-isomorphism atom mapping, and the coordinate pose without numerical movement. ProLIF uses a sequence-template RDKit protein graph and canonical UniProt residue numbering. PoseBusters uses the CYP protein chain only as `mol_cond`; heme, waters, additives, and other cofactors are deliberately excluded from the receptor condition. Ir and Ru use explicit 2.00 and 2.05 A ProLIF van der Waals radii. Legacy PDB `CONECT` records are used when distance inference cannot recover a deposited organometallic graph.

The automated 6CSB nearest-heme choice was role-curated: the experimental ground truth uses complete RTZ chain A residue 602 rather than a partial HEGA-10 (`2CV`) detergent fragment. Existing generated 2CV structures are retained because they are real campaign outputs, but they carry `ligand_role=crystallization_additive_automated_selection` and `ost_target_role_valid=false`.

## Quality states

ProLIF and PoseBusters executed successfully for all **{summary['prolif_status_counts']['ok']:,}** selected-ligand rows; the eight apo/heme-only rows are explicitly not applicable. PoseBusters reports an all-check scientific pass for {posebusters_pass:,} rows and at least one failed check for {posebusters_fail:,} rows. These are scientific results, not execution statuses.

- `qc_status=complete`: ProLIF and PoseBusters both executed. A complete row can still fail one or more PoseBusters scientific checks; inspect `posebusters_overall_pass` and `posebusters_failed_checks_json`.
- `qc_status=not_applicable_no_selected_ligand`: experimental apo/heme-only structure.
- `ost_status` records lDDT-PLI assignment or the reason no value is available.
- No invalid/silent ligand-omission coordinate is included. The 40 excluded Chai-1 attempts are listed in `provenance/excluded_attempts.parquet`.

## Files and verification

- `data/structures.parquet`: canonical one-row-per-structure table.
- `structures/experimental/`: 149 rigidly aligned PDB ground truths.
- `structures/predicted/`: 16,590 verified aligned predictions.
- `provenance/`: coverage, exact transforms, anchor map, protocols, excluded attempts, field groups, and checksums.
- `scripts/verify_release.py`: full hash, schema, row/path, and coordinate-parse verification.

Run `python scripts/verify_release.py .` from the dataset root. `provenance/SHA256SUMS` covers every published file except itself and the post-verification `provenance/verification.json` report.

## Scope and limitations

Training-window fields classify PDB release dates relative to documented method cutoffs; they do **not** prove training-set membership. Prediction method names are campaign labels; exact predictor commit and container digests were not consistently captured in the source manifests. Ligand role selection is mostly an automated nearest-heme heuristic, with the explicit 6CSB correction above. PoseBusters internal-energy checks can be unavailable/fail for organometallic ligands because UFF lacks parameters; this is preserved as a failed scientific check, not treated as a pipeline failure. OpenADMET-generated structures, derived tables, documentation, and release scripts are provided under Apache-2.0. Experimental PDB coordinate data retain their upstream CC0 status; see `LICENSES.md`.

## Citation

Use `CITATION.cff` for this dataset snapshot and cite the original PDB entries and structure-prediction methods appropriate to your analysis.
"""


def checksums(release: Path) -> None:
    target = release / "provenance" / "SHA256SUMS"
    files = sorted(path for path in release.rglob("*") if path.is_file() and path != target and path != release / "provenance" / "verification.json")
    temporary = target.with_suffix(".building")
    with temporary.open("w") as handle:
        for index, path in enumerate(files, start=1):
            handle.write(f"{sha256(path)}  {path.relative_to(release).as_posix()}\n")
            if index % 2500 == 0:
                print(json.dumps({"event": "checksum_progress", "done": index, "total": len(files)}), flush=True)
    os.replace(temporary, target)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, default=DEFAULT_WORK / "base_structures.parquet")
    parser.add_argument("--qc", type=Path, default=DEFAULT_WORK / "qc" / "structure_qc.parquet")
    parser.add_argument("--release", type=Path, default=DEFAULT_RELEASE)
    args = parser.parse_args()
    release = args.release.resolve()
    base = pd.read_parquet(args.base)
    qc = pd.read_parquet(args.qc)
    if len(base) != 16739 or len(qc) != 16739:
        raise ValueError(f"input_universe_gate_failed:base={len(base)}:qc={len(qc)}")
    if qc.structure_id.nunique() != 16739:
        raise ValueError("qc_structure_id_uniqueness_gate_failed")
    duplicate = set(base.columns) & (set(qc.columns) - {"structure_id"})
    if duplicate:
        raise ValueError(f"unexpected_join_column_overlap:{sorted(duplicate)}")
    frame = base.merge(qc, on="structure_id", validate="one_to_one")
    frame.insert(0, "dataset_version", DATASET_VERSION)
    frame.insert(1, "record_index", range(len(frame)))
    if frame.method.value_counts().to_dict() != EXPECTED_METHOD_COUNTS:
        raise ValueError(f"method_count_gate_failed:{frame.method.value_counts().to_dict()}")
    if not frame.coordinate_path.map(lambda value: (release / str(value)).is_file()).all():
        raise ValueError("missing_coordinate_file_gate_failed")
    if not frame.alignment_status.eq("PASS").all():
        raise ValueError("alignment_status_gate_failed")
    if int(frame.qc_status.eq("not_applicable_no_selected_ligand").sum()) != 8:
        raise ValueError("apo_qc_status_gate_failed")

    output = release / "data" / "structures.parquet"
    atomic_parquet(frame, output)
    method_coverage = coverage_tables(frame, release)
    excluded = excluded_attempts(release)
    package_versions = {
        "python": __import__("sys").version.split()[0],
        "pandas": pd.__version__,
        "pyarrow": pa.__version__,
        "gemmi": __import__("gemmi").__version__,
        "rdkit": __import__("rdkit").__version__,
        "prolif": __import__("prolif").__version__,
        "posebusters": __import__("posebusters").__version__,
    }
    summary = {
        "dataset_version": DATASET_VERSION,
        "created_at": now(),
        "rows": len(frame),
        "columns": len(frame.columns),
        "experimental_rows": int(frame.source_kind.eq("experimental").sum()),
        "predicted_rows": int(frame.source_kind.eq("predicted").sum()),
        "unique_pdb_ids": int(frame.pdb_id.nunique()),
        "ligand_bearing_ground_truths": int(frame.source_kind.eq("experimental").mul(frame.ligand_selected_ccd.notna()).sum()),
        "apo_or_heme_only_ground_truths": int(frame.source_kind.eq("experimental").mul(frame.ligand_selected_ccd.isna()).sum()),
        "method_counts": frame.method.value_counts().sort_index().to_dict(),
        "qc_status_counts": frame.qc_status.value_counts(dropna=False).to_dict(),
        "prolif_status_counts": frame.prolif_status.value_counts(dropna=False).to_dict(),
        "posebusters_status_counts": frame.posebusters_status.value_counts(dropna=False).to_dict(),
        "posebusters_overall_pass_counts": {str(key): int(value) for key, value in frame.posebusters_overall_pass.value_counts(dropna=False).items()},
        "ost_status_counts": frame.ost_status.value_counts(dropna=False).to_dict(),
        "method_coverage": method_coverage,
        "excluded_alignment_attempts": len(excluded),
        "alignment_invariant": "cyp3a4_1tqn_conserved_scaffold_ca_v1",
        "qc_protocol": frame.qc_protocol.unique().tolist(),
        "package_versions": package_versions,
        "parquet_sha256": sha256(output),
        "parquet_bytes": output.stat().st_size,
    }
    atomic_json(summary, release / "provenance" / "release_summary.json")
    atomic_json(
        {
            "dataset_version": DATASET_VERSION,
            "coordinate_policy": "proper rigid Kabsch alignment; one transform applied to every atom; no internal coordinate changes",
            "ligand_topology_policy": "selected CCD SMILES authoritative; heavy-atom element and graph-isomorphism mapping; legacy PDB CONECT fallback for deposited organometallic graphs",
            "prolif_receptor": "selected CYP protein chain only; heme/water/additives/cofactors excluded",
            "posebusters_receptor": "selected CYP protein chain only; PoseBusters config=dock full_report=True",
            "apo_policy": "included experimental coordinates; ligand QC explicitly not applicable",
            "role_curation": {"6CSB": "experimental target changed from partial HEGA-10 2CV additive to complete RTZ chain A residue 602; generated 2CV outputs retained and flagged"},
            "package_versions": package_versions,
        },
        release / "provenance" / "protocol.json",
    )
    source_protocol = json.loads((REPO / "analysis" / "cyp_cofold_postprocess_20260804" / "protocol.json").read_text())
    source_protocol["coordinate_alignment"]["aligned_coordinate_root"] = "structures/predicted"
    atomic_json(source_protocol, release / "provenance" / "source_postprocess_protocol.json")
    atomic_json(field_groups(list(frame.columns)), release / "provenance" / "field_groups.json")
    (release / ".gitattributes").write_text("*.cif.gz filter=lfs diff=lfs merge=lfs -text\n*.parquet filter=lfs diff=lfs merge=lfs -text\n")
    (release / "README.md").write_text(build_readme(summary))
    (release / "LICENSES.md").write_text(
        """# Component licensing\n\nThe Hugging Face dataset card uses `license: apache-2.0`. The full license is in `LICENSE`.\n\n## Apache-2.0 components\n\nCopyright 2026 OpenADMET contributors. OpenADMET-generated coordinate structures, derived Parquet and provenance data, documentation, and release scripts are licensed under the Apache License, Version 2.0, to the extent the contributors hold rights in those materials. The structure-prediction software and model weights themselves are not redistributed in this dataset.\n\n## Experimental PDB coordinates: CC0-1.0\n\nThe [RCSB PDB usage policy](https://www.rcsb.org/pages/usage-policy) states that data files in the PDB archive and data returned by RCSB programmatic APIs are available under the CC0 1.0 Universal Public Domain Dedication. RCSB encourages attribution of the original structure authors. Cite the PDB IDs and corresponding structure publications used in an analysis. Policy checked 2026-08-09.\n\nThe CC0 designation for those upstream files is preserved; Apache-2.0 does not restrict the underlying CC0 material.\n"""
    )
    shutil.copy2(APACHE_LICENSE_SOURCE, release / "LICENSE")
    (release / "NOTICE").write_text(
        "OpenADMET CYP PDB cofolding structures\nCopyright 2026 OpenADMET contributors\n\nExperimental PDB coordinate files originate from the Protein Data Bank and are identified by row-level PDB IDs. Cite the original structure authors and corresponding publications where practical.\n"
    )
    (release / "CITATION.cff").write_text(
        f"""cff-version: 1.2.0\nmessage: \"If you use this dataset, cite this dataset snapshot, the relevant PDB entries, and the prediction methods.\"\ntitle: \"OpenADMET CYP PDB cofolding structures\"\nversion: \"{DATASET_VERSION}\"\ndate-released: \"{now()[:10]}\"\nauthors:\n  - name: \"OpenADMET CYP dataset contributors\"\ntype: dataset\n"""
    )
    release_scripts = release / "scripts"
    release_scripts.mkdir(parents=True, exist_ok=True)
    shutil.copy2(REPO / "scripts" / "verify_hf_release.py", release_scripts / "verify_release.py")
    pipeline_scripts = release / "provenance" / "pipeline_scripts"
    pipeline_scripts.mkdir(parents=True, exist_ok=True)
    for name in ["build_hf_release.py", "compute_release_qc.py", "finalize_hf_release.py"]:
        source = REPO / "scripts" / name
        if source.is_file():
            shutil.copy2(source, pipeline_scripts / name)
    checksums(release)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
