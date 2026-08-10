#!/usr/bin/env python3
"""Build the aligned coordinate universe and base one-row-per-structure table.

This script intentionally stops before ProLIF/PoseBusters.  The expensive QC is
resumable and is implemented in ``compute_release_qc.py``; the two tables are
joined by ``finalize_hf_release.py``.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import gemmi
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[1]
DEFAULT_RELEASE = REPO / "release" / "huggingface_cyp_cofold_v1"
DEFAULT_WORK = REPO / "work" / "cyp_hf_release_v1"
ALIGNMENT_INVARIANT = "cyp3a4_1tqn_conserved_scaffold_ca_v1"
INTERVAL = {
    "CYP1A2": (27, 516),
    "CYP2C9": (30, 490),
    "CYP2D6": (34, 497),
    "CYP3A4": (30, 496),
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


def position(atom: gemmi.Atom) -> np.ndarray:
    return np.array([atom.pos.x, atom.pos.y, atom.pos.z], dtype=float)


def residue_atom(residue: gemmi.Residue, name: str) -> gemmi.Atom | None:
    for atom in residue:
        if atom.name.strip() == name:
            return atom
    return None


def chain_ca(structure: gemmi.Structure, chain_name: str = "A") -> list[tuple[int, np.ndarray]]:
    chain = structure[0].find_chain(chain_name)
    if not chain:
        raise ValueError(f"missing_chain:{chain_name}")
    output: list[tuple[int, np.ndarray]] = []
    for residue in chain:
        atom = residue_atom(residue, "CA")
        if atom is not None:
            output.append((residue.seqid.num, position(atom)))
    return output


def kabsch(moving: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    moving_center = moving.mean(axis=0)
    reference_center = reference.mean(axis=0)
    u, _singular, vt = np.linalg.svd((moving - moving_center).T @ (reference - reference_center))
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1] *= -1
        rotation = vt.T @ u.T
    translation = reference_center - rotation @ moving_center
    transformed = moving @ rotation.T + translation
    rmsd = float(np.sqrt(np.mean(np.sum((transformed - reference) ** 2, axis=1))))
    return rotation, translation, rmsd


class CanonicalAligner:
    def __init__(self, reference: Path, anchor_map: Path) -> None:
        self.reference_path = reference
        self.anchor_map_path = anchor_map
        self.reference_sha256 = sha256(reference)
        self.anchor_map_sha256 = sha256(anchor_map)
        self.reference_ca = dict(chain_ca(gemmi.read_structure(str(reference))))
        with anchor_map.open(newline="") as handle:
            self.rows = list(csv.DictReader(handle))

    def align(self, source: Path, gene: str, output: Path) -> dict[str, Any]:
        structure = gemmi.read_structure(str(source))
        residues = chain_ca(structure)
        start, end = INTERVAL[gene]
        numbering_mode = (
            "core-index"
            if len(residues) == end - start + 1 and [item[0] for item in residues[:3]] == [1, 2, 3]
            else "uniprot"
        )
        if numbering_mode == "core-index":
            coordinates = {start + index: xyz for index, (_number, xyz) in enumerate(residues)}
        else:
            coordinates = dict(residues)

        pairs: list[tuple[int, int]] = []
        found_rows: list[dict[str, str]] = []
        for row in self.rows:
            reference_number = int(row["reference_resnum"])
            target_number = int(row[f"{gene}_resnum"])
            if reference_number in self.reference_ca and target_number in coordinates:
                pairs.append((reference_number, target_number))
                found_rows.append(row)
        blocks = sorted({row["block"] for row in self.rows})
        coverage_by_block = {block: sum(row["block"] == block for row in found_rows) for block in blocks}
        if len(pairs) < 180 or any(count == 0 for count in coverage_by_block.values()):
            raise ValueError(
                f"anchor_coverage:{len(pairs)}/{len(self.rows)}:blocks={coverage_by_block}"
            )
        moving = np.array([coordinates[target] for _reference, target in pairs])
        reference = np.array([self.reference_ca[ref] for ref, _target in pairs])
        rotation, translation, rmsd = kabsch(moving, reference)

        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        transformed = rotation @ position(atom) + translation
                        atom.pos = gemmi.Position(*(float(value) for value in transformed))

        output.parent.mkdir(parents=True, exist_ok=True)
        coordinate_bytes = structure.make_mmcif_document().as_string().encode("utf-8")
        temporary = output.with_suffix(output.suffix + ".building")
        with temporary.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0) as compressed:
                compressed.write(coordinate_bytes)
        os.replace(temporary, output)
        return {
            "invariant": ALIGNMENT_INVARIANT,
            "source_path": str(source),
            "source_sha256": sha256(source),
            "gene": gene,
            "numbering_mode": numbering_mode,
            "reference_path": str(self.reference_path),
            "reference_sha256": self.reference_sha256,
            "anchor_map_path": str(self.anchor_map_path),
            "anchor_map_sha256": self.anchor_map_sha256,
            "anchor_pairs": len(pairs),
            "anchor_pairs_expected": len(self.rows),
            "anchor_coverage_by_block": coverage_by_block,
            "anchor_rmsd_A": rmsd,
            "rotation_3x3": rotation.tolist(),
            "translation_A": translation.tolist(),
            "coordinate_action": "x_canonical = rotation_3x3 @ x_source + translation_A; one rigid transform applied to every atom",
            "aligned_relative_path": output.as_posix(),
            "aligned_sha256": sha256(output),
        }


def consistent_truth_selection(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, low_memory=False)
    selected = frame[frame.selected_ccd.notna()].copy()
    columns = [
        "pdb_id",
        "selected_ccd",
        "truth_ligand_chain_selected",
        "truth_ligand_resnum_selected",
        "truth_ligand_atoms",
        "status",
        "ost_lddt_pli",
        "ost_lddt_pli_n_contacts",
        "error",
    ]
    selected = selected[columns]
    for pdb_id, group in selected.groupby("pdb_id"):
        identities = group[
            ["selected_ccd", "truth_ligand_chain_selected", "truth_ligand_resnum_selected"]
        ].drop_duplicates()
        if len(identities) != 1:
            raise ValueError(f"inconsistent_truth_selection:{pdb_id}:{identities.to_dict('records')}")
    return selected.sort_values(["pdb_id", "status"]).drop_duplicates("pdb_id", keep="first")


def release_training_columns() -> list[str]:
    return [
        "deposit_date",
        "initial_release_date",
        "chai1_training_cutoff",
        "chai1_training_window_label",
        "chai1_outside_training_by_release_date",
        "chai1_cutoff_confidence",
        "boltz2_training_cutoff",
        "boltz2_training_window_label",
        "boltz2_outside_training_by_release_date",
        "boltz2_cutoff_confidence",
        "protenix_training_cutoff",
        "protenix_training_window_label",
        "protenix_outside_training_by_release_date",
        "protenix_cutoff_confidence",
        "openfold3_preview_training_cutoff",
        "openfold3_preview_training_window_label",
        "openfold3_preview_outside_training_by_release_date",
        "openfold3_preview_cutoff_confidence",
    ]


def author_ligand_location(label: Any, ccd: str | None) -> tuple[str | None, int | None]:
    """Parse canonical crystal-manifest labels such as ``A:BHF800``."""
    if ccd is None or label is None or pd.isna(label):
        return None, None
    chain, separator, residue = str(label).partition(":")
    if not separator or not residue.startswith(ccd):
        raise ValueError(f"invalid_ligand_label:{label}:ccd={ccd}")
    number = residue[len(ccd) :]
    if not number.lstrip("-").isdigit():
        raise ValueError(f"invalid_ligand_residue_number:{label}:ccd={ccd}")
    return chain, int(number)


def ligand_metadata(system: pd.DataFrame, inventory: pd.DataFrame) -> pd.DataFrame:
    extra = inventory.rename(
        columns={
            "ccd": "selected_ccd",
            "name": "inventory_ligand_name",
            "smiles": "inventory_smiles",
            "formal_charge": "inventory_formal_charge",
            "heavy_atoms": "inventory_heavy_atoms",
            "closest_fe_A": "inventory_closest_fe_A",
        }
    )[
        [
            "pdb_id",
            "selected_ccd",
            "inventory_ligand_name",
            "inventory_smiles",
            "inventory_formal_charge",
            "inventory_heavy_atoms",
            "formula",
            "mw",
            "chemistry_eligible",
            "chemistry_reason",
            "pdb_heavy_atoms",
            "inventory_closest_fe_A",
        ]
    ]
    if extra.duplicated(["pdb_id", "selected_ccd"]).any():
        raise ValueError("ligand_inventory_has_duplicate_pdb_ccd_keys")
    merged = system.merge(extra, on=["pdb_id", "selected_ccd"], how="left", validate="one_to_one")
    for primary, fallback in [
        ("selected_name", "inventory_ligand_name"),
        ("smiles", "inventory_smiles"),
        ("formal_charge", "inventory_formal_charge"),
        ("heavy_atoms", "inventory_heavy_atoms"),
        ("closest_fe_A", "inventory_closest_fe_A"),
    ]:
        merged[primary] = merged[primary].combine_first(merged[fallback])
    return merged.drop(
        columns=["inventory_ligand_name", "inventory_smiles", "inventory_formal_charge", "inventory_heavy_atoms", "inventory_closest_fe_A"]
    )


def hydrate_completed_pilot_systems(system: pd.DataFrame, truth_scores: pd.DataFrame) -> pd.DataFrame:
    """Restore ligand identities deliberately omitted from the remaining-139 campaign manifest."""
    output = system.copy()
    selected = truth_scores.set_index("pdb_id")["selected_ccd"].to_dict()
    mask = output.status.eq("already_complete_panel") & output.selected_ccd.isna()
    for index in output.index[mask]:
        pdb_id = str(output.at[index, "pdb_id"])
        ccd = selected.get(pdb_id)
        if not ccd or pd.isna(ccd):
            raise ValueError(f"missing_completed_pilot_ligand:{pdb_id}")
        gene = str(output.at[index, "gene"])
        output.at[index, "selected_ccd"] = str(ccd)
        output.at[index, "system_id"] = f"{gene.lower()}_{pdb_id.lower()}_{str(ccd).lower()}"
        output.at[index, "status"] = "selected_nearest_heme_completed_pilot"
        output.at[index, "candidate_ccds"] = str(ccd)
        output.at[index, "eligible_candidate_ccds"] = str(ccd)
        output.at[index, "selection_rule"] = "completed 10-PDB pilot ligand; deposited selected CCD restored from truth-validation ledger"
    if int(mask.sum()) != 10:
        raise ValueError(f"unexpected_completed_pilot_count:{int(mask.sum())}")
    return output


def curate_experimental_target_roles(systems: pd.DataFrame, inventory: pd.DataFrame) -> pd.DataFrame:
    """Apply the single obvious target-vs-additive correction in the CYP universe."""
    output = systems.copy()
    target = inventory[(inventory.pdb_id == "6CSB") & (inventory.ccd == "RTZ")]
    if len(target) != 1:
        raise ValueError("6CSB_RTZ_inventory_gate_failed")
    ligand = target.iloc[0]
    index = output.index[output.pdb_id.eq("6CSB")]
    if len(index) != 1 or output.loc[index[0], "selected_ccd"] != "2CV":
        raise ValueError("6CSB_automated_2CV_selection_gate_failed")
    i = index[0]
    output.at[i, "system_id"] = "cyp2d6_6csb_rtz"
    output.at[i, "selected_ccd"] = "RTZ"
    output.at[i, "selected_name"] = ligand["name"]
    output.at[i, "smiles"] = ligand["smiles"]
    output.at[i, "heavy_atoms"] = ligand["heavy_atoms"]
    output.at[i, "formal_charge"] = ligand["formal_charge"]
    output.at[i, "closest_fe_A"] = 8.50342284024498
    output.at[i, "formula"] = ligand["formula"]
    output.at[i, "mw"] = ligand["mw"]
    output.at[i, "chemistry_eligible"] = ligand["chemistry_eligible"]
    output.at[i, "chemistry_reason"] = ligand["chemistry_reason"]
    output.at[i, "pdb_heavy_atoms"] = ligand["pdb_heavy_atoms"]
    output.at[i, "status"] = "role_curated_primary_ligand"
    output.at[i, "selection_rule"] = "manual target-role curation: complete RTZ chain A residue 602; automated 2CV choice is a partial HEGA-10 detergent fragment"
    return output


def copy_verified(source: Path, destination: Path, expected_sha256: str) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and sha256(destination) == expected_sha256:
        return
    temporary = destination.with_suffix(destination.suffix + ".building")
    shutil.copyfile(source, temporary)
    observed = sha256(temporary)
    if observed != expected_sha256:
        temporary.unlink(missing_ok=True)
        raise ValueError(f"copy_hash_mismatch:{source}:{observed}!={expected_sha256}")
    os.replace(temporary, destination)


def build_experimental(
    release: Path,
    systems: pd.DataFrame,
    release_dates: pd.DataFrame,
    training: pd.DataFrame,
    crystal_metrics: pd.DataFrame,
    crystal_manifest: pd.DataFrame,
    truth_scores: pd.DataFrame,
    aligner: CanonicalAligner,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    metadata = systems.merge(
        release_dates[
            ["pdb_id", "method", "resolution", "coordinate_format", "source_url", "bytes", "sha256", "note"]
        ].rename(columns={"method": "experimental_method", "source_url": "rcsb_source_url", "sha256": "deposited_sha256", "bytes": "deposited_bytes"}),
        on="pdb_id",
        validate="one_to_one",
    )
    metadata = metadata.merge(
        training[["pdb_id", *release_training_columns()]], on="pdb_id", validate="one_to_one"
    )
    metadata = metadata.merge(crystal_metrics, on=["pdb_id", "gene"], validate="one_to_one")
    cyp_manifest = crystal_manifest[crystal_manifest.dataset.str.startswith("CYP")][
        ["structure_id", "selected_chain", "spacegroup", "unit_cell", "ligand_label", "ligand_atom_count"]
    ].rename(columns={"structure_id": "pdb_id"})
    metadata = metadata.merge(cyp_manifest, on="pdb_id", validate="one_to_one")
    metadata = metadata.merge(truth_scores.add_prefix("truth_score_"), left_on="pdb_id", right_on="truth_score_pdb_id", how="left", validate="one_to_one")

    rows: list[dict[str, Any]] = []
    transforms: list[dict[str, Any]] = []
    for index, row in enumerate(metadata.to_dict("records"), start=1):
        pdb_id = str(row["pdb_id"]).upper()
        gene = str(row["gene"])
        source = Path(str(row["truth_path"]))
        relative = Path("structures") / "experimental" / gene.lower() / f"{pdb_id.lower()}.cif.gz"
        transform = aligner.align(source, gene, release / relative)
        transform["source_path"] = source.relative_to(REPO).as_posix()
        transform["reference_path"] = aligner.reference_path.relative_to(REPO).as_posix()
        transform["anchor_map_path"] = aligner.anchor_map_path.relative_to(REPO).as_posix()
        transform["pdb_id"] = pdb_id
        transform["aligned_relative_path"] = relative.as_posix()
        transforms.append(transform)
        selected_ccd = None if pd.isna(row.get("selected_ccd")) else str(row["selected_ccd"])
        if pdb_id == "6CSB" and selected_ccd == "RTZ":
            ligand_chain, ligand_residue_number = "A", 602
        else:
            ligand_chain, ligand_residue_number = author_ligand_location(row.get("ligand_label"), selected_ccd)
        truth_status = row.get("truth_score_status")
        truth_score_matches_selection = row.get("truth_score_selected_ccd") == selected_ccd
        record: dict[str, Any] = {
            "structure_id": f"exp:{pdb_id.lower()}",
            "source_kind": "experimental",
            "method": "experimental",
            "pdb_id": pdb_id,
            "gene": gene,
            "accession": row["accession"],
            "source_system_id": row["system_id"],
            "seed": pd.NA,
            "sample_index": pd.NA,
            "coordinate_path": relative.as_posix(),
            "coordinate_sha256": transform["aligned_sha256"],
            "coordinate_source_sha256": row["truth_sha256"],
            "coordinate_source_url": row["rcsb_source_url"],
            "coordinate_format": "mmCIF.gz",
            "experimental_source_format": row["coordinate_format"],
            "experimental_method": row["experimental_method"],
            "resolution": row["resolution"],
            "deposited_bytes": row["deposited_bytes"],
            "ligand_selected_ccd": selected_ccd,
            "ligand_role": "primary_target_role_curated" if pdb_id == "6CSB" else ("none" if selected_ccd is None else "primary_selected_ligand"),
            "ground_truth_target_ligand_ccd": selected_ccd,
            "ost_target_role_valid": True,
            "ligand_name": row.get("selected_name"),
            "ligand_smiles": row.get("smiles"),
            "ligand_formula": row.get("formula"),
            "ligand_molecular_weight": row.get("mw"),
            "ligand_formal_charge": row.get("formal_charge"),
            "ligand_expected_heavy_atoms": row.get("heavy_atoms"),
            "ligand_deposited_heavy_atoms": row.get("ligand_atom_count"),
            "ligand_chemistry_eligible": row.get("chemistry_eligible"),
            "ligand_chemistry_reason": row.get("chemistry_reason"),
            "ligand_selection_status": row.get("status"),
            "ligand_selection_rule": row.get("selection_rule"),
            "ligand_candidate_ccds": row.get("candidate_ccds"),
            "ligand_eligible_candidate_ccds": row.get("eligible_candidate_ccds"),
            "ligand_closest_heme_distance_A": row.get("closest_fe_A"),
            "ligand_chain": ligand_chain,
            "ligand_residue_number": ligand_residue_number,
            "protein_chain": row.get("selected_chain"),
            "alignment_invariant": ALIGNMENT_INVARIANT,
            "alignment_status": "PASS",
            "alignment_numbering_mode": transform["numbering_mode"],
            "alignment_anchor_pairs": transform["anchor_pairs"],
            "alignment_anchor_rmsd_A": transform["anchor_rmsd_A"],
            "ost_status": (
                "NOT_APPLICABLE_NO_SELECTED_LIGAND" if selected_ccd is None
                else truth_status if truth_score_matches_selection
                else "NOT_CALCULATED_ROLE_CURATED_TARGET_CHANGED"
            ),
            "ost_lddt_pli": row.get("truth_score_ost_lddt_pli") if truth_score_matches_selection else math.nan,
            "ost_lddt_pli_n_contacts": row.get("truth_score_ost_lddt_pli_n_contacts") if truth_score_matches_selection else math.nan,
            "ost_error": None if selected_ccd is None else (row.get("truth_score_error") if truth_score_matches_selection else "prior truth self-score used automated 2CV additive selection"),
            "pl_iptm": math.nan,
            "native_ligand_iptm": math.nan,
            "native_iptm": math.nan,
            "spacegroup": row.get("spacegroup"),
            "unit_cell": row.get("unit_cell"),
            "heme_fe_to_1tqn_A": row.get("fe_to_1tqn_A"),
            "heme_normal_to_1tqn_deg": row.get("heme_normal_to_1tqn_deg"),
            "axial_cys_sg_to_1tqn_A": row.get("axial_cys_sg_to_1tqn_A"),
            "axial_cys_resnum": row.get("axial_cys_resnum"),
            "fe_cys_A": row.get("fe_cys_A"),
            "coverage_status": "ground_truth_all_current_uniprot_pdb_crossrefs",
        }
        for column in release_training_columns():
            record[column] = row.get(column)
        rows.append(record)
        if index % 25 == 0 or index == len(metadata):
            print(json.dumps({"event": "experimental_alignment_progress", "done": index, "total": len(metadata)}), flush=True)
    return rows, transforms


def build_predictions(
    release: Path,
    systems: pd.DataFrame,
    release_dates: pd.DataFrame,
    training: pd.DataFrame,
) -> list[dict[str, Any]]:
    model_path = REPO / "analysis" / "cyp_cofold_postprocess_20260804" / "model_manifest.csv"
    alignment_path = REPO / "analysis" / "cyp_cofold_postprocess_20260804" / "alignment_results.csv"
    scores_path = REPO / "analysis" / "cyp_cofold_postprocess_20260804" / "ost_lddt_pli_all_models.csv"
    models = pd.read_csv(model_path, low_memory=False)
    eligible = models[models.score_eligible.fillna(False)].copy()
    if len(eligible) != 16630:
        raise ValueError(f"unexpected_score_eligible_count:{len(eligible)}")
    alignment = pd.read_csv(alignment_path, low_memory=False).rename(
        columns={
            "status": "alignment_status",
            "error": "alignment_error",
            "anchor_pairs": "alignment_anchor_pairs",
            "anchor_rmsd_A": "alignment_anchor_rmsd_A",
            "model_ligand_heavy_atoms": "aligned_model_ligand_heavy_atoms",
            "source_sha256": "alignment_source_sha256",
            "aligned_sha256": "coordinate_sha256",
            "aligned_path": "verified_aligned_path",
        }
    )
    valid_alignment = alignment[alignment.alignment_status.eq("PASS")].copy()
    combined = eligible.merge(valid_alignment, on="model_path", validate="one_to_one", suffixes=("", "_alignment"))
    scores = pd.read_csv(scores_path, low_memory=False)[
        [
            "model_path",
            "status",
            "ost_state",
            "ost_lddt_pli",
            "ost_lddt_pli_n_contacts",
            "chain_mapping",
            "truth_ligand_chain_selected",
            "truth_ligand_resnum_selected",
            "model_ligand_atoms",
            "error",
        ]
    ].rename(columns={"status": "ost_status", "error": "ost_error"})
    combined = combined.merge(scores, on="model_path", validate="one_to_one")
    if len(combined) != 16590:
        raise ValueError(f"unexpected_valid_prediction_count:{len(combined)}")

    ligand_columns = [
        "pdb_id",
        "accession",
        "selected_ccd",
        "selected_name",
        "smiles",
        "formula",
        "mw",
        "formal_charge",
        "heavy_atoms",
        "chemistry_eligible",
        "chemistry_reason",
        "status",
        "selection_rule",
        "candidate_ccds",
        "eligible_candidate_ccds",
        "closest_fe_A",
    ]
    combined = combined.merge(systems[ligand_columns], on=["pdb_id", "selected_ccd"], validate="many_to_one", suffixes=("", "_system"))
    combined = combined.merge(
        release_dates[["pdb_id", "source_url", "resolution", "method"]].rename(
            columns={"method": "experimental_method"}
        ),
        on="pdb_id",
        validate="many_to_one",
    )
    combined = combined.merge(
        training[["pdb_id", *release_training_columns()]], on="pdb_id", validate="many_to_one"
    )
    if combined.duplicated(["method", "pdb_id", "selected_ccd", "seed", "sample"]).any():
        duplicate = combined.loc[
            combined.duplicated(["method", "pdb_id", "selected_ccd", "seed", "sample"], keep=False),
            ["method", "pdb_id", "selected_ccd", "seed", "sample", "system_id"],
        ]
        raise ValueError(f"non_unique_prediction_identity:{duplicate.head().to_dict('records')}")

    rows: list[dict[str, Any]] = []
    for index, row in enumerate(combined.to_dict("records"), start=1):
        method = str(row["method"])
        gene = str(row["gene"])
        pdb_id = str(row["pdb_id"]).upper()
        ccd = str(row["selected_ccd"])
        seed = int(row["seed"])
        sample = int(row["sample"])
        system_slug = f"{gene.lower()}_{pdb_id.lower()}_{ccd.lower()}"
        relative = (
            Path("structures")
            / "predicted"
            / method
            / system_slug
            / f"seed_{seed}_sample_{sample}.cif.gz"
        )
        source = Path(str(row["verified_aligned_path"]))
        expected = str(row["coordinate_sha256"])
        copy_verified(source, release / relative, expected)
        record: dict[str, Any] = {
            "structure_id": f"pred:{method}:{pdb_id.lower()}:{ccd.lower()}:seed{seed}:sample{sample}",
            "source_kind": "predicted",
            "method": method,
            "pdb_id": pdb_id,
            "gene": gene,
            "accession": row["accession"],
            "source_system_id": row["system_id"],
            "seed": seed,
            "sample_index": sample,
            "coordinate_path": relative.as_posix(),
            "coordinate_sha256": expected,
            "coordinate_source_sha256": row["model_sha256"],
            "coordinate_source_url": None,
            "coordinate_format": "mmCIF.gz",
            "experimental_source_format": None,
            "experimental_method": row["experimental_method"],
            "resolution": row["resolution"],
            "deposited_bytes": pd.NA,
            "ligand_selected_ccd": ccd,
            "ligand_role": "crystallization_additive_automated_selection" if pdb_id == "6CSB" and ccd == "2CV" else "primary_selected_ligand",
            "ground_truth_target_ligand_ccd": "RTZ" if pdb_id == "6CSB" and ccd == "2CV" else ccd,
            "ost_target_role_valid": not (pdb_id == "6CSB" and ccd == "2CV"),
            "ligand_name": row.get("selected_name"),
            "ligand_smiles": row.get("smiles"),
            "ligand_formula": row.get("formula"),
            "ligand_molecular_weight": row.get("mw"),
            "ligand_formal_charge": row.get("formal_charge"),
            "ligand_expected_heavy_atoms": row.get("heavy_atoms"),
            "ligand_deposited_heavy_atoms": row.get("aligned_model_ligand_heavy_atoms"),
            "ligand_chemistry_eligible": row.get("chemistry_eligible"),
            "ligand_chemistry_reason": row.get("chemistry_reason"),
            "ligand_selection_status": row.get("status_system"),
            "ligand_selection_rule": row.get("selection_rule"),
            "ligand_candidate_ccds": row.get("candidate_ccds"),
            "ligand_eligible_candidate_ccds": row.get("eligible_candidate_ccds"),
            "ligand_closest_heme_distance_A": row.get("closest_fe_A"),
            "ligand_chain": row.get("model_ligand_chain"),
            "ligand_residue_number": 1,
            "protein_chain": "A",
            "alignment_invariant": ALIGNMENT_INVARIANT,
            "alignment_status": row["alignment_status"],
            "alignment_numbering_mode": "core-index",
            "alignment_anchor_pairs": row["alignment_anchor_pairs"],
            "alignment_anchor_rmsd_A": row["alignment_anchor_rmsd_A"],
            "ost_status": row["ost_status"],
            "ost_lddt_pli": row["ost_lddt_pli"],
            "ost_lddt_pli_n_contacts": row["ost_lddt_pli_n_contacts"],
            "ost_error": row["ost_error"],
            "pl_iptm": row.get("pl_iptm"),
            "native_ligand_iptm": row.get("native_ligand_iptm"),
            "native_iptm": row.get("native_iptm"),
            "spacegroup": None,
            "unit_cell": None,
            "heme_fe_to_1tqn_A": math.nan,
            "heme_normal_to_1tqn_deg": math.nan,
            "axial_cys_sg_to_1tqn_A": math.nan,
            "axial_cys_resnum": math.nan,
            "fe_cys_A": math.nan,
            "coverage_status": "valid_aligned_ligand_prediction",
        }
        for column in release_training_columns():
            record[column] = row.get(column)
        rows.append(record)
        if index % 1000 == 0 or index == len(combined):
            print(json.dumps({"event": "prediction_copy_progress", "done": index, "total": len(combined)}), flush=True)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release", type=Path, default=DEFAULT_RELEASE)
    parser.add_argument("--work", type=Path, default=DEFAULT_WORK)
    args = parser.parse_args()
    release = args.release.resolve()
    work = args.work.resolve()
    release.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)

    system_path = REPO / "runs" / "cyp_all_remaining139_four_engine_seed1_5sample_20260730" / "system_manifest.csv"
    inventory_path = REPO / "runs" / "cyp_all_remaining139_four_engine_seed1_5sample_20260730" / "ligand_inventory.csv"
    release_dates_path = REPO / "intro" / "pdbs" / "manifest_with_release_dates.csv"
    training_path = REPO / "analysis" / "cyp_crystal_pocket_variability_audit_20260802" / "cyp_crystal_deposition_training_status.csv"
    crystal_metrics_path = REPO / "analysis" / "canonical_alignment_20260731" / "crystal_alignment_metrics.csv"
    crystal_manifest_path = REPO / "analysis" / "cyp_crystal_pocket_variability_audit_20260802" / "structure_manifest.csv"
    truth_scores_path = REPO / "analysis" / "cyp_cofold_postprocess_20260804" / "truth_self_validation.csv"
    anchor_map_path = REPO / "analysis" / "canonical_alignment_20260731" / "anchor_map.csv"
    reference_path = REPO / "intro" / "pdbs" / "CYP3A4" / "1TQN.pdb"

    truth_scores = consistent_truth_selection(truth_scores_path)
    raw_systems = hydrate_completed_pilot_systems(
        pd.read_csv(system_path, low_memory=False), truth_scores
    )
    inventory = pd.read_csv(inventory_path, low_memory=False)
    systems = ligand_metadata(raw_systems, inventory)
    experimental_systems = curate_experimental_target_roles(systems, inventory)
    release_dates = pd.read_csv(release_dates_path, low_memory=False)
    training = pd.read_csv(training_path, low_memory=False)
    crystal_metrics = pd.read_csv(crystal_metrics_path, low_memory=False)
    crystal_manifest = pd.read_csv(crystal_manifest_path, low_memory=False)
    if len(systems) != 149 or systems.pdb_id.nunique() != 149:
        raise ValueError("system_manifest_universe_gate_failed")

    aligner = CanonicalAligner(reference_path, anchor_map_path)
    experimental, transforms = build_experimental(
        release,
        experimental_systems,
        release_dates,
        training,
        crystal_metrics,
        crystal_manifest,
        truth_scores,
        aligner,
    )
    predictions = build_predictions(release, systems, release_dates, training)
    base = pd.DataFrame(experimental + predictions)
    base["seed"] = base["seed"].astype("Int64")
    base["sample_index"] = base["sample_index"].astype("Int64")
    base = base.sort_values(["source_kind", "method", "gene", "pdb_id", "seed", "sample_index"], na_position="first").reset_index(drop=True)
    if len(base) != 16739 or base.structure_id.nunique() != 16739 or base.coordinate_path.nunique() != 16739:
        raise ValueError(
            f"base_universe_gate_failed:rows={len(base)}:ids={base.structure_id.nunique()}:paths={base.coordinate_path.nunique()}"
        )
    missing = [path for path in base.coordinate_path if not (release / path).is_file()]
    if missing:
        raise ValueError(f"missing_release_coordinates:{missing[:10]}")

    atomic_parquet(base, work / "base_structures.parquet")
    transforms_path = release / "provenance" / "experimental_alignment_transforms.jsonl"
    transforms_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = transforms_path.with_suffix(".jsonl.building")
    with temporary.open("w") as handle:
        for transform in transforms:
            handle.write(json.dumps(transform, sort_keys=True) + "\n")
    os.replace(temporary, transforms_path)
    shutil.copy2(anchor_map_path, release / "provenance" / "anchor_map.csv")
    summary = {
        "created_at": now(),
        "rows": len(base),
        "experimental_rows": int(base.source_kind.eq("experimental").sum()),
        "predicted_rows": int(base.source_kind.eq("predicted").sum()),
        "unique_pdb_ids": int(base.pdb_id.nunique()),
        "method_counts": base.method.value_counts().sort_index().to_dict(),
        "selected_ligand_pdb_ids": int(systems.selected_ccd.notna().sum()),
        "apo_or_heme_only_pdb_ids": int(systems.selected_ccd.isna().sum()),
        "alignment_invariant": ALIGNMENT_INVARIANT,
        "base_parquet": str(work / "base_structures.parquet"),
    }
    atomic_json(summary, work / "build_summary.json")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
