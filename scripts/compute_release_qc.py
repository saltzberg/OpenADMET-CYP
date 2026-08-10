#!/usr/bin/env python3
"""Compute authoritative-topology ProLIF and PoseBusters QC for every release row.

The calculation is partitioned into deterministic Parquet shards and can be
resumed safely.  Coordinates are never optimized or otherwise moved.
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import math
import multiprocessing as mp
import os
import tempfile
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import gemmi
import networkx as nx
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger


REPO = Path(__file__).resolve().parents[1]
DEFAULT_RELEASE = REPO / "release" / "huggingface_cyp_cofold_v1"
DEFAULT_WORK = REPO / "work" / "cyp_hf_release_v1"
QC_PROTOCOL = "cyp_release_prolif_posebusters_authoritative_ccd_topology_v2"
PROLIF_PROTOCOL = "ProLIF_2.1.0_default_fingerprint_count_true_authoritative_CCD_SMILES_RDKit_sequence_template_CYP_protein_only_explicit_Ir2.00_Ru2.05_vdw_radii"
POSEBUSTERS_PROTOCOL = "PoseBusters_0.6.5_dock_full_report_authoritative_CCD_SMILES_pose_CYP_protein_only"
AA3 = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL MSE SEC PYL ASX GLX UNK".split())
AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "MSE": "M",
    "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
INTERACTIONS = [
    "Hydrophobic", "HBDonor", "HBAcceptor", "PiStacking", "Anionic", "Cationic", "CationPi", "PiCation", "VdWContact"
]
POSEBUSTERS_CORE_CHECKS = [
    "mol_pred_loaded",
    "mol_cond_loaded",
    "sanitization",
    "inchi_convertible",
    "all_atoms_connected",
    "no_radicals",
    "bond_lengths",
    "bond_angles",
    "internal_steric_clash",
    "aromatic_ring_flatness",
    "non-aromatic_ring_non-flatness",
    "double_bond_flatness",
    "internal_energy",
    "protein-ligand_maximum_distance",
    "minimum_distance_to_protein",
    "volume_overlap_with_protein",
]
RADII = {
    "H": 0.31, "B": 0.85, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07,
    "S": 1.05, "CL": 1.02, "BR": 1.20, "I": 1.39, "SE": 1.20, "IR": 1.41, "RU": 1.46,
}
UNIPROT_START = {"CYP1A2": 27, "CYP2C9": 30, "CYP2D6": 34, "CYP3A4": 30}
_RELEASE: Path | None = None
_BUSTER: Any = None


@dataclass(frozen=True)
class Atom:
    record: str
    name: str
    residue: str
    chain: str
    residue_number: int
    insertion_code: str
    element: str
    coord: np.ndarray


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


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


def normalize(value: Any) -> Any:
    if value is pd.NA:
        return None
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if math.isnan(float(value)) else float(value)
    if isinstance(value, float) and math.isnan(value):
        return None
    if isinstance(value, Path):
        return str(value)
    return value


def initialize_worker(release: str) -> None:
    global _RELEASE, _BUSTER
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    RDLogger.DisableLog("rdApp.*")
    warnings.filterwarnings("ignore")
    _RELEASE = Path(release)
    from posebusters import PoseBusters

    _BUSTER = PoseBusters(config="dock", max_workers=0, chunk_size=10)


def altloc_priority(altloc: str) -> int | None:
    if altloc in {"\x00", " ", "", ".", "?"}:
        return 0
    if altloc == "A":
        return 1
    if altloc == "1":
        return 2
    return None


def parse_atoms(path: Path) -> list[Atom]:
    structure = gemmi.read_structure(str(path))
    if not len(structure):
        raise ValueError("coordinate_file_has_no_models")
    chosen: dict[tuple[str, int, str, str, str], tuple[int, Atom]] = {}
    for chain in structure[0]:
        for residue in chain:
            is_polymer = residue.entity_type == gemmi.EntityType.Polymer or str(residue.het_flag) == "A"
            record = "ATOM" if is_polymer and residue.name.upper() in AA3 else "HETATM"
            insertion = str(residue.seqid.icode).strip()
            for atom in residue:
                priority = altloc_priority(str(atom.altloc))
                if priority is None or atom.element.atomic_number <= 1:
                    continue
                item = Atom(
                    record=record,
                    name=atom.name.strip(),
                    residue=residue.name.upper(),
                    chain=chain.name,
                    residue_number=int(residue.seqid.num) if residue.seqid.num is not None else 1,
                    insertion_code=insertion,
                    element=atom.element.name.upper(),
                    coord=np.array([atom.pos.x, atom.pos.y, atom.pos.z], dtype=float),
                )
                key = (item.chain, item.residue_number, item.insertion_code, item.residue, item.name)
                previous = chosen.get(key)
                if previous is None or priority < previous[0]:
                    chosen[key] = (priority, item)
    return [entry[1] for entry in chosen.values()]


def select_components(row: dict[str, Any], atoms: list[Atom]) -> tuple[list[Atom], list[Atom]]:
    protein_chain = str(row["protein_chain"])
    protein = [atom for atom in atoms if atom.record == "ATOM" and atom.chain == protein_chain and atom.residue in AA1]
    if not protein:
        raise ValueError(f"no_protein_heavy_atoms:chain={protein_chain}")

    ligand_chain = str(row["ligand_chain"])
    if str(row["source_kind"]) == "experimental":
        ligand_number = int(row["ligand_residue_number"])
        ligand_ccd = str(row["ligand_selected_ccd"])
        ligand = [
            atom for atom in atoms
            if atom.record == "HETATM" and atom.chain == ligand_chain
            and atom.residue_number == ligand_number and atom.residue == ligand_ccd
        ]
    else:
        ligand = [atom for atom in atoms if atom.record == "HETATM" and atom.chain == ligand_chain]
    if not ligand:
        raise ValueError(
            f"no_selected_ligand_heavy_atoms:chain={ligand_chain}:resnum={row.get('ligand_residue_number')}:ccd={row.get('ligand_selected_ccd')}"
        )
    return protein, ligand


def rdkit_elements_adjacency(smiles: str) -> tuple[Chem.Mol, list[str], dict[int, set[int]]]:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("authoritative_smiles_parse_failed")
    molecule = Chem.RemoveAllHs(molecule)
    elements = [atom.GetSymbol().upper() for atom in molecule.GetAtoms()]
    adjacency = {index: set() for index in range(molecule.GetNumAtoms())}
    for bond in molecule.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        adjacency[begin].add(end)
        adjacency[end].add(begin)
    return molecule, elements, adjacency


def infer_adjacency(elements: list[str], coordinates: np.ndarray) -> dict[int, set[int]]:
    adjacency = {index: set() for index in range(len(elements))}
    for i in range(len(elements)):
        for j in range(i + 1, len(elements)):
            distance = float(np.linalg.norm(coordinates[i] - coordinates[j]))
            upper = RADII.get(elements[i], 0.77) + RADII.get(elements[j], 0.77) + 0.45
            if 0.45 <= distance <= upper:
                adjacency[i].add(j)
                adjacency[j].add(i)
    return adjacency


def bond_upper(element_a: str, element_b: str) -> float:
    pair = frozenset((element_a.upper(), element_b.upper()))
    if pair & {"IR", "RU"}:
        return 2.70
    if pair == frozenset(("C", "I")):
        return 2.25
    if pair == frozenset(("P", "S")):
        return 2.15
    if "I" in pair or "BR" in pair:
        return 2.20
    return 2.05


def bonded_distance_audit(
    coordinates: np.ndarray, elements: list[str], adjacency: dict[int, set[int]]
) -> tuple[int, float, float]:
    distances: list[float] = []
    failures = 0
    for i, neighbors in adjacency.items():
        for j in neighbors:
            if i >= j:
                continue
            distance = float(np.linalg.norm(coordinates[i] - coordinates[j]))
            distances.append(distance)
            if distance < 0.75 or distance > bond_upper(elements[i], elements[j]):
                failures += 1
    return failures, min(distances, default=math.nan), max(distances, default=math.nan)


def authoritative_mapping(
    elements: list[str], coordinates: np.ndarray, smiles: str, coordinate_adjacency: dict[int, set[int]] | None = None
) -> tuple[Chem.Mol, list[int], str, dict[int, set[int]], float, float]:
    molecule, reference_elements, reference_adjacency = rdkit_elements_adjacency(smiles)
    if len(reference_elements) != len(elements):
        raise ValueError(
            f"heavy_atom_count_mismatch:authoritative={len(reference_elements)}:coordinate={len(elements)}"
        )
    if sorted(reference_elements) != sorted(elements):
        raise ValueError("heavy_atom_element_multiset_mismatch")
    inferred = coordinate_adjacency if coordinate_adjacency is not None else infer_adjacency(elements, coordinates)
    reference_graph = nx.Graph()
    coordinate_graph = nx.Graph()
    for index, element in enumerate(reference_elements):
        reference_graph.add_node(index, element=element)
    for index, element in enumerate(elements):
        coordinate_graph.add_node(index, element=element)
    for i, neighbors in reference_adjacency.items():
        for j in neighbors:
            if i < j:
                reference_graph.add_edge(i, j)
    for i, neighbors in inferred.items():
        for j in neighbors:
            if i < j:
                coordinate_graph.add_edge(i, j)
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        reference_graph,
        coordinate_graph,
        node_match=lambda left, right: left["element"] == right["element"],
    )
    candidates: list[tuple[int, float, float, list[int]]] = []
    for mapping in matcher.isomorphisms_iter():
        order = [mapping[index] for index in range(len(reference_elements))]
        bad, minimum, maximum = bonded_distance_audit(
            coordinates[order], reference_elements, reference_adjacency
        )
        candidates.append((bad, maximum, minimum, order))
        if len(candidates) >= 256:
            break
    if candidates:
        candidates.sort(key=lambda item: (item[0], item[1], -item[2], item[3]))
        bad, maximum, minimum, order = candidates[0]
        if bad == 0:
            return molecule, order, "inferred_graph_isomorphism", reference_adjacency, minimum, maximum
    if elements == reference_elements:
        bad, minimum, maximum = bonded_distance_audit(coordinates, reference_elements, reference_adjacency)
        if bad == 0:
            return (
                molecule,
                list(range(len(elements))),
                "authoritative_element_order_identity",
                reference_adjacency,
                minimum,
                maximum,
            )
        return (
            molecule,
            list(range(len(elements))),
            f"authoritative_element_order_identity_with_{bad}_bond_geometry_outliers",
            reference_adjacency,
            minimum,
            maximum,
        )
    if candidates:
        return (
            molecule,
            order,
            f"inferred_graph_isomorphism_with_{bad}_bond_geometry_outliers",
            reference_adjacency,
            minimum,
            maximum,
        )
    raise ValueError("no_geometry_valid_authoritative_atom_mapping")


def deposited_pdb_conect_adjacency(row: dict[str, Any], ligand: list[Atom]) -> dict[int, set[int]]:
    """Recover exact selected-ligand connectivity from legacy PDB CONECT records."""
    source = REPO / "intro" / "pdbs" / str(row["gene"]) / f"{str(row['pdb_id']).upper()}.pdb"
    if not source.is_file():
        raise ValueError(f"deposited_pdb_not_available:{source}")
    chain = str(row["ligand_chain"])
    ccd = str(row["ligand_selected_ccd"])
    residue_number = int(row["ligand_residue_number"])
    serial_by_name: dict[str, int] = {}
    conect_lines: list[str] = []
    with source.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("CONECT"):
                conect_lines.append(line)
                continue
            if not line.startswith("HETATM") or line[16].strip() not in {"", "A", "1"}:
                continue
            try:
                number = int(line[22:26])
                serial = int(line[6:11])
            except ValueError:
                continue
            if line[21].strip() == chain and number == residue_number and line[17:20].strip().upper() == ccd:
                element = line[76:78].strip().upper()
                if element not in {"H", "D"}:
                    serial_by_name[line[12:16].strip()] = serial
    if set(serial_by_name) != {atom.name for atom in ligand}:
        raise ValueError("deposited_PDB_CONECT_atom_name_mismatch")
    index_by_serial = {serial_by_name[atom.name]: index for index, atom in enumerate(ligand)}
    adjacency = {index: set() for index in range(len(ligand))}
    for line in conect_lines:
        try:
            serials = [int(line[index : index + 5]) for index in range(6, len(line), 5) if line[index : index + 5].strip()]
        except ValueError:
            continue
        if not serials or serials[0] not in index_by_serial:
            continue
        begin = index_by_serial[serials[0]]
        for serial in serials[1:]:
            if serial in index_by_serial:
                end = index_by_serial[serial]
                adjacency[begin].add(end)
                adjacency[end].add(begin)
    if not any(adjacency.values()):
        raise ValueError("selected_ligand_has_no_deposited_PDB_CONECT_bonds")
    return adjacency


def ligand_molecule(smiles: str, ligand: list[Atom], coordinate_adjacency: dict[int, set[int]] | None = None) -> tuple[Chem.Mol, dict[str, Any]]:
    coordinates = np.vstack([atom.coord for atom in ligand])
    elements = [atom.element for atom in ligand]
    molecule, order, mapping_mode, adjacency, minimum, maximum = authoritative_mapping(
        elements, coordinates, smiles, coordinate_adjacency
    )
    ordered = coordinates[order]
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for index, point in enumerate(ordered):
        conformer.SetAtomPosition(index, tuple(float(value) for value in point))
    molecule.RemoveAllConformers()
    molecule.AddConformer(conformer)
    metadata = {
        "ligand_topology_status": "ok",
        "ligand_atom_mapping_mode": mapping_mode,
        "ligand_coordinate_heavy_atoms": len(ligand),
        "ligand_authoritative_heavy_atoms": molecule.GetNumAtoms(),
        "ligand_authoritative_bonds": sum(len(neighbors) for neighbors in adjacency.values()) // 2,
        "ligand_authoritative_min_bond_distance_A": minimum,
        "ligand_authoritative_max_bond_distance_A": maximum,
        "ligand_coordinate_max_abs_delta_A": 0.0,
    }
    return molecule, metadata


def protein_rdkit_molecule(
    protein: list[Atom], source_kind: str, gene: str
) -> tuple[Chem.Mol, str, int]:
    residue_order: list[tuple[int, str, str]] = []
    by_residue: dict[tuple[int, str, str], dict[str, Atom]] = {}
    for atom in protein:
        key = (atom.residue_number, atom.insertion_code, atom.residue)
        if key not in by_residue:
            residue_order.append(key)
            by_residue[key] = {}
        by_residue[key][atom.name] = atom
    sequence = "".join(AA1.get(key[2], "") for key in residue_order)
    if len(sequence) != len(residue_order):
        unsupported = sorted({key[2] for key in residue_order if key[2] not in AA1})
        raise ValueError(f"unsupported_protein_residue_for_sequence_template:{unsupported}")
    template = Chem.MolFromFASTA(sequence)
    if template is None:
        raise ValueError("protein_sequence_template_failed")
    observed_by_position = {
        position: by_residue[key] for position, key in enumerate(residue_order, start=1)
    }
    if source_kind == "predicted":
        canonical_numbers = {
            position: UNIPROT_START[gene] + position - 1 for position in observed_by_position
        }
        numbering_mode = "core_index_to_uniprot_interval"
    else:
        canonical_numbers = {
            position: residue_order[position - 1][0] for position in observed_by_position
        }
        numbering_mode = "deposited_author_numbering"
    remove: list[int] = []
    for atom in template.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            remove.append(atom.GetIdx())
            continue
        position = info.GetResidueNumber()
        if position not in observed_by_position or info.GetName().strip() not in observed_by_position[position]:
            remove.append(atom.GetIdx())
    editable = Chem.RWMol(template)
    for atom_index in sorted(remove, reverse=True):
        editable.RemoveAtom(atom_index)
    molecule = editable.GetMol()
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    chain = protein[0].chain
    for atom in molecule.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            raise ValueError("protein_template_atom_missing_residue_info")
        position = info.GetResidueNumber()
        observed = observed_by_position[position][info.GetName().strip()]
        conformer.SetAtomPosition(atom.GetIdx(), tuple(float(value) for value in observed.coord))
        info.SetResidueNumber(canonical_numbers[position])
        info.SetChainId(chain)
    molecule.RemoveAllConformers()
    molecule.AddConformer(conformer)
    sanitize_status = Chem.SanitizeMol(molecule, catchErrors=True)
    if int(sanitize_status) != 0:
        raise ValueError(f"protein_sequence_template_sanitize_failed:{int(sanitize_status)}")
    return molecule, numbering_mode, molecule.GetNumAtoms()


def prolif_contacts(ligand: Chem.Mol, protein: Chem.Mol) -> dict[str, Any]:
    import prolif as plf

    prepared_ligand = Chem.AddHs(Chem.Mol(ligand), addCoords=True)
    fingerprint = plf.Fingerprint(count=True, parameters={"VdWContact": {"vdwradii": {"Ir": 2.00, "Ru": 2.05}}})
    fingerprint.run_from_iterable(
        [plf.Molecule(prepared_ligand)], plf.Molecule(protein), n_jobs=1, progress=False
    )
    frame = fingerprint.to_dataframe()
    counts = {name: 0 for name in INTERACTIONS}
    contacts: list[dict[str, Any]] = []
    residues: set[str] = set()
    if len(frame):
        for column, value in frame.iloc[0].items():
            try:
                count = int(value)
            except Exception:
                continue
            if count <= 0:
                continue
            if isinstance(column, tuple) and len(column) >= 3:
                protein_residue = str(column[1])
                interaction = str(column[2])
            else:
                parts = str(column).split("|")
                protein_residue = parts[-2] if len(parts) >= 2 else str(column)
                interaction = parts[-1]
            counts[interaction] = counts.get(interaction, 0) + count
            residues.add(protein_residue)
            contacts.append(
                {"protein_residue": protein_residue, "interaction": interaction, "count": count}
            )
    result: dict[str, Any] = {
        "prolif_status": "ok",
        "prolif_protocol": PROLIF_PROTOCOL,
        "prolif_total_interactions": sum(counts.values()),
        "prolif_contact_residue_count": len(residues),
        "prolif_contacts_json": json.dumps(
            sorted(contacts, key=lambda item: (item["protein_residue"], item["interaction"])),
            separators=(",", ":"),
        ),
    }
    for interaction in INTERACTIONS:
        result[f"prolif_count_{interaction.lower()}"] = counts.get(interaction, 0)
    return result


def pdb_atom_name(name: str, element: str) -> str:
    value = str(name).strip()[:4]
    if len(value) >= 4:
        return value
    return f" {value:<3}" if len(str(element).strip()) == 1 else f"{value:<4}"


def write_receptor_pdb(protein: list[Atom], path: Path) -> float:
    lines: list[str] = []
    max_rounding = 0.0
    for serial, atom in enumerate(protein, start=1):
        x, y, z = (float(value) for value in atom.coord)
        max_rounding = max(
            max_rounding,
            abs(x - round(x, 3)),
            abs(y - round(y, 3)),
            abs(z - round(z, 3)),
        )
        chain = atom.chain[0] if atom.chain and atom.chain[0].isalnum() else "A"
        name = pdb_atom_name(atom.name, atom.element)
        lines.append(
            f"ATOM  {serial:5d} {name} {atom.residue[:3]:>3} {chain}{atom.residue_number:4d}{atom.insertion_code[:1]:1}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {atom.element:>2}"
        )
    lines.extend(["TER", "END"])
    path.write_text("\n".join(lines) + "\n")
    return max_rounding


def posebusters_report(ligand: Chem.Mol, protein: list[Atom], temporary: Path) -> dict[str, Any]:
    receptor = temporary / "receptor.pdb"
    rounding = write_receptor_pdb(protein, receptor)
    with open(os.devnull, "w") as sink, contextlib.redirect_stderr(sink):
        report = _BUSTER.bust(ligand, mol_cond=receptor, full_report=True)
    if report.empty:
        raise ValueError("posebusters_empty_report")
    raw = {str(key): normalize(value) for key, value in report.iloc[0].to_dict().items()}
    checks = {check: raw.get(check) for check in POSEBUSTERS_CORE_CHECKS}
    failed = [check for check, value in checks.items() if value is not True]
    result: dict[str, Any] = {
        "posebusters_status": "ok",
        "posebusters_protocol": POSEBUSTERS_PROTOCOL,
        "posebusters_overall_pass": not failed,
        "posebusters_n_checks": len(POSEBUSTERS_CORE_CHECKS),
        "posebusters_n_failed_checks": len(failed),
        "posebusters_failed_checks_json": json.dumps(failed, separators=(",", ":")),
        "posebusters_ligand_chemistry_valid": all(
            checks.get(name) is True
            for name in ["sanitization", "inchi_convertible", "all_atoms_connected", "no_radicals"]
        ),
        "posebusters_ligand_geometry_valid": all(
            checks.get(name) is True
            for name in [
                "bond_lengths", "bond_angles", "internal_steric_clash", "aromatic_ring_flatness",
                "non-aromatic_ring_non-flatness", "double_bond_flatness", "internal_energy",
            ]
        ),
        "posebusters_protein_ligand_clash_free": all(
            checks.get(name) is True
            for name in ["minimum_distance_to_protein", "volume_overlap_with_protein"]
        ),
        "posebusters_in_pocket": checks.get("protein-ligand_maximum_distance") is True,
        "posebusters_receptor_pdb_rounding_max_abs_delta_A": rounding,
    }
    for key, value in raw.items():
        result[f"posebusters_{key}"] = value
    return result


def empty_metric_fields(status: str) -> dict[str, Any]:
    result: dict[str, Any] = {
        "prolif_status": status,
        "prolif_protocol": PROLIF_PROTOCOL,
        "prolif_total_interactions": None,
        "prolif_contact_residue_count": None,
        "prolif_contacts_json": "[]",
        "posebusters_status": status,
        "posebusters_protocol": POSEBUSTERS_PROTOCOL,
        "posebusters_overall_pass": None,
        "posebusters_n_checks": len(POSEBUSTERS_CORE_CHECKS),
        "posebusters_n_failed_checks": None,
        "posebusters_failed_checks_json": "[]",
    }
    for interaction in INTERACTIONS:
        result[f"prolif_count_{interaction.lower()}"] = None
    return result


def compute_one(row: dict[str, Any]) -> dict[str, Any]:
    started = time.time()
    base: dict[str, Any] = {
        "structure_id": str(row["structure_id"]),
        "qc_protocol": QC_PROTOCOL,
        "qc_started_at": now(),
    }
    if not row.get("ligand_selected_ccd") or pd.isna(row.get("ligand_selected_ccd")):
        base.update(
            {
                "qc_status": "not_applicable_no_selected_ligand",
                "component_extraction_status": "not_applicable_no_selected_ligand",
                "ligand_topology_status": "not_applicable_no_selected_ligand",
                "protein_heavy_atoms": None,
                "coordinate_heavy_atoms": None,
                **empty_metric_fields("not_applicable_no_selected_ligand"),
            }
        )
        base["qc_seconds"] = time.time() - started
        return base
    try:
        if _RELEASE is None:
            raise RuntimeError("worker_not_initialized")
        coordinate = _RELEASE / str(row["coordinate_path"])
        atoms = parse_atoms(coordinate)
        protein, ligand_atoms = select_components(row, atoms)
        base.update(
            {
                "component_extraction_status": "ok",
                "coordinate_heavy_atoms": len(atoms),
                "protein_heavy_atoms": len(protein),
            }
        )
        try:
            ligand, topology = ligand_molecule(str(row["ligand_smiles"]), ligand_atoms)
        except Exception as primary_error:
            if str(row["source_kind"]) == "experimental" and str(row.get("experimental_source_format")) == "PDB":
                adjacency = deposited_pdb_conect_adjacency(row, ligand_atoms)
                ligand, topology = ligand_molecule(str(row["ligand_smiles"]), ligand_atoms, adjacency)
                topology["ligand_atom_mapping_mode"] = "deposited_PDB_CONECT_graph_isomorphism"
                topology["ligand_distance_mapping_error"] = f"{type(primary_error).__name__}:{str(primary_error)[:1000]}"
            else:
                raise
        base.update(topology)
    except Exception as exc:
        base.update(
            {
                "qc_status": "failed_component_or_topology",
                "qc_error": f"{type(exc).__name__}:{str(exc)[:1000]}",
                "ligand_topology_status": "failed",
                **empty_metric_fields("not_calculated_component_or_topology_failure"),
            }
        )
        base["qc_seconds"] = time.time() - started
        return base

    try:
        protein_molecule, numbering_mode, protein_rdkit_atoms = protein_rdkit_molecule(
            protein, str(row["source_kind"]), str(row["gene"])
        )
        base["prolif_protein_numbering_mode"] = numbering_mode
        base["prolif_protein_rdkit_atoms"] = protein_rdkit_atoms
        base.update(prolif_contacts(ligand, protein_molecule))
    except Exception as exc:
        base.update(
            {
                "prolif_status": "failed",
                "prolif_error": f"{type(exc).__name__}:{str(exc)[:1000]}",
                "prolif_protocol": PROLIF_PROTOCOL,
            }
        )

    try:
        with tempfile.TemporaryDirectory(prefix="cyp_release_qc_") as temporary:
            base.update(posebusters_report(ligand, protein, Path(temporary)))
    except Exception as exc:
        base.update(
            {
                "posebusters_status": "failed",
                "posebusters_error": f"{type(exc).__name__}:{str(exc)[:1000]}",
                "posebusters_protocol": POSEBUSTERS_PROTOCOL,
            }
        )
    prolif_ok = base.get("prolif_status") == "ok"
    posebusters_ok = base.get("posebusters_status") == "ok"
    base["qc_status"] = "complete" if prolif_ok and posebusters_ok else "partial_metric_failure"
    base["qc_seconds"] = time.time() - started
    return {key: normalize(value) for key, value in base.items()}


def valid_shard(path: Path, expected_ids: list[str]) -> bool:
    try:
        frame = pd.read_parquet(path, columns=["structure_id", "qc_protocol"])
        return (
            len(frame) == len(expected_ids)
            and frame.structure_id.tolist() == expected_ids
            and frame.qc_protocol.eq(QC_PROTOCOL).all()
        )
    except Exception:
        return False


def merge_shards(base: pd.DataFrame, shard_paths: list[Path], output: Path) -> pd.DataFrame:
    frames = [pd.read_parquet(path) for path in shard_paths]
    qc = pd.concat(frames, ignore_index=True, sort=False)
    if len(qc) != len(base) or qc.structure_id.nunique() != len(base):
        raise ValueError(f"qc_universe_gate_failed:rows={len(qc)}:ids={qc.structure_id.nunique()}")
    qc = base[["structure_id"]].merge(qc, on="structure_id", validate="one_to_one")
    atomic_parquet(qc, output)
    return qc


def json_value_counts(series: pd.Series) -> dict[str, int]:
    """Return stable, JSON-safe counts including missing values."""
    counts = series.value_counts(dropna=False)
    return {
        ("null" if pd.isna(key) else str(key).lower() if isinstance(key, (bool, np.bool_)) else str(key)): int(value)
        for key, value in counts.items()
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, default=DEFAULT_WORK / "base_structures.parquet")
    parser.add_argument("--release", type=Path, default=DEFAULT_RELEASE)
    parser.add_argument("--work", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--shard-size", type=int, default=100)
    parser.add_argument("--limit", type=int, help="Smoke-test N rows without writing production shards")
    args = parser.parse_args()
    base = pd.read_parquet(args.base)
    if len(base) != 16739 or base.structure_id.nunique() != 16739:
        raise ValueError("base_release_universe_gate_failed")

    if args.limit is not None:
        selected = pd.concat(
            [
                base[base.ligand_selected_ccd.notna()].head(max(1, args.limit - 1)),
                base[base.ligand_selected_ccd.isna()].head(1),
            ],
            ignore_index=True,
        ).head(args.limit)
        initialize_worker(str(args.release.resolve()))
        smoke = pd.DataFrame([compute_one(row) for row in selected.to_dict("records")])
        print(smoke.to_json(orient="records", indent=2))
        return

    shard_dir = args.work / "qc" / "shards"
    shard_dir.mkdir(parents=True, exist_ok=True)
    partitions: list[tuple[int, pd.DataFrame, Path]] = []
    for start in range(0, len(base), args.shard_size):
        number = start // args.shard_size
        shard = base.iloc[start : start + args.shard_size].copy()
        path = shard_dir / f"shard_{number:04d}.parquet"
        partitions.append((number, shard, path))

    context = mp.get_context("spawn")
    with ProcessPoolExecutor(
        max_workers=args.workers,
        mp_context=context,
        initializer=initialize_worker,
        initargs=(str(args.release.resolve()),),
    ) as pool:
        for number, shard, path in partitions:
            expected_ids = shard.structure_id.astype(str).tolist()
            if valid_shard(path, expected_ids):
                print(json.dumps({"event": "qc_shard_skip", "shard": number, "rows": len(shard)}), flush=True)
                continue
            records = shard.where(pd.notna(shard), None).to_dict("records")
            results = list(pool.map(compute_one, records, chunksize=20))
            result = pd.DataFrame(results)
            result = shard[["structure_id"]].merge(result, on="structure_id", validate="one_to_one")
            atomic_parquet(result, path)
            counts = result.qc_status.value_counts(dropna=False).to_dict()
            print(
                json.dumps(
                    {"event": "qc_shard_complete", "shard": number, "rows": len(result), "status_counts": counts}
                ),
                flush=True,
            )

    qc = merge_shards(base, [partition[2] for partition in partitions], args.work / "qc" / "structure_qc.parquet")
    summary = {
        "created_at": now(),
        "protocol": QC_PROTOCOL,
        "rows": len(qc),
        "status_counts": json_value_counts(qc.qc_status),
        "prolif_status_counts": json_value_counts(qc.prolif_status),
        "posebusters_status_counts": json_value_counts(qc.posebusters_status),
        "posebusters_overall_pass_counts": json_value_counts(qc.posebusters_overall_pass),
        "atom_mapping_mode_counts": json_value_counts(qc.ligand_atom_mapping_mode),
        "total_qc_seconds": float(qc.qc_seconds.sum()),
        "workers": args.workers,
        "shard_size": args.shard_size,
    }
    atomic_json(summary, args.work / "qc" / "summary.json")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
