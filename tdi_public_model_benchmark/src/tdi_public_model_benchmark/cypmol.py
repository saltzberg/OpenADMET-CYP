"""Repaired, resumable CYPMol substrate/inhibitor inference.

Run this module only in the isolated ``cypmol`` environment. It loads the
published full-model checkpoints strictly, uses the authors' deterministic
RDKit conformer recipe, and scores CYP2D6/CYP3A4 with fixed endpoint-specific
protein contexts.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.metadata
import json
import os
import pickle
import platform
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

BENCHMARK = Path(__file__).resolve().parents[2]
VENDOR = BENCHMARK / "vendor"
for source_root in [VENDOR / "Uni-Core-0.0.1", VENDOR / "Uni-Mol-v0.1"]:
    sys.path.insert(0, str(source_root))

import torch
import torch.nn as nn
import esm
from esm.data import Alphabet
from esm.model.esm2 import ESM2
from unicore.data import Dictionary
from unicore.modules import init_bert_params
from unimol.models.transformer_encoder_with_pair import TransformerEncoderWithPair
from unimol.models.unimol import GaussianLayer, NonLinearHead

HF_REVISION = "3810a93e04d1dbab4cde73792a1e3830b1f17e41"
GIT_REVISION = "0c6657dc882f23293ff9a43977c8f31818f8a6f8"
ENDPOINTS = ["CYP2D6", "CYP3A4"]
CHECKPOINTS = {
    "substrate": [
        "new_sequence_srs_pocket_sub1.pt", "new_sequence_srs_pocket_sub2.pt",
        "new_sequence_srs_pocket_sub3.pt", "new_sequence_srs_pocket_sub5.pt",
        "new_sequence_srs_pocket_sub6.pt", "new_sequence_srs_pocket_sub10.pt",
        "new_sequence_srs_pocket_sub11.pt", "new_sequence_srs_pocket_sub13.pt",
        "new_sequence_srs_pocket_sub15.pt", "new_sequence_srs_pocket_sub16.pt",
    ],
    "inhibitor": [
        "inhi_seq_pocket_sub0.pt", "inhi_seq_pocket_sub1.pt",
        "inhi_seq_pocket_sub2.pt", "inhi_seq_pocket_sub3.pt",
        "inhi_seq_pocket_sub4.pt", "inhi_seq_pocket_sub7.pt",
        "inhi_seq_pocket_sub9.pt", "inhi_seq_pocket_sub10.pt",
        "inhi_seq_pocket_sub12.pt", "inhi_seq_pocket_sub15.pt",
    ],
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class UniMolModel(nn.Module):
    def __init__(self, dictionary: Dictionary, initialize: bool = True):
        super().__init__()
        self.padding_idx = dictionary.pad()
        self.embed_tokens = nn.Embedding(len(dictionary), 512, self.padding_idx)
        self.encoder = TransformerEncoderWithPair(
            encoder_layers=15,
            embed_dim=512,
            ffn_embed_dim=2048,
            attention_heads=64,
            emb_dropout=0.1,
            dropout=0.1,
            attention_dropout=0.1,
            activation_dropout=0.0,
            max_seq_len=512,
            activation_fn="gelu",
            no_final_head_layer_norm=True,
        )
        n_edge_type = len(dictionary) * len(dictionary)
        self.gbf_proj = NonLinearHead(128, 64, "gelu")
        self.gbf = GaussianLayer(128, n_edge_type)
        if initialize:
            self.apply(init_bert_params)

    def forward(self, src_tokens: torch.Tensor, src_distance: torch.Tensor,
                src_edge_type: torch.Tensor) -> torch.Tensor:
        padding_mask = src_tokens.eq(self.padding_idx)
        if not padding_mask.any():
            padding_mask = None
        x = self.embed_tokens(src_tokens)
        n_node = src_distance.size(-1)
        graph_attn_bias = self.gbf_proj(self.gbf(src_distance, src_edge_type))
        graph_attn_bias = graph_attn_bias.permute(0, 3, 1, 2).contiguous()
        graph_attn_bias = graph_attn_bias.view(-1, n_node, n_node)
        encoder_rep, _, _, _, _ = self.encoder(
            x, padding_mask=padding_mask, attn_mask=graph_attn_bias
        )
        return encoder_rep


class CrossAttentionLayer(nn.Module):
    def __init__(self):
        super().__init__()
        self.layer_normalization_cross_attention_1 = nn.LayerNorm(1280)
        self.cross_attention = nn.MultiheadAttention(
            embed_dim=512, num_heads=8, kdim=1280, vdim=1280, batch_first=True
        )
        self.layer_normalization_cross_attention_2 = nn.LayerNorm(512)
        self.feed_forward_cross_attention = nn.Sequential(
            nn.Linear(512, 512), nn.GELU(), nn.Linear(512, 512)
        )

    def forward(self, x: torch.Tensor, context: torch.Tensor) -> torch.Tensor:
        context = self.layer_normalization_cross_attention_1(context)
        y, _ = self.cross_attention(x, context, context, key_padding_mask=None)
        old = y
        y = self.layer_normalization_cross_attention_2(y)
        return self.feed_forward_cross_attention(y) + old


class SelfAttentionLayer(nn.Module):
    def __init__(self):
        super().__init__()
        self.layer_normalization_self_attention_1 = nn.LayerNorm(512)
        self.self_attention = nn.MultiheadAttention(
            embed_dim=512, num_heads=8, kdim=512, vdim=512, batch_first=True
        )
        self.layer_normalization_self_attention_2 = nn.LayerNorm(512)
        self.feed_forward_self_attention = nn.Sequential(
            nn.Linear(512, 512), nn.GELU(), nn.Linear(512, 512)
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        old = x
        x = self.layer_normalization_self_attention_1(x)
        x, _ = self.self_attention(x, x, x)
        x = x + old
        old = x
        x = self.layer_normalization_self_attention_2(x)
        return self.feed_forward_self_attention(x) + old


class CYPMolClassifier(nn.Module):
    def __init__(self, task: str, dictionary: Dictionary, initialize: bool = True):
        super().__init__()
        self.task = task
        self.molecule_encoder = UniMolModel(dictionary, initialize=initialize)
        alphabet = Alphabet.from_architecture("ESM-1b")
        self.protein_encoder = ESM2(
            num_layers=33,
            embed_dim=1280,
            attention_heads=20,
            alphabet=alphabet,
            token_dropout=True,
        )
        self.alphabet = alphabet
        self.batch_converter = alphabet.get_batch_converter(truncation_seq_length=2048)
        self.transformer_layer_cross_attention = CrossAttentionLayer()
        self.transformer_layer_self_attention = SelfAttentionLayer()
        if task == "substrate":
            self.mlp = nn.Sequential(
                nn.Linear(512, 128), nn.ReLU(), nn.Linear(128, 32), nn.ReLU(),
                nn.Linear(32, 32), nn.ReLU(), nn.Linear(32, 2),
            )
        elif task == "inhibitor":
            self.mlp = nn.Sequential(nn.Linear(512, 128), nn.ReLU(), nn.Linear(128, 2))
        else:
            raise ValueError(task)

    @torch.inference_mode()
    def protein_context(self, sequence: str, pocket: list[int],
                        srs_intervals: list[int] | None, device: torch.device) -> torch.Tensor:
        _, _, tokens = self.batch_converter([("", sequence)])
        tokens = tokens.to(device)
        output = self.protein_encoder(tokens, repr_layers=[33], return_contacts=False)
        embedding = output["representations"][33]
        pieces = [embedding]
        if self.task == "substrate":
            assert srs_intervals is not None and len(srs_intervals) % 2 == 0
            selected = []
            for index in range(0, len(srs_intervals), 2):
                begin = srs_intervals[index] + 1
                end = srs_intervals[index + 1] + 2
                selected.extend(embedding[0, begin:end])
            pieces.append(torch.stack(selected).unsqueeze(0))
        pocket_embedding = torch.stack([embedding[0, site + 1] for site in pocket]).unsqueeze(0)
        pieces.append(pocket_embedding)
        return torch.cat(pieces, dim=1)

    @torch.inference_mode()
    def score(self, batch: dict[str, torch.Tensor], context: torch.Tensor) -> torch.Tensor:
        molecule = self.molecule_encoder(
            batch["src_tokens"], batch["src_distance"], batch["src_edge_type"]
        )
        expanded = context.expand(molecule.size(0), -1, -1)
        x = self.transformer_layer_cross_attention(molecule, expanded)
        x = self.transformer_layer_self_attention(x)
        return torch.softmax(self.mlp(x[:, 0, :]), dim=1)[:, 1]


def load_dictionary(path: Path) -> Dictionary:
    dictionary = Dictionary.load(str(path))
    dictionary.add_symbol("[MASK]", is_special=True)
    if len(dictionary) != 31:
        raise ValueError(f"Unexpected CYPMol dictionary size: {len(dictionary)}")
    return dictionary


def load_model(task: str, dictionary: Dictionary, checkpoint: Path,
               device: torch.device) -> CYPMolClassifier:
    with torch.device("meta"):
        model = CYPMolClassifier(task, dictionary, initialize=False)
    state = torch.load(
        checkpoint, map_location="cpu", weights_only=True, mmap=True
    )
    expected = set(model.state_dict())
    observed = set(state)
    if expected != observed:
        raise ValueError(
            f"Checkpoint architecture mismatch: missing={sorted(expected-observed)[:20]}, "
            f"unexpected={sorted(observed-expected)[:20]}"
        )
    mismatches = {
        key: (tuple(model.state_dict()[key].shape), tuple(state[key].shape))
        for key in expected
        if tuple(model.state_dict()[key].shape) != tuple(state[key].shape)
    }
    if mismatches:
        raise ValueError(f"Checkpoint shape mismatch: {dict(list(mismatches.items())[:10])}")
    model.load_state_dict(state, strict=True, assign=True)
    model = model.to(device)
    model.eval()
    return model


def parse_tuple(value: str) -> list[int]:
    text = str(value).strip().strip("()")
    return [] if not text else [int(item.strip()) for item in text.split(",")]


def load_endpoint_contexts(task: str, source_repo: Path) -> dict[str, dict[str, Any]]:
    if task == "substrate":
        frame = pd.read_csv(source_repo / "Substrate/substrate_exammple.csv")
        sequence_column = "sequence"
    else:
        with zipfile.ZipFile(source_repo / "Inhibitor/inhibitor.zip") as archive:
            with archive.open("inhibitor.csv") as handle:
                frame = pd.read_csv(handle, encoding="gb18030")
        sequence_column = "protein_seq"
    frame["CYP"] = frame["CYP"].astype(str).str.upper()
    contexts = {}
    for endpoint in ENDPOINTS:
        rows = frame[frame["CYP"] == endpoint.removeprefix("CYP")]
        if rows.empty:
            raise ValueError(f"No {endpoint} context in {task} example")
        fields = [sequence_column, "pocket_site"] + (["protein_SRS"] if task == "substrate" else [])
        for field in fields:
            if rows[field].nunique() != 1:
                raise ValueError(f"Conflicting {task} {endpoint} {field} values")
        row = rows.iloc[0]
        contexts[endpoint] = {
            "sequence": str(row[sequence_column]),
            "pocket": parse_tuple(row["pocket_site"]),
            "srs": parse_tuple(row["protein_SRS"]) if task == "substrate" else None,
        }
    return contexts


def prepare_conformer(smiles: str) -> dict[str, Any]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("rdkit_parse_failure")
    mol = AllChem.AddHs(mol)
    result = AllChem.EmbedMolecule(
        mol, randomSeed=42, useRandomCoords=True, maxAttempts=1000
    )
    if result != 0:
        raise ValueError(f"embed_failure_{result}")
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception as error:
        raise ValueError(f"mmff_exception_{type(error).__name__}") from error
    atoms = np.asarray([atom.GetSymbol() for atom in mol.GetAtoms()])
    coordinates = mol.GetConformer().GetPositions().astype(np.float32)
    keep = atoms != "H"
    atoms = atoms[keep]
    coordinates = coordinates[keep]
    if len(atoms) > 256:
        seed = int(hash((1, None, 0)) % 1e6)
        rng = np.random.RandomState(seed)
        selected = rng.choice(len(atoms), 256, replace=False)
        atoms = atoms[selected]
        coordinates = coordinates[selected]
    coordinates = (coordinates - coordinates.mean(axis=0)).astype(np.float32)
    return {"atoms": atoms.tolist(), "coordinates": coordinates}


def prepare_cache(input_csv: Path, cache_path: Path, failure_path: Path,
                  limit: int | None = None) -> tuple[pd.DataFrame, dict[str, Any]]:
    source = pd.read_csv(input_csv)
    if limit is not None:
        source = source.head(limit).copy()
    source = source[["Molecule_Name", "SMILES"]].rename(columns={"SMILES": "model_smiles"})
    input_hash = sha256_file(input_csv)
    metadata_path = cache_path.with_suffix(".manifest.json")
    if cache_path.is_file() and metadata_path.is_file():
        metadata = json.loads(metadata_path.read_text())
        if metadata.get("input_sha256") == input_hash and metadata.get("limit") == limit:
            with cache_path.open("rb") as handle:
                return source, pickle.load(handle)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache = {}
    failures = []
    for position, row in enumerate(source.itertuples(index=False), start=1):
        try:
            cache[row.Molecule_Name] = prepare_conformer(row.model_smiles)
        except Exception as error:
            failures.append({
                "Molecule_Name": row.Molecule_Name,
                "model_smiles": row.model_smiles,
                "failure": str(error),
            })
        if position == 1 or position % 250 == 0:
            print(f"conformers {position}/{len(source)} failures={len(failures)}", flush=True)
    with cache_path.open("wb") as handle:
        pickle.dump(cache, handle, protocol=pickle.HIGHEST_PROTOCOL)
    pd.DataFrame(failures, columns=["Molecule_Name", "model_smiles", "failure"]).to_csv(
        failure_path, index=False
    )
    metadata_path.write_text(json.dumps({
        "input_sha256": input_hash,
        "limit": limit,
        "rows": len(source),
        "valid": len(cache),
        "failures": len(failures),
        "cache_sha256": sha256_file(cache_path),
    }, indent=2) + "\n")
    return source, cache


def molecule_batch(records: list[tuple[str, dict[str, Any]]], dictionary: Dictionary,
                   device: torch.device) -> dict[str, torch.Tensor]:
    token_rows = []
    distance_rows = []
    edge_rows = []
    for _, record in records:
        atom_tokens = [dictionary.index(atom) for atom in record["atoms"]]
        tokens = np.asarray([dictionary.bos(), *atom_tokens, dictionary.eos()], dtype=np.int64)
        coordinates = np.concatenate(
            [np.zeros((1, 3), dtype=np.float32), record["coordinates"],
             np.zeros((1, 3), dtype=np.float32)], axis=0
        )
        distances = np.linalg.norm(coordinates[:, None, :] - coordinates[None, :, :], axis=-1)
        edge_types = tokens[:, None] * len(dictionary) + tokens[None, :]
        token_rows.append(tokens)
        distance_rows.append(distances.astype(np.float32))
        edge_rows.append(edge_types.astype(np.int64))
    length = max(len(row) for row in token_rows)
    tokens = np.full((len(records), length), dictionary.pad(), dtype=np.int64)
    distances = np.zeros((len(records), length, length), dtype=np.float32)
    edges = np.zeros((len(records), length, length), dtype=np.int64)
    for index, (token_row, distance_row, edge_row) in enumerate(
        zip(token_rows, distance_rows, edge_rows)
    ):
        n = len(token_row)
        tokens[index, :n] = token_row
        distances[index, :n, :n] = distance_row
        edges[index, :n, :n] = edge_row
    return {
        "src_tokens": torch.from_numpy(tokens).to(device),
        "src_distance": torch.from_numpy(distances).to(device),
        "src_edge_type": torch.from_numpy(edges).to(device),
    }


def package_versions() -> dict[str, str]:
    result = {}
    for name in ["torch", "fair-esm", "rdkit", "numpy", "pandas"]:
        try:
            result[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            result[name] = "not-installed"
    return result


def run_task(task: str, input_csv: Path, source_repo: Path, assets_dir: Path,
             dictionary_path: Path, output_csv: Path, manifest_path: Path,
             failure_path: Path, work_dir: Path, batch_size: int, limit: int | None,
             max_checkpoints: int | None) -> dict[str, Any]:
    if task not in CHECKPOINTS:
        raise ValueError(task)
    if not torch.cuda.is_available() or torch.cuda.get_device_capability()[0] < 12:
        raise RuntimeError("CYPMol requires the verified Blackwell CUDA runtime")
    device = torch.device("cuda")
    torch.backends.cuda.matmul.allow_tf32 = False
    dictionary = load_dictionary(dictionary_path)
    contexts = load_endpoint_contexts(task, source_repo)
    source, cache = prepare_cache(
        input_csv, work_dir / "conformers.pkl", failure_path, limit=limit
    )
    valid_ids = [identifier for identifier in source["Molecule_Name"] if identifier in cache]
    checkpoints = CHECKPOINTS[task][:max_checkpoints]
    task_asset_dir = assets_dir / task.capitalize()
    missing = [name for name in checkpoints if not (task_asset_dir / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing CYPMol {task} checkpoints: {missing}")
    partial_dir = work_dir / task
    partial_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_records = []
    score_frames = []
    for checkpoint_index, name in enumerate(checkpoints, start=1):
        checkpoint = task_asset_dir / name
        digest = sha256_file(checkpoint)
        partial_csv = partial_dir / f"{checkpoint.stem}.csv"
        partial_manifest = partial_csv.with_suffix(".json")
        if partial_csv.is_file() and partial_manifest.is_file():
            metadata = json.loads(partial_manifest.read_text())
            if metadata.get("checkpoint_sha256") == digest and metadata.get("ids") == len(valid_ids):
                frame = pd.read_csv(partial_csv)
                score_frames.append(frame)
                checkpoint_records.append(metadata)
                print(f"reuse {task} checkpoint {checkpoint_index}/{len(checkpoints)} {name}", flush=True)
                continue
        print(f"load {task} checkpoint {checkpoint_index}/{len(checkpoints)} {name}", flush=True)
        model = load_model(task, dictionary, checkpoint, device)
        endpoint_context = {
            endpoint: model.protein_context(
                contexts[endpoint]["sequence"], contexts[endpoint]["pocket"],
                contexts[endpoint]["srs"], device
            )
            for endpoint in ENDPOINTS
        }
        rows = []
        for start in range(0, len(valid_ids), batch_size):
            ids = valid_ids[start:start + batch_size]
            records = [(identifier, cache[identifier]) for identifier in ids]
            batch = molecule_batch(records, dictionary, device)
            scores = {
                endpoint: model.score(batch, endpoint_context[endpoint]).float().cpu().numpy()
                for endpoint in ENDPOINTS
            }
            for row_index, identifier in enumerate(ids):
                rows.append({
                    "Molecule_Name": identifier,
                    **{f"{endpoint}_{task}_{checkpoint.stem}": float(scores[endpoint][row_index])
                       for endpoint in ENDPOINTS},
                })
            if start == 0 or (start // batch_size) % 100 == 0:
                print(f"{task} {name} rows={min(start+batch_size,len(valid_ids))}/{len(valid_ids)}", flush=True)
        frame = pd.DataFrame(rows)
        frame.to_csv(partial_csv, index=False)
        metadata = {
            "checkpoint": name,
            "checkpoint_sha256": digest,
            "checkpoint_size": checkpoint.stat().st_size,
            "ids": len(valid_ids),
            "output_sha256": sha256_file(partial_csv),
        }
        partial_manifest.write_text(json.dumps(metadata, indent=2) + "\n")
        checkpoint_records.append(metadata)
        score_frames.append(frame)
        del model, endpoint_context, batch
        gc.collect()
        torch.cuda.empty_cache()
    combined = source.copy()
    task_columns = {endpoint: [] for endpoint in ENDPOINTS}
    for frame in score_frames:
        combined = combined.merge(frame, on="Molecule_Name", how="left", validate="one_to_one")
        for endpoint in ENDPOINTS:
            task_columns[endpoint].extend(
                column for column in frame.columns if column.startswith(f"{endpoint}_{task}_")
            )
    for endpoint in ENDPOINTS:
        values = combined[task_columns[endpoint]]
        combined[f"{endpoint}_{task}_score"] = values.mean(axis=1)
        combined[f"{endpoint}_{task}_variance"] = values.var(axis=1, ddof=0)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_csv, index=False)
    score_columns = [f"{endpoint}_{task}_score" for endpoint in ENDPOINTS]
    scored = combined[score_columns].notna().all(axis=1)
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "model": f"CYPMol {task} ensemble",
        "task": task,
        "git_revision": GIT_REVISION,
        "hf_revision": HF_REVISION,
        "score_semantics": f"mean class-1 {task} score across {len(checkpoints)} published checkpoints",
        "input_csv": str(input_csv.resolve()),
        "input_sha256": sha256_file(input_csv),
        "input_rows": len(source),
        "scored_rows": int(scored.sum()),
        "failure_rows": int((~scored).sum()),
        "representation": "original supplied SMILES; authors' AddHs/ETKDG randomSeed=42/MMFF/heavy-atom Uni-Mol preprocessing",
        "endpoint_context": contexts,
        "batch_size": batch_size,
        "checkpoint_count": len(checkpoints),
        "checkpoints": checkpoint_records,
        "runtime": {
            "python": platform.python_version(),
            "packages": package_versions(),
            "cuda": torch.version.cuda,
            "gpu": torch.cuda.get_device_name(),
        },
        "outputs": {
            "predictions": {"path": str(output_csv), "sha256": sha256_file(output_csv)},
            "failures": {"path": str(failure_path), "sha256": sha256_file(failure_path)},
            "conformer_cache": {
                "path": str(work_dir / "conformers.pkl"),
                "sha256": sha256_file(work_dir / "conformers.pkl"),
            },
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", choices=sorted(CHECKPOINTS), required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--source-repo", type=Path, required=True)
    parser.add_argument("--assets-dir", type=Path, required=True)
    parser.add_argument("--dictionary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--failures", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=6)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--max-checkpoints", type=int)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    result = run_task(
        task=args.task,
        input_csv=args.input,
        source_repo=args.source_repo,
        assets_dir=args.assets_dir,
        dictionary_path=args.dictionary,
        output_csv=args.output,
        manifest_path=args.manifest,
        failure_path=args.failures,
        work_dir=args.work_dir,
        batch_size=args.batch_size,
        limit=args.limit,
        max_checkpoints=args.max_checkpoints,
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
