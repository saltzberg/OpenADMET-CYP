"""Resumable DeepP450 inference using complete published checkpoints.

This adapter reuses the verified CYPMol Uni-Mol/ESM-2 runtime because the two
public projects share the same encoder/fusion architecture. DeepP450 differs by
attending to the full CYP sequence without CYPMol pocket/SRS extensions.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch

from .cypmol import (
    CYPMolClassifier,
    ENDPOINTS,
    load_dictionary,
    load_endpoint_contexts,
    molecule_batch,
    package_versions,
    prepare_cache,
    sha256_file,
)

GIT_REVISION = "dea269a228177491aa02b02f5f3c396f415778ac"
CHECKPOINT_NAMES = [f"{index}.pt" for index in range(1, 21)]


def load_model(dictionary, checkpoint: Path, device: torch.device) -> CYPMolClassifier:
    """Build the shared architecture and require exact DeepP450 state parity."""
    with torch.device("meta"):
        model = CYPMolClassifier("substrate", dictionary, initialize=False)
    state = torch.load(checkpoint, map_location="cpu", weights_only=True, mmap=True)
    expected_state = model.state_dict()
    if set(expected_state) != set(state):
        raise ValueError(
            f"Checkpoint architecture mismatch: missing={sorted(set(expected_state)-set(state))[:20]}, "
            f"unexpected={sorted(set(state)-set(expected_state))[:20]}"
        )
    mismatches = {
        key: (tuple(expected_state[key].shape), tuple(state[key].shape))
        for key in expected_state
        if tuple(expected_state[key].shape) != tuple(state[key].shape)
    }
    if mismatches:
        raise ValueError(f"Checkpoint shape mismatch: {dict(list(mismatches.items())[:10])}")
    model.load_state_dict(state, strict=True, assign=True)
    model = model.to(device)
    model.eval()
    return model


@torch.inference_mode()
def sequence_context(model: CYPMolClassifier, sequence: str,
                     device: torch.device) -> torch.Tensor:
    _, _, tokens = model.batch_converter([("", sequence)])
    tokens = tokens.to(device)
    output = model.protein_encoder(tokens, repr_layers=[33], return_contacts=False)
    return output["representations"][33]


@torch.inference_mode()
def score(model: CYPMolClassifier, batch: dict[str, torch.Tensor],
          context: torch.Tensor) -> torch.Tensor:
    molecule = model.molecule_encoder(
        batch["src_tokens"], batch["src_distance"], batch["src_edge_type"]
    )
    expanded = context.expand(molecule.size(0), -1, -1)
    x = model.transformer_layer_cross_attention(molecule, expanded)
    x = model.transformer_layer_self_attention(x)
    return torch.softmax(model.mlp(x[:, 0, :]), dim=1)[:, 1]


def run_inference(input_csv: Path, checkpoint_dir: Path, dictionary_path: Path,
                  context_source_repo: Path, output_csv: Path, manifest_path: Path,
                  failure_path: Path, work_dir: Path, batch_size: int = 8,
                  limit: int | None = None, max_checkpoints: int | None = None) -> dict[str, Any]:
    if not torch.cuda.is_available() or torch.cuda.get_device_capability()[0] < 12:
        raise RuntimeError("DeepP450 requires the verified Blackwell CUDA runtime")
    device = torch.device("cuda")
    torch.backends.cuda.matmul.allow_tf32 = False
    dictionary = load_dictionary(dictionary_path)
    published_contexts = load_endpoint_contexts("substrate", context_source_repo)
    sequences = {endpoint: published_contexts[endpoint]["sequence"] for endpoint in ENDPOINTS}
    source, conformers = prepare_cache(
        input_csv, work_dir / "conformers.pkl", failure_path, limit=limit
    )
    valid_ids = [identifier for identifier in source["Molecule_Name"] if identifier in conformers]
    checkpoint_names = CHECKPOINT_NAMES[:max_checkpoints]
    missing = [name for name in checkpoint_names if not (checkpoint_dir / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing DeepP450 checkpoints: {missing}")

    partial_dir = work_dir / "deepp450"
    partial_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_records = []
    frames = []
    for checkpoint_index, name in enumerate(checkpoint_names, start=1):
        checkpoint = checkpoint_dir / name
        digest = sha256_file(checkpoint)
        partial_csv = partial_dir / f"{checkpoint.stem}.csv"
        partial_manifest = partial_csv.with_suffix(".json")
        if partial_csv.is_file() and partial_manifest.is_file():
            metadata = json.loads(partial_manifest.read_text())
            if metadata.get("checkpoint_sha256") == digest and metadata.get("ids") == len(valid_ids):
                frames.append(pd.read_csv(partial_csv))
                checkpoint_records.append(metadata)
                print(f"reuse checkpoint {checkpoint_index}/{len(checkpoint_names)} {name}", flush=True)
                continue
        print(f"load checkpoint {checkpoint_index}/{len(checkpoint_names)} {name}", flush=True)
        model = load_model(dictionary, checkpoint, device)
        contexts = {
            endpoint: sequence_context(model, sequences[endpoint], device)
            for endpoint in ENDPOINTS
        }
        rows = []
        for start in range(0, len(valid_ids), batch_size):
            ids = valid_ids[start:start + batch_size]
            batch = molecule_batch(
                [(identifier, conformers[identifier]) for identifier in ids], dictionary, device
            )
            scores = {
                endpoint: score(model, batch, contexts[endpoint]).float().cpu().numpy()
                for endpoint in ENDPOINTS
            }
            for row_index, identifier in enumerate(ids):
                rows.append({
                    "Molecule_Name": identifier,
                    **{f"{endpoint}_model_{checkpoint.stem}": float(scores[endpoint][row_index])
                       for endpoint in ENDPOINTS},
                })
            if start == 0 or (start // batch_size) % 100 == 0:
                print(f"{name} rows={min(start+batch_size,len(valid_ids))}/{len(valid_ids)}", flush=True)
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
        frames.append(frame)
        checkpoint_records.append(metadata)
        del model, contexts, batch
        gc.collect()
        torch.cuda.empty_cache()

    combined = source.copy()
    endpoint_columns = {endpoint: [] for endpoint in ENDPOINTS}
    for frame in frames:
        combined = combined.merge(frame, on="Molecule_Name", how="left", validate="one_to_one")
        for endpoint in ENDPOINTS:
            endpoint_columns[endpoint].extend(
                column for column in frame if column.startswith(f"{endpoint}_model_")
            )
    for endpoint in ENDPOINTS:
        values = combined[endpoint_columns[endpoint]]
        combined[f"{endpoint}_substrate_score"] = values.mean(axis=1)
        combined[f"{endpoint}_ensemble_variance"] = values.var(axis=1, ddof=0)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_csv, index=False)
    score_columns = [f"{endpoint}_substrate_score" for endpoint in ENDPOINTS]
    scored = combined[score_columns].notna().all(axis=1)
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "model": "DeepP450 substrate ensemble",
        "git_revision": GIT_REVISION,
        "score_semantics": f"mean class-1 substrate score across {len(checkpoint_names)} published checkpoints",
        "input_csv": str(input_csv.resolve()),
        "input_sha256": sha256_file(input_csv),
        "input_rows": len(source),
        "scored_rows": int(scored.sum()),
        "failure_rows": int((~scored).sum()),
        "representation": "original supplied SMILES; authors' AddHs/ETKDG randomSeed=42/MMFF/heavy-atom Uni-Mol preprocessing",
        "endpoint_sequences": sequences,
        "batch_size": batch_size,
        "checkpoint_count": len(checkpoint_names),
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
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--checkpoint-dir", type=Path, required=True)
    parser.add_argument("--dictionary", type=Path, required=True)
    parser.add_argument("--context-source-repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--failures", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--max-checkpoints", type=int)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    result = run_inference(
        input_csv=args.input,
        checkpoint_dir=args.checkpoint_dir,
        dictionary_path=args.dictionary,
        context_source_repo=args.context_source_repo,
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
