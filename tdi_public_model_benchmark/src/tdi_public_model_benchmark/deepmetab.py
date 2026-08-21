"""Pinned DeepMetab substrate inference with ID-preserving outputs.

DeepMetab's public checkpoints use ``readout.*`` state-dict keys while the
compatible ChemProp 1.x model calls the same layers ``ffn.*``. The upstream
``load_checkpoint`` helper silently skips those output-layer weights. This
adapter maps the names explicitly and then requires a strict full-state load.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

TASKS = [
    "CYP1A2",
    "CYP2A6",
    "CYP2B6",
    "CYP2C8",
    "CYP2C9",
    "CYP2C19",
    "CYP2D6",
    "CYP2E1",
    "CYP3A4",
]
EXPECTED_CHECKPOINT_COUNT = 5


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def map_checkpoint_state_dict(state_dict: dict[str, Any]) -> dict[str, Any]:
    """Map DeepMetab's output-head names to ChemProp 1.x names."""
    mapped: dict[str, Any] = {}
    for key, value in state_dict.items():
        new_key = "ffn." + key[len("readout.") :] if key.startswith("readout.") else key
        if new_key in mapped:
            raise ValueError(f"Checkpoint key collision after mapping: {new_key}")
        mapped[new_key] = value
    return mapped


def load_strict_checkpoint(path: Path, device: str = "cpu"):
    """Load a public DeepMetab checkpoint without silently randomizing its head."""
    import torch
    from chemprop.args import TrainArgs
    from chemprop.models import MoleculeModel

    target_device = torch.device(device)
    checkpoint = torch.load(str(path), map_location=target_device)
    args = TrainArgs()
    args.from_dict(vars(checkpoint["args"]), skip_unsettable=True)
    args.device = target_device
    model = MoleculeModel(args)
    mapped = map_checkpoint_state_dict(checkpoint["state_dict"])
    expected = model.state_dict()
    missing = sorted(set(expected) - set(mapped))
    unexpected = sorted(set(mapped) - set(expected))
    shape_mismatches = {
        key: {"checkpoint": list(mapped[key].shape), "model": list(expected[key].shape)}
        for key in sorted(set(mapped) & set(expected))
        if tuple(mapped[key].shape) != tuple(expected[key].shape)
    }
    if missing or unexpected or shape_mismatches:
        raise ValueError(
            "Checkpoint is not architecture-compatible: "
            f"missing={missing}, unexpected={unexpected}, shapes={shape_mismatches}"
        )
    model.load_state_dict(mapped, strict=True)
    model.to(target_device)
    model.eval()
    return model


def _prepare_chemprop_input(frame: pd.DataFrame, path: Path) -> None:
    model_input = pd.DataFrame({"Smiles": frame["model_smiles"]})
    for task in TASKS:
        model_input[task] = 0
    model_input.to_csv(path, index=False)


def _score_checkpoint(
    checkpoint: Path,
    model_input_path: Path,
    batch_size: int,
    device: str,
) -> tuple[list[str], np.ndarray]:
    from chemprop.args import TrainArgs
    from chemprop.data import MoleculeDataLoader, get_data
    from chemprop.train import predict

    arguments = [
        "--data_path",
        str(model_input_path),
        "--dataset_type",
        "classification",
        "--no_cuda",
    ]
    args = TrainArgs().parse_args(arguments)
    args.target_columns = TASKS
    args.smiles_columns = ["Smiles"]
    args.task_names = TASKS
    args.batch_size = batch_size
    args.num_workers = 0

    dataset = get_data(path=str(model_input_path), args=args)
    loader = MoleculeDataLoader(
        dataset=dataset,
        batch_size=batch_size,
        num_workers=0,
        shuffle=False,
    )
    model = load_strict_checkpoint(checkpoint, device=device)
    scores = np.asarray(predict(model=model, data_loader=loader), dtype=float)
    smiles = [entry[0] for entry in dataset.smiles()]
    if scores.shape != (len(smiles), len(TASKS)):
        raise ValueError(
            f"Unexpected prediction shape for {checkpoint}: {scores.shape}; "
            f"expected {(len(smiles), len(TASKS))}"
        )
    if not np.isfinite(scores).all():
        raise ValueError(f"Non-finite predictions from {checkpoint}")
    return smiles, scores


def _package_versions() -> dict[str, str]:
    packages = ["chemprop", "torch", "rdkit-pypi", "numpy", "pandas", "scikit-learn"]
    versions = {}
    for package in packages:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not-installed"
    return versions


def run_inference(
    input_csv: Path,
    checkpoint_dir: Path,
    output_csv: Path,
    manifest_path: Path,
    failure_csv: Path,
    id_column: str = "Molecule_Name",
    smiles_column: str = "SMILES",
    batch_size: int = 64,
    device: str = "cpu",
    limit: int | None = None,
) -> dict[str, Any]:
    input_csv = input_csv.resolve()
    checkpoint_dir = checkpoint_dir.resolve()
    frame = pd.read_csv(input_csv)
    missing = sorted({id_column, smiles_column} - set(frame.columns))
    if missing:
        raise ValueError(f"Input is missing required columns: {missing}")
    if frame[id_column].isna().any() or frame[id_column].duplicated().any():
        raise ValueError(f"{id_column} must be present and unique")
    if frame[smiles_column].isna().any():
        raise ValueError(f"{smiles_column} must be present for every input row")
    if limit is not None:
        frame = frame.head(limit).copy()
    frame = frame[[id_column, smiles_column]].rename(
        columns={id_column: "Molecule_Name", smiles_column: "model_smiles"}
    )
    if frame["model_smiles"].duplicated().any():
        raise ValueError("DeepMetab adapter requires unique input SMILES for stable ID recovery")

    checkpoints = [checkpoint_dir / f"model_{index}.pt" for index in range(1, 6)]
    missing_checkpoints = [str(path) for path in checkpoints if not path.is_file()]
    if missing_checkpoints:
        raise FileNotFoundError(f"Missing DeepMetab checkpoints: {missing_checkpoints}")

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    failure_csv.parent.mkdir(parents=True, exist_ok=True)
    model_input_path = output_csv.with_suffix(".chemprop_input.csv")
    _prepare_chemprop_input(frame, model_input_path)

    per_model: list[np.ndarray] = []
    scored_smiles: list[str] | None = None
    for index, checkpoint in enumerate(checkpoints, start=1):
        current_smiles, scores = _score_checkpoint(
            checkpoint, model_input_path, batch_size=batch_size, device=device
        )
        if scored_smiles is None:
            scored_smiles = current_smiles
        elif current_smiles != scored_smiles:
            raise ValueError(f"Prediction row order changed for checkpoint model_{index}.pt")
        per_model.append(scores)

    assert scored_smiles is not None
    score_stack = np.stack(per_model, axis=0)
    mean_scores = score_stack.mean(axis=0)
    variances = score_stack.var(axis=0, ddof=0)
    scored = pd.DataFrame({"model_smiles": scored_smiles})
    for task_index, task in enumerate(TASKS):
        for model_index in range(EXPECTED_CHECKPOINT_COUNT):
            scored[f"{task}_model_{model_index + 1}"] = score_stack[
                model_index, :, task_index
            ]
        scored[f"{task}_substrate_score"] = mean_scores[:, task_index]
        scored[f"{task}_ensemble_variance"] = variances[:, task_index]

    output = frame.merge(scored, on="model_smiles", how="left", validate="one_to_one")
    score_columns = [f"{task}_substrate_score" for task in TASKS]
    failed = output[score_columns].isna().any(axis=1)
    failures = output.loc[failed, ["Molecule_Name", "model_smiles"]].copy()
    failures["failure"] = "not_returned_by_chemprop"
    failures.to_csv(failure_csv, index=False)
    output.to_csv(output_csv, index=False)

    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "model": "DeepMetab substrate multitask ensemble",
        "upstream_repository": "https://github.com/YilingZhou/DeepMetab",
        "upstream_revision": "78c7511327f1a4042b61a64c44d35abb4c4b6b9c",
        "score_semantics": "unthresholded mean of five ChemProp classification scores; higher means more likely substrate",
        "variance_semantics": "population variance across five checkpoint scores (ddof=0)",
        "input_csv": str(input_csv),
        "input_sha256": sha256_file(input_csv),
        "input_rows": int(len(frame)),
        "scored_rows": int((~failed).sum()),
        "failure_rows": int(failed.sum()),
        "id_column": "Molecule_Name",
        "smiles_column": smiles_column,
        "representation": "original supplied SMILES; no TDI-label-dependent transformation",
        "tasks": TASKS,
        "batch_size": batch_size,
        "device": device,
        "checkpoint_key_adapter": "readout.* -> ffn.* followed by strict full-state load",
        "checkpoints": [
            {
                "path": str(path),
                "size": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in checkpoints
        ],
        "runtime": {
            "python": platform.python_version(),
            "packages": _package_versions(),
        },
        "outputs": {
            "predictions": {"path": str(output_csv), "sha256": sha256_file(output_csv)},
            "failures": {"path": str(failure_csv), "sha256": sha256_file(failure_csv)},
            "chemprop_input": {
                "path": str(model_input_path),
                "sha256": sha256_file(model_input_path),
            },
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--checkpoint-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--failures", type=Path, required=True)
    parser.add_argument("--id-column", default="Molecule_Name")
    parser.add_argument("--smiles-column", default="SMILES")
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--device", default="cpu", choices=["cpu", "cuda"])
    parser.add_argument("--limit", type=int)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    manifest = run_inference(
        input_csv=args.input,
        checkpoint_dir=args.checkpoint_dir,
        output_csv=args.output,
        manifest_path=args.manifest,
        failure_csv=args.failures,
        id_column=args.id_column,
        smiles_column=args.smiles_column,
        batch_size=args.batch_size,
        device=args.device,
        limit=args.limit,
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
