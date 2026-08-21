"""Command-line interface for the TDI benchmark."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from . import __version__
from .data import prepare_labels, sha256_file
from .evaluation import Comparison, join_comparison, summarize_comparison

REGISTRY_COLUMNS = {
    "model", "score_name", "endpoint", "prediction_path", "id_column",
    "score_column", "score_direction", "status",
}


def _comparison_from_row(row: pd.Series, registry_dir: Path) -> Comparison:
    prediction_path = Path(str(row["prediction_path"]))
    if not prediction_path.is_absolute():
        prediction_path = (registry_dir / prediction_path).resolve()
    return Comparison(
        model=str(row["model"]),
        score_name=str(row["score_name"]),
        endpoint=str(row["endpoint"]),
        prediction_path=prediction_path,
        id_column=str(row["id_column"]),
        score_column=str(row["score_column"]),
        score_direction=str(row["score_direction"]).lower(),
        source_revision="" if pd.isna(row.get("source_revision")) else str(row.get("source_revision")),
        notes="" if pd.isna(row.get("notes")) else str(row.get("notes")),
    )


def evaluate_registry(
    labels_path: Path,
    registry_path: Path,
    output_dir: Path,
    bootstrap_replicates: int,
    seed: int,
    confidence: float,
    overwrite: bool,
) -> dict[str, object]:
    labels_path = labels_path.resolve()
    registry_path = registry_path.resolve()
    if output_dir.exists() and any(output_dir.iterdir()) and not overwrite:
        raise FileExistsError(f"Output directory is not empty: {output_dir}; use --overwrite")
    output_dir.mkdir(parents=True, exist_ok=True)

    labels = pd.read_csv(labels_path)
    registry = pd.read_csv(registry_path, keep_default_na=True)
    missing = sorted(REGISTRY_COLUMNS - set(registry.columns))
    if missing:
        raise ValueError(f"Model registry is missing columns: {missing}")
    duplicate_keys = registry.duplicated(["model", "score_name", "endpoint"])
    if duplicate_keys.any():
        raise ValueError("Model registry has duplicate model/score_name/endpoint rows")

    ready = registry.loc[registry["status"].astype(str).str.lower() == "ready"]
    if ready.empty:
        raise ValueError("Model registry has no rows with status=ready")
    if "source_revision" not in ready.columns or ready["source_revision"].isna().any():
        raise ValueError("Every status=ready row must record a source_revision")
    if ready["source_revision"].astype(str).str.strip().eq("").any():
        raise ValueError("Every status=ready row must record a source_revision")

    metric_frames = []
    enrichment_frames = []
    bootstrap_frames = []
    joined_frames = []
    input_records = []
    for index, row in ready.iterrows():
        comparison = _comparison_from_row(row, registry_path.parent)
        joined = join_comparison(labels, comparison)
        metrics, enrichment, bootstrap = summarize_comparison(
            joined,
            comparison,
            n_replicates=bootstrap_replicates,
            seed=seed + int(index),
            confidence=confidence,
        )
        metric_frames.append(metrics)
        enrichment_frames.append(enrichment)
        bootstrap_frames.append(bootstrap)
        joined_frames.append(joined)
        input_records.append(
            {
                "model": comparison.model,
                "score_name": comparison.score_name,
                "endpoint": comparison.endpoint,
                "prediction_path": str(comparison.prediction_path),
                "prediction_sha256": sha256_file(comparison.prediction_path),
                "score_column": comparison.score_column,
                "score_direction": comparison.score_direction,
                "source_revision": comparison.source_revision,
                "notes": comparison.notes,
            }
        )

    outputs = {
        "metrics.csv": pd.concat(metric_frames, ignore_index=True),
        "enrichment.csv": pd.concat(enrichment_frames, ignore_index=True),
        "bootstrap_metrics.csv": pd.concat(bootstrap_frames, ignore_index=True),
        "joined_predictions.csv": pd.concat(joined_frames, ignore_index=True),
    }
    for name, frame in outputs.items():
        frame.to_csv(output_dir / name, index=False)

    registry_status = registry.copy()
    registry_status["evaluated"] = registry_status["status"].astype(str).str.lower() == "ready"
    registry_status.to_csv(output_dir / "registry_status.csv", index=False)
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "package_version": __version__,
        "labels_path": str(labels_path),
        "labels_sha256": sha256_file(labels_path),
        "registry_path": str(registry_path),
        "registry_sha256": sha256_file(registry_path),
        "bootstrap": {
            "unit": "Bemis-Murcko scaffold cluster over the full labelled endpoint universe",
            "replicates": bootstrap_replicates,
            "seed": seed,
            "confidence": confidence,
            "interval": "percentile",
        },
        "metric_definitions": {
            "primary": "pr_auc_average_precision (sklearn average_precision_score)",
            "roc_auc": "sklearn roc_auc_score",
            "coverage": "finite scored labelled molecules / all labelled molecules for endpoint",
            "enrichment": "positive rate in ceil(fraction * n_scored) / scored-set prevalence",
            "ranking_ties": "fractional expected inclusion when the top-k boundary cuts a score tie",
        },
        "comparisons": input_records,
        "outputs": {name: sha256_file(output_dir / name) for name in outputs},
    }
    (output_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="tdi-benchmark")
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare-labels", help="Build long-form TDI labels")
    prepare.add_argument("--source", type=Path, required=True)
    prepare.add_argument("--output-dir", type=Path, required=True)

    evaluate = subparsers.add_parser("evaluate-registry", help="Evaluate ready registry rows")
    evaluate.add_argument("--labels", type=Path, required=True)
    evaluate.add_argument("--registry", type=Path, required=True)
    evaluate.add_argument("--output-dir", type=Path, required=True)
    evaluate.add_argument("--bootstrap-replicates", type=int, default=2000)
    evaluate.add_argument("--seed", type=int, default=20260820)
    evaluate.add_argument("--confidence", type=float, default=0.95)
    evaluate.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    if args.command == "prepare-labels":
        result = prepare_labels(args.source, args.output_dir)
        print(json.dumps(result, indent=2))
    elif args.command == "evaluate-registry":
        result = evaluate_registry(
            args.labels,
            args.registry,
            args.output_dir,
            args.bootstrap_replicates,
            args.seed,
            args.confidence,
            args.overwrite,
        )
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
