"""Metrics and scaffold-cluster bootstrap evaluation."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

SUPPORTED_ENDPOINTS = {"CYP2D6", "CYP3A4"}
TOP_FRACTIONS = (0.05, 0.10, 0.20)


@dataclass(frozen=True)
class Comparison:
    model: str
    score_name: str
    endpoint: str
    prediction_path: Path
    id_column: str
    score_column: str
    score_direction: str
    source_revision: str = ""
    notes: str = ""


def _validate_labels(labels: pd.DataFrame) -> None:
    required = {"Molecule_Name", "endpoint", "is_tdi", "scaffold_id"}
    missing = sorted(required - set(labels.columns))
    if missing:
        raise ValueError(f"Label table is missing columns: {missing}")
    if labels.duplicated(["Molecule_Name", "endpoint"]).any():
        raise ValueError("Label table has duplicate molecule/endpoint rows")
    values = set(pd.to_numeric(labels["is_tdi"], errors="coerce").dropna().unique())
    if not values <= {0, 1}:
        raise ValueError(f"Labels must be binary 0/1; observed {sorted(values)}")


def load_predictions(comparison: Comparison) -> pd.DataFrame:
    if comparison.endpoint not in SUPPORTED_ENDPOINTS:
        raise ValueError(f"Unsupported endpoint: {comparison.endpoint}")
    if comparison.score_direction not in {"higher", "lower"}:
        raise ValueError("score_direction must be 'higher' or 'lower'")
    if not comparison.prediction_path.is_file():
        raise FileNotFoundError(comparison.prediction_path)

    raw = pd.read_csv(comparison.prediction_path)
    missing = [c for c in (comparison.id_column, comparison.score_column) if c not in raw.columns]
    if missing:
        raise ValueError(f"Prediction file {comparison.prediction_path} is missing {missing}")
    pred = raw[[comparison.id_column, comparison.score_column]].rename(
        columns={comparison.id_column: "Molecule_Name", comparison.score_column: "raw_score"}
    )
    if pred["Molecule_Name"].isna().any() or pred["Molecule_Name"].duplicated().any():
        raise ValueError(f"Prediction IDs must be present and unique: {comparison.prediction_path}")

    original_nonblank = pred["raw_score"].notna()
    numeric = pd.to_numeric(pred["raw_score"], errors="coerce")
    bad_text = original_nonblank & numeric.isna()
    if bad_text.any():
        examples = pred.loc[bad_text, "raw_score"].head(5).tolist()
        raise ValueError(f"Nonnumeric prediction scores in {comparison.prediction_path}: {examples}")
    if np.isinf(numeric.dropna()).any():
        raise ValueError(f"Infinite prediction scores in {comparison.prediction_path}")
    pred["raw_score"] = numeric.astype(float)
    multiplier = 1.0 if comparison.score_direction == "higher" else -1.0
    pred["oriented_score"] = pred["raw_score"] * multiplier
    return pred


def join_comparison(labels: pd.DataFrame, comparison: Comparison) -> pd.DataFrame:
    _validate_labels(labels)
    endpoint_labels = labels.loc[labels["endpoint"] == comparison.endpoint].copy()
    if endpoint_labels.empty:
        raise ValueError(f"No labels found for endpoint {comparison.endpoint}")
    predictions = load_predictions(comparison)
    joined = endpoint_labels.merge(predictions, on="Molecule_Name", how="left", validate="one_to_one")
    joined.insert(0, "model", comparison.model)
    joined.insert(1, "score_name", comparison.score_name)
    joined["score_direction"] = comparison.score_direction
    joined["source_revision"] = comparison.source_revision
    return joined


def _enrichment(scored: pd.DataFrame, fraction: float) -> dict[str, float | int]:
    n = len(scored)
    if n == 0:
        return {
            "fraction": fraction, "k": 0, "hits": 0.0, "positive_rate": np.nan,
            "background_prevalence": np.nan, "enrichment": np.nan,
            "boundary_score": np.nan, "boundary_tie_size": 0,
            "boundary_fraction": np.nan,
        }
    ranked = scored.sort_values("oriented_score", ascending=False, kind="mergesort")
    k = max(1, math.ceil(fraction * n))
    boundary_score = float(ranked.iloc[k - 1]["oriented_score"])
    above = ranked[ranked["oriented_score"] > boundary_score]
    boundary = ranked[ranked["oriented_score"] == boundary_score]
    needed_from_boundary = k - len(above)
    boundary_fraction = needed_from_boundary / len(boundary)
    # When a top-k boundary cuts through a score tie, no member of that tie is
    # rank-identifiable. Use fractional expected inclusion rather than an
    # arbitrary molecule-ID ordering.
    hits = float(above["is_tdi"].sum()) + boundary_fraction * float(boundary["is_tdi"].sum())
    prevalence = float(scored["is_tdi"].mean())
    top_rate = hits / k
    lift = top_rate / prevalence if prevalence > 0 else np.nan
    return {
        "fraction": fraction,
        "k": k,
        "hits": hits,
        "positive_rate": top_rate,
        "background_prevalence": prevalence,
        "enrichment": lift,
        "boundary_score": boundary_score,
        "boundary_tie_size": int(len(boundary)),
        "boundary_fraction": float(boundary_fraction),
    }


def point_metrics(joined: pd.DataFrame) -> tuple[dict[str, float], list[dict[str, float | int]]]:
    finite = np.isfinite(joined["oriented_score"].to_numpy(dtype=float, na_value=np.nan))
    scored = joined.loc[finite].copy()
    coverage = float(len(scored) / len(joined)) if len(joined) else np.nan
    metrics = {"coverage": coverage, "pr_auc_average_precision": np.nan, "roc_auc": np.nan}
    if not scored.empty and scored["is_tdi"].nunique() == 2:
        metrics["pr_auc_average_precision"] = float(
            average_precision_score(scored["is_tdi"], scored["oriented_score"])
        )
        metrics["roc_auc"] = float(roc_auc_score(scored["is_tdi"], scored["oriented_score"]))
    enrichments = [_enrichment(scored, fraction) for fraction in TOP_FRACTIONS]
    for row in enrichments:
        metrics[f"enrichment_top_{int(row['fraction'] * 100)}pct"] = float(row["enrichment"])
    return metrics, enrichments


def scaffold_bootstrap(
    joined: pd.DataFrame,
    n_replicates: int,
    seed: int,
) -> pd.DataFrame:
    if n_replicates < 1:
        raise ValueError("n_replicates must be at least 1")
    group_codes, group_ids = pd.factorize(joined["scaffold_id"], sort=True)
    n_groups = len(group_ids)
    if n_groups == 0:
        raise ValueError("Cannot bootstrap an empty comparison")
    rng = np.random.default_rng(seed)
    rows = []
    for replicate in range(n_replicates):
        sampled = rng.integers(0, n_groups, size=n_groups)
        group_weights = np.bincount(sampled, minlength=n_groups)
        row_weights = group_weights[group_codes]
        boot_indices = np.repeat(np.arange(len(joined)), row_weights)
        boot = joined.iloc[boot_indices].reset_index(drop=True)
        metrics, _ = point_metrics(boot)
        rows.extend(
            {"replicate": replicate, "metric": metric, "estimate": value}
            for metric, value in metrics.items()
        )
    return pd.DataFrame(rows)


def _percentile_interval(values: Iterable[float], confidence: float) -> tuple[float, float, int]:
    array = np.asarray(list(values), dtype=float)
    array = array[np.isfinite(array)]
    if len(array) == 0:
        return np.nan, np.nan, 0
    alpha = (1.0 - confidence) / 2.0
    low, high = np.quantile(array, [alpha, 1.0 - alpha])
    return float(low), float(high), int(len(array))


def summarize_comparison(
    joined: pd.DataFrame,
    comparison: Comparison,
    n_replicates: int,
    seed: int,
    confidence: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if not 0 < confidence < 1:
        raise ValueError("confidence must lie strictly between 0 and 1")
    estimates, enrichments = point_metrics(joined)
    boot = scaffold_bootstrap(joined, n_replicates=n_replicates, seed=seed)
    finite = joined["oriented_score"].notna() & np.isfinite(joined["oriented_score"].fillna(0))
    common = {
        "model": comparison.model,
        "score_name": comparison.score_name,
        "endpoint": comparison.endpoint,
        "score_direction": comparison.score_direction,
        "n_labeled": int(len(joined)),
        "n_scored": int(finite.sum()),
        "n_positive_labeled": int(joined["is_tdi"].sum()),
        "n_positive_scored": int(joined.loc[finite, "is_tdi"].sum()),
        "n_scaffolds": int(joined["scaffold_id"].nunique()),
        "bootstrap_replicates": n_replicates,
        "confidence_level": confidence,
    }
    metric_rows = []
    for metric, estimate in estimates.items():
        low, high, valid = _percentile_interval(
            boot.loc[boot["metric"] == metric, "estimate"], confidence
        )
        metric_rows.append(
            {**common, "metric": metric, "estimate": estimate, "ci_low": low,
             "ci_high": high, "valid_bootstrap_replicates": valid}
        )
    metrics_frame = pd.DataFrame(metric_rows)

    enrichment_rows = []
    for row in enrichments:
        metric = f"enrichment_top_{int(row['fraction'] * 100)}pct"
        metric_summary = metrics_frame.loc[metrics_frame["metric"] == metric].iloc[0]
        enrichment_rows.append(
            {
                **common,
                **row,
                "ci_low": metric_summary["ci_low"],
                "ci_high": metric_summary["ci_high"],
                "valid_bootstrap_replicates": int(metric_summary["valid_bootstrap_replicates"]),
            }
        )
    enrichment_frame = pd.DataFrame(enrichment_rows)
    boot.insert(0, "model", comparison.model)
    boot.insert(1, "score_name", comparison.score_name)
    boot.insert(2, "endpoint", comparison.endpoint)
    return metrics_frame, enrichment_frame, boot
