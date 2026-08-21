import numpy as np
import pandas as pd
import pytest

from tdi_public_model_benchmark.evaluation import (
    Comparison,
    join_comparison,
    point_metrics,
    scaffold_bootstrap,
    summarize_comparison,
)


def label_frame(n=20):
    return pd.DataFrame(
        {
            "Molecule_Name": [f"M{i:02d}" for i in range(n)],
            "endpoint": ["CYP2D6"] * n,
            "is_tdi": [0] * (n // 2) + [1] * (n - n // 2),
            "scaffold_id": [f"S{i // 2:02d}" for i in range(n)],
            "SMILES": ["CC"] * n,
        }
    )


def write_predictions(tmp_path, scores, name="score"):
    path = tmp_path / "predictions.csv"
    pd.DataFrame(
        {"Molecule_Name": [f"M{i:02d}" for i in range(len(scores))], name: scores}
    ).to_csv(path, index=False)
    return path


def comparison(path, direction="higher"):
    return Comparison(
        model="Synthetic",
        score_name="synthetic_score",
        endpoint="CYP2D6",
        prediction_path=path,
        id_column="Molecule_Name",
        score_column="score",
        score_direction=direction,
    )


def test_perfect_ranking_and_enrichment(tmp_path):
    labels = label_frame()
    path = write_predictions(tmp_path, np.arange(20), "score")
    joined = join_comparison(labels, comparison(path))
    metrics, enrichment = point_metrics(joined)
    assert metrics["pr_auc_average_precision"] == pytest.approx(1.0)
    assert metrics["roc_auc"] == pytest.approx(1.0)
    assert metrics["coverage"] == pytest.approx(1.0)
    assert [row["k"] for row in enrichment] == [1, 2, 4]
    assert all(row["enrichment"] == pytest.approx(2.0) for row in enrichment)


def test_lower_direction_is_oriented(tmp_path):
    labels = label_frame()
    path = write_predictions(tmp_path, -np.arange(20), "score")
    joined = join_comparison(labels, comparison(path, direction="lower"))
    metrics, _ = point_metrics(joined)
    assert metrics["pr_auc_average_precision"] == pytest.approx(1.0)
    assert metrics["roc_auc"] == pytest.approx(1.0)


def test_partial_coverage_is_not_imputed(tmp_path):
    labels = label_frame()
    scores = np.arange(20, dtype=float)
    scores[::2] = np.nan
    path = write_predictions(tmp_path, scores, "score")
    joined = join_comparison(labels, comparison(path))
    metrics, _ = point_metrics(joined)
    assert metrics["coverage"] == pytest.approx(0.5)
    assert joined["raw_score"].isna().sum() == 10


def test_enrichment_uses_fractional_inclusion_at_score_ties(tmp_path):
    labels = label_frame()
    path = write_predictions(tmp_path, np.ones(20), "score")
    joined = join_comparison(labels, comparison(path))
    metrics, enrichment = point_metrics(joined)
    assert all(row["enrichment"] == pytest.approx(1.0) for row in enrichment)
    assert enrichment[0]["k"] == 1
    assert enrichment[0]["boundary_tie_size"] == 20
    assert enrichment[0]["boundary_fraction"] == pytest.approx(0.05)
    assert metrics["enrichment_top_5pct"] == pytest.approx(1.0)


def test_scaffold_bootstrap_and_summary(tmp_path):
    labels = label_frame(40)
    path = write_predictions(tmp_path, np.arange(40), "score")
    comp = comparison(path)
    joined = join_comparison(labels, comp)
    boot = scaffold_bootstrap(joined, n_replicates=30, seed=7)
    assert boot["replicate"].nunique() == 30
    assert set(boot["metric"]) >= {
        "pr_auc_average_precision", "roc_auc", "coverage",
        "enrichment_top_5pct", "enrichment_top_10pct", "enrichment_top_20pct",
    }
    metrics, enrichment, _ = summarize_comparison(
        joined, comp, n_replicates=30, seed=7, confidence=0.95
    )
    assert (metrics["ci_low"] <= metrics["ci_high"]).all()
    assert len(enrichment) == 3


def test_duplicate_prediction_ids_fail(tmp_path):
    labels = label_frame(4)
    path = tmp_path / "predictions.csv"
    pd.DataFrame({"Molecule_Name": ["M00", "M00"], "score": [0.1, 0.2]}).to_csv(
        path, index=False
    )
    with pytest.raises(ValueError, match="unique"):
        join_comparison(labels, comparison(path))
