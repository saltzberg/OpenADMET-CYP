#!/usr/bin/env python3
"""Verify all native-affine experiment artifacts and leakage/scope invariants."""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_experiment as R


def main() -> int:
    manifest = json.loads((R.EXP / "artifacts/run_manifest.json").read_text())
    assert manifest["status"] == "completed_experiment_not_submission"
    assert manifest["experiment_level"] == "adapted_native_model_affine_only"
    assert manifest["parameters"]["bootstrap_replicates"] == 2000
    assert manifest["parameters"]["clipping"] is False
    assert manifest["parameters"]["blind_label_use"] is False
    for relative, digest in manifest["source_files"].items():
        assert R.sha256_file(R.ROOT / relative) == digest, relative
    for relative, digest in manifest["outputs"].items():
        assert R.sha256_file(R.EXP / relative) == digest, relative

    source = pd.read_parquet(R.TRAIN_SOURCE)
    folds = pd.read_csv(R.FOLD_SOURCE)
    oof = pd.read_parquet(R.EXP / "artifacts/oof_predictions.parquet")
    assert len(oof) == 6525
    assert not oof.duplicated(["compound_id", "endpoint"]).any()
    assert np.isfinite(oof[["native_prediction", "calibrated_prediction"]]).all().all()
    joined = source.merge(oof, on=["compound_id", "endpoint"], suffixes=("_source", "_oof"), validate="one_to_one")
    assert np.array_equal(joined.prediction.to_numpy(), joined.native_prediction.to_numpy())
    for column in ["observed", "conf_low", "conf_high"]:
        assert np.array_equal(joined[f"{column}_source"].to_numpy(), joined[f"{column}_oof"].to_numpy())
    fold_check = oof[["compound_id", "fold", "scaffold_id"]].drop_duplicates().merge(
        folds[["compound_id", "fold", "scaffold_id"]], on="compound_id", suffixes=("_oof", "_source"), validate="one_to_one")
    assert np.array_equal(fold_check.fold_oof.to_numpy(), fold_check.fold_source.to_numpy())
    assert np.array_equal(fold_check.scaffold_id_oof.to_numpy(), fold_check.scaffold_id_source.to_numpy())
    assert oof.groupby(["endpoint", "fold"]).size().gt(0).all()

    coefficients = pd.read_csv(R.EXP / "artifacts/affine_coefficients.csv")
    outer = coefficients[coefficients.fit_scope == "outer_training"]
    final = coefficients[coefficients.fit_scope == "all_training_labels_for_blind_application"]
    assert len(outer) == 20 and outer.outer_fold.nunique() == 5
    assert len(final) == 4 and final.outer_fold.isna().all()
    assert (coefficients.slope > 0).all()
    for row in outer.itertuples(index=False):
        frame = oof[(oof.endpoint == row.endpoint) & (oof.fold == row.outer_fold)]
        expected = row.intercept + row.slope * frame.native_prediction.to_numpy()
        assert np.allclose(expected, frame.calibrated_prediction.to_numpy(), rtol=0, atol=1e-12)
        assert row.n_fit == int(((oof.endpoint == row.endpoint) & (oof.fold != row.outer_fold)).sum())
    for (endpoint, fold), frame in oof.groupby(["endpoint", "fold"]):
        raw_order = frame.sort_values(["native_prediction", "compound_id"]).compound_id.tolist()
        calibrated_order = frame.sort_values(["calibrated_prediction", "compound_id"]).compound_id.tolist()
        assert raw_order == calibrated_order, (endpoint, fold)

    blind_source = pd.read_csv(R.BLIND_SOURCE)
    blind = pd.read_csv(R.EXP / "artifacts/calibrated_blind_predictions.csv")
    assert len(blind) == 3000
    assert list(blind.columns) == ["compound_id", "SMILES", "endpoint", "native_prediction", "calibrated_prediction"]
    blind_join = blind_source.merge(blind, on=["compound_id", "SMILES", "endpoint"], validate="one_to_one")
    assert np.array_equal(blind_join.prediction.to_numpy(), blind_join.native_prediction.to_numpy())
    for row in final.itertuples(index=False):
        frame = blind[blind.endpoint == row.endpoint]
        assert np.allclose(row.intercept + row.slope * frame.native_prediction, frame.calibrated_prediction, rtol=0, atol=1e-12)

    metrics = pd.read_csv(R.EXP / "artifacts/metrics.csv")
    expected_metrics = R.compute_metrics(oof).sort_values(["model", "endpoint"]).reset_index(drop=True)
    actual_metrics = metrics.sort_values(["model", "endpoint"]).reset_index(drop=True)
    assert list(actual_metrics.columns) == list(expected_metrics.columns)
    assert np.allclose(actual_metrics.select_dtypes("number"), expected_metrics.select_dtypes("number"), rtol=1e-12, atol=1e-12)
    paired = pd.read_csv(R.EXP / "artifacts/paired_bootstrap_summary.csv")
    replicate = pd.read_parquet(R.EXP / "artifacts/paired_bootstrap_replicates.parquet")
    assert len(paired) == 5 and (paired.replicates == 2000).all()
    assert len(replicate) == 10000 and replicate.replicate.nunique() == 2000
    assert set(replicate.endpoint) == {*R.ENDPOINTS, "MA"}
    macro_delta = float(paired.loc[paired.endpoint == "MA", "delta_point_estimate"].iloc[0])
    endpoint_regression = bool((paired.loc[paired.endpoint.isin(R.ENDPOINTS), "ci_low"] > 0).any())
    assert manifest["decision_rule_met"] == (macro_delta < 0 and not endpoint_regression)

    readme = (R.EXP / "README.md").read_text()
    assert readme.count("<!-- EXECUTED_RESULTS_START -->") == 1
    assert "2,000-replicate paired scaffold-bootstrap" in readme
    assert f"decision rule was **{'MET' if manifest['decision_rule_met'] else 'NOT MET'}**" in readme
    html = (R.EXP / "index.html").read_text()
    assert html.count("<title>") == html.count("<main>") == html.count("<h1>") == 1
    assert "Created 2026-08-21 UTC" in html and "Last edited" in html
    assert "Only the affine stage is challenge-fitted" in html
    assert "Training predictions versus truth" in html and "Raw and linear-fit predictions" in html and html.count("<figure>") == 3
    image = Image.open(R.EXP / "figures/01_prediction_flow.png")
    image.verify()
    assert (R.EXP / "figures/01_prediction_flow.svg").stat().st_size > 1000
    scatter = pd.read_csv(R.EXP / "artifacts/linear_fit_scatter_source.csv")
    assert len(scatter) == 6525 and not scatter.duplicated(["compound_id", "endpoint"]).any()
    assert np.isfinite(scatter[["native_prediction", "calibrated_prediction"]]).all().all()
    scatter_image = Image.open(R.EXP / "figures/02_linear_fit_vs_observed.png")
    scatter_image.verify()
    assert (R.EXP / "figures/02_linear_fit_vs_observed.svg").stat().st_size > 1000
    overlay_image = Image.open(R.EXP / "figures/03_raw_and_linear_fit_vs_observed.png")
    overlay_image.verify()
    assert (R.EXP / "figures/03_raw_and_linear_fit_vs_observed.svg").stat().st_size > 1000
    print(f"PASS hashes, 6,525 OOF rows, 3,000 blind rows, 24 affine fits, 2,000 grouped paired bootstraps; decision_rule_met={manifest['decision_rule_met']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
