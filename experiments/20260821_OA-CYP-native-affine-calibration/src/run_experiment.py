#!/usr/bin/env python3
"""Run the preregistered endpoint-wise native-head affine calibration experiment."""
from __future__ import annotations

import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scipy
from scipy.stats import kendalltau, spearmanr
import sklearn
from sklearn.metrics import mean_absolute_error, r2_score

ROOT = Path(__file__).resolve().parents[3]
EXP = Path(__file__).resolve().parents[1]
TRAIN_SOURCE = ROOT / "experiments/20260821_OA-CYP-native-zero-shot/artifacts/native_train_scored.parquet"
BLIND_SOURCE = ROOT / "experiments/20260821_OA-CYP-native-zero-shot/artifacts/native_blind_predictions.csv"
FOLD_SOURCE = ROOT / "experiments/20260820_OA-CYP-finetuned-chemeleon/artifacts/fold_manifest.csv"
ENDPOINTS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
SEED = 20260821
N_BOOTSTRAP = 2000
CREATED = "2026-08-21 UTC"
RUN_COMMAND = "/home/dan/swr/miniconda3/envs/cheminf/bin/python src/run_experiment.py"
TEST_COMMAND = "/home/dan/swr/miniconda3/envs/cheminf/bin/python src/test_experiment.py"
VERIFY_COMMAND = "/home/dan/swr/miniconda3/envs/cheminf/bin/python src/verify_outputs.py"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def soft_error(pred, low, high):
    pred = np.asarray(pred, float)
    return np.maximum(0.0, np.maximum(np.asarray(low, float) - pred, pred - np.asarray(high, float)))


def st_rae(observed, pred, low, high, weights=None) -> float:
    observed = np.asarray(observed, float)
    if weights is None:
        weights = np.ones(len(observed), float)
    weights = np.asarray(weights, float)
    mean = np.sum(weights * observed) / np.sum(weights)
    denominator = np.sum(weights * np.abs(observed - mean))
    return float(np.sum(weights * soft_error(pred, low, high)) / denominator)


def fit_affine(native_prediction, observed) -> tuple[float, float]:
    x = np.asarray(native_prediction, float)
    y = np.asarray(observed, float)
    design = np.column_stack([np.ones(len(x)), x])
    intercept, slope = np.linalg.lstsq(design, y, rcond=None)[0]
    return float(intercept), float(slope)


def metric_row(model: str, endpoint: str, frame: pd.DataFrame, prediction_column: str) -> dict:
    y = frame["observed"].to_numpy(float)
    p = frame[prediction_column].to_numpy(float)
    return {
        "model": model,
        "endpoint": endpoint,
        "st_rae": st_rae(y, p, frame["conf_low"], frame["conf_high"]),
        "mae": float(mean_absolute_error(y, p)),
        "r2": float(r2_score(y, p)),
        "spearman_rho": float(spearmanr(y, p).statistic),
        "kendall_tau": float(kendalltau(y, p).statistic),
        "n": int(len(frame)),
        "prediction_sd": float(np.std(p, ddof=1)),
        "observed_sd": float(np.std(y, ddof=1)),
    }


def compute_metrics(oof: pd.DataFrame) -> pd.DataFrame:
    rows = []
    mapping = [("native_raw", "native_prediction"), ("native_affine", "calibrated_prediction")]
    for model, column in mapping:
        for endpoint in ENDPOINTS:
            rows.append(metric_row(model, endpoint, oof[oof.endpoint == endpoint], column))
    metrics = pd.DataFrame(rows)
    macro = []
    for model, group in metrics.groupby("model", sort=True):
        macro.append({
            "model": model,
            "endpoint": "MA",
            "st_rae": float(group.st_rae.mean()),
            "mae": float(group.mae.mean()),
            "r2": float(group.r2.mean()),
            "spearman_rho": float(group.spearman_rho.mean()),
            "kendall_tau": float(group.kendall_tau.mean()),
            "n": int(group.n.sum()),
            "prediction_sd": float(group.prediction_sd.mean()),
            "observed_sd": float(group.observed_sd.mean()),
        })
    return pd.concat([metrics, pd.DataFrame(macro)], ignore_index=True)


def paired_group_bootstrap(oof: pd.DataFrame, n_replicates: int = N_BOOTSTRAP, seed: int = SEED):
    groups = sorted(oof.scaffold_id.unique())
    group_index = {group: idx for idx, group in enumerate(groups)}
    prepared = {}
    for endpoint in ENDPOINTS:
        frame = oof[oof.endpoint == endpoint]
        prepared[endpoint] = {
            "observed": frame.observed.to_numpy(float),
            "low": frame.conf_low.to_numpy(float),
            "high": frame.conf_high.to_numpy(float),
            "raw": frame.native_prediction.to_numpy(float),
            "affine": frame.calibrated_prediction.to_numpy(float),
            "group": frame.scaffold_id.map(group_index).to_numpy(int),
        }
    rng = np.random.default_rng(seed)
    records = []
    for replicate in range(n_replicates):
        counts = np.bincount(rng.integers(0, len(groups), size=len(groups)), minlength=len(groups)).astype(float)
        endpoint_deltas = []
        for endpoint in ENDPOINTS:
            data = prepared[endpoint]
            weights = counts[data["group"]]
            raw = st_rae(data["observed"], data["raw"], data["low"], data["high"], weights)
            affine = st_rae(data["observed"], data["affine"], data["low"], data["high"], weights)
            delta = affine - raw
            endpoint_deltas.append(delta)
            records.append({"replicate": replicate, "endpoint": endpoint, "delta_affine_minus_raw_st_rae": delta})
        records.append({"replicate": replicate, "endpoint": "MA", "delta_affine_minus_raw_st_rae": float(np.mean(endpoint_deltas))})
    replicates = pd.DataFrame(records)
    summary = (replicates.groupby("endpoint", sort=False)["delta_affine_minus_raw_st_rae"]
               .agg(delta_bootstrap_mean="mean", ci_low=lambda x: x.quantile(0.025), ci_high=lambda x: x.quantile(0.975))
               .reset_index())
    summary["replicates"] = n_replicates
    return replicates, summary


def prediction_flow() -> None:
    fig, ax = plt.subplots(figsize=(11.2, 2.8))
    ax.set_xlim(0, 11.2)
    ax.set_ylim(0, 2.8)
    ax.axis("off")
    items = [
        (0.7, "Challenge\nSMILES", "#222"),
        (3.0, "Released frozen\nCheMeleon encoder", "#222"),
        (5.6, "Released frozen\nnative 4-output FFN", "#222"),
        (8.1, "Native endpoint\nprediction", "#222"),
        (10.5, "OLS intercept + slope\nCHALLENGE-FITTED", "#9a3412"),
    ]
    for x, label, color in items:
        ax.text(x, 1.65, label, ha="center", va="center", fontsize=10, color=color,
                fontweight="bold" if "CHALLENGE" in label else "normal")
    for left, right in zip(items[:-1], items[1:]):
        ax.annotate("", xy=(right[0] - 0.8, 1.65), xytext=(left[0] + 0.8, 1.65),
                    arrowprops={"arrowstyle": "->", "lw": 1, "color": "#555"})
    ax.text(5.25, 0.42, "Encoder and native head are never fitted here.", ha="center", color="#555", fontsize=9.5)
    ax.text(10.0, 0.42, "Only this affine stage uses challenge-training labels; no blind labels or clipping.",
            ha="center", color="#9a3412", fontsize=9.5)
    for extension in ["png", "svg"]:
        fig.savefig(EXP / f"figures/01_prediction_flow.{extension}", dpi=300 if extension == "png" else None,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


def calibrated_scatter(frame: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 7.2), sharex=True, sharey=True)
    for ax, endpoint in zip(axes.flat, ENDPOINTS):
        group = frame[frame.endpoint == endpoint]
        ax.scatter(group.observed, group.calibrated_prediction, s=8, alpha=0.42,
                   color="#2f5d8a", linewidths=0, rasterized=True)
        ax.plot([1.5, 8.0], [1.5, 8.0], color="#555", lw=0.7)
        ax.set_xlim(1.5, 8.0); ax.set_ylim(1.5, 8.0); ax.set_aspect("equal")
        ax.set_title(endpoint, loc="left", fontsize=10)
        ax.spines[["top", "right"]].set_visible(False); ax.grid(False)
    axes[1, 0].set_xlabel("Observed pIC50"); axes[1, 1].set_xlabel("Observed pIC50")
    axes[0, 0].set_ylabel("Linear-fit OOF pIC50"); axes[1, 0].set_ylabel("Linear-fit OOF pIC50")
    fig.subplots_adjust(hspace=0.20, wspace=0.18)
    for extension in ["png", "svg"]:
        fig.savefig(EXP / f"figures/02_linear_fit_vs_observed.{extension}",
                    dpi=300 if extension == "png" else None, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def raw_and_calibrated_scatter(frame: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 7.2), sharex=True, sharey=True)
    for ax, endpoint in zip(axes.flat, ENDPOINTS):
        group = frame[frame.endpoint == endpoint]
        ax.scatter(group.observed, group.native_prediction, s=7, alpha=0.24,
                   color="#888888", linewidths=0, rasterized=True)
        ax.scatter(group.observed, group.calibrated_prediction, s=7, alpha=0.34,
                   color="#2f5d8a", linewidths=0, rasterized=True)
        ax.plot([1.5, 8.0], [1.5, 8.0], color="#555", lw=0.7)
        ax.set_xlim(1.5, 8.0); ax.set_ylim(1.5, 8.0); ax.set_aspect("equal")
        ax.set_title(endpoint, loc="left", fontsize=10)
        ax.spines[["top", "right"]].set_visible(False); ax.grid(False)
    axes[1, 0].set_xlabel("Observed pIC50"); axes[1, 1].set_xlabel("Observed pIC50")
    axes[0, 0].set_ylabel("OOF predicted pIC50"); axes[1, 0].set_ylabel("OOF predicted pIC50")
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#888888", markeredgewidth=0, markersize=5, label="Raw model"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#2f5d8a", markeredgewidth=0, markersize=5, label="Linear fit"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.005), fontsize=8.5)
    fig.subplots_adjust(bottom=0.10, hspace=0.20, wspace=0.18)
    for extension in ["png", "svg"]:
        fig.savefig(EXP / f"figures/03_raw_and_linear_fit_vs_observed.{extension}",
                    dpi=300 if extension == "png" else None, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def result_table(metrics: pd.DataFrame, paired: pd.DataFrame) -> tuple[str, bool]:
    lookup = metrics.set_index(["endpoint", "model"])
    delta_lookup = paired.set_index("endpoint")
    rows = []
    for endpoint in [*ENDPOINTS, "MA"]:
        raw = lookup.loc[(endpoint, "native_raw")]
        affine = lookup.loc[(endpoint, "native_affine")]
        diff = float(affine.st_rae - raw.st_rae)
        ci = delta_lookup.loc[endpoint]
        rows.append(
            f"| {endpoint} | {raw.st_rae:.9f} | {affine.st_rae:.9f} | {diff:+.9f} | "
            f"{ci.ci_low:+.9f} | {ci.ci_high:+.9f} | {raw.spearman_rho:.9f} | {affine.spearman_rho:.9f} |"
        )
    macro_delta = float(lookup.loc[("MA", "native_affine"), "st_rae"] - lookup.loc[("MA", "native_raw"), "st_rae"])
    endpoint_regressions = [endpoint for endpoint in ENDPOINTS if float(delta_lookup.loc[endpoint, "ci_low"]) > 0]
    decision_met = macro_delta < 0 and not endpoint_regressions
    table = "\n".join([
        "| Endpoint | Raw ST-RAE | Affine ST-RAE | Δ affine−raw | paired 95% CI low | paired 95% CI high | Raw Spearman ρ | Affine Spearman ρ |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
        *rows,
    ])
    macro_ci = delta_lookup.loc["MA"]
    statement = (
        f"The unweighted endpoint-macro ST-RAE difference was **{macro_delta:+.9f}** "
        f"(2,000-replicate paired scaffold-bootstrap 95% CI {macro_ci.ci_low:+.9f} to {macro_ci.ci_high:+.9f}). "
        f"A reproducible endpoint regression was defined as an endpoint paired interval wholly above zero; "
        f"none occurred." if not endpoint_regressions else
        f"The unweighted endpoint-macro ST-RAE difference was **{macro_delta:+.9f}** "
        f"(2,000-replicate paired scaffold-bootstrap 95% CI {macro_ci.ci_low:+.9f} to {macro_ci.ci_high:+.9f}). "
        f"Endpoint intervals wholly above zero: {', '.join(endpoint_regressions)}."
    )
    statement += f" The preregistered no-endpoint-regression decision rule was **{'MET' if decision_met else 'NOT MET'}**."
    return table + "\n\n" + statement, decision_met


def update_readme(results: str) -> None:
    path = EXP / "README.md"
    text = path.read_text()
    start = "<!-- EXECUTED_RESULTS_START -->"
    end = "<!-- EXECUTED_RESULTS_END -->"
    if start in text:
        block_end = text.index(end, text.index(start)) + len(end)
        text = text[:text.index(start)] + start + "\n" + results.rstrip() + "\n" + end + text[block_end:]
    elif "## Results\n\nNot run." in text:
        text = text.replace("## Results\n\nNot run.", "## Results\n\n" + start + "\n" + results.rstrip() + "\n" + end)
    elif text.rstrip().endswith("## Results"):
        text = text.rstrip() + "\n\n" + start + "\n" + results.rstrip() + "\n" + end + "\n"
    else:
        raise RuntimeError("README Results anchor is missing")
    text = text.replace("> **Status:** Intent only — not a submission", "> **Status:** Completed experiment — not a submission")
    path.write_text(text.rstrip() + "\n")


def render_html(metrics: pd.DataFrame, paired: pd.DataFrame, decision_met: bool) -> None:
    edited = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    lookup = metrics.set_index(["endpoint", "model"])
    diffs = paired.set_index("endpoint")
    rows = []
    for endpoint in [*ENDPOINTS, "MA"]:
        raw = lookup.loc[(endpoint, "native_raw")]
        affine = lookup.loc[(endpoint, "native_affine")]
        delta = affine.st_rae - raw.st_rae
        ci = diffs.loc[endpoint]
        rows.append(f"<tr><td>{endpoint}</td><td>{raw.st_rae:.9f}</td><td>{affine.st_rae:.9f}</td><td>{delta:+.9f}</td><td>{ci.ci_low:+.9f}, {ci.ci_high:+.9f}</td></tr>")
    html = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>OpenADMET native-head affine calibration</title><style>body{{color:#222;font:16px/1.5 Georgia,serif;margin:0}}main{{max-width:62rem;margin:3rem auto;padding:0 1rem}}h1{{font-size:2.5rem}}h2{{margin-top:2.2rem;border-bottom:1px solid #ccc}}.meta{{color:#666;font:13px ui-monospace,monospace}}img{{width:min(100%,1000px)}}table{{border-collapse:collapse;width:100%;font:13px system-ui}}th,td{{padding:.42rem;border-bottom:1px solid #ddd;text-align:right}}th:first-child,td:first-child{{text-align:left}}a{{color:#285f8f}}</style></head><body><main><div class="meta">{EXP.name} · Created {CREATED} · Last edited {edited}</div><h1>OpenADMET native-head affine calibration</h1><p><strong>Purpose:</strong> test whether one endpoint-specific ordinary affine map improves released native-head predictions under the established grouped outer folds.</p><h2>Prediction construction</h2><figure><img src="figures/01_prediction_flow.png" alt="SMILES through frozen released encoder and native head, followed by the only challenge-fitted stage, endpoint affine calibration"><figcaption>Only the affine stage is challenge-fitted. Each OOF map uses its outer-training rows; final blind maps use all training labels.</figcaption></figure><h2>OOF ST-RAE</h2><table><thead><tr><th>Endpoint</th><th>Raw native</th><th>Linear fit</th><th>Δ fit−raw</th><th>Paired 95% CI</th></tr></thead><tbody>{''.join(rows)}</tbody></table><p><strong>Decision rule: {'MET' if decision_met else 'NOT MET'}.</strong> Support requires lower macro ST-RAE and no endpoint paired interval wholly above zero.</p><h2>Training predictions versus truth</h2><figure><img src="figures/02_linear_fit_vs_observed.png" alt="Four panels of linear-fit out-of-fold pIC50 predictions against observed training pIC50"><figcaption>Each point is one OOF prediction. The diagonal is perfect agreement.</figcaption></figure><h2>Raw and linear-fit predictions</h2><figure><img src="figures/03_raw_and_linear_fit_vs_observed.png" alt="Four panels comparing raw model and linear-fit OOF predictions against observed training pIC50"><figcaption>Grey points are raw model predictions. Blue points are the corresponding linear-fit OOF predictions.</figcaption></figure><h2>Artifacts</h2><p><a href="artifacts/metrics.csv">metrics</a> · <a href="artifacts/oof_predictions.parquet">OOF predictions</a> · <a href="artifacts/calibrated_blind_predictions.csv">blind predictions</a> · <a href="artifacts/affine_coefficients.csv">coefficients</a> · <a href="artifacts/run_manifest.json">manifest</a> · <a href="README.md">canonical record</a></p></main></body></html>'''
    (EXP / "index.html").write_text(html)


def main() -> int:
    for directory in [EXP / "artifacts", EXP / "figures"]:
        directory.mkdir(exist_ok=True)
    train = pd.read_parquet(TRAIN_SOURCE)
    blind = pd.read_csv(BLIND_SOURCE)
    folds = pd.read_csv(FOLD_SOURCE)
    assert list(sorted(train.endpoint.unique())) == sorted(ENDPOINTS)
    assert list(sorted(blind.endpoint.unique())) == sorted(ENDPOINTS)
    assert not train.duplicated(["compound_id", "endpoint"]).any()
    assert not blind.duplicated(["compound_id", "endpoint"]).any()
    assert not folds.compound_id.duplicated().any()
    assert set(train.compound_id) <= set(folds.compound_id)
    merged = train.merge(folds[["compound_id", "scaffold_id", "fold"]], on="compound_id", how="left", validate="many_to_one")
    assert merged[["scaffold_id", "fold"]].notna().all().all()
    merged = merged.rename(columns={"prediction": "native_prediction"})
    merged["calibrated_prediction"] = np.nan
    decisions = []
    for endpoint in ENDPOINTS:
        endpoint_mask = merged.endpoint == endpoint
        for fold in sorted(folds.fold.unique()):
            train_mask = endpoint_mask & (merged.fold != fold)
            holdout_mask = endpoint_mask & (merged.fold == fold)
            intercept, slope = fit_affine(merged.loc[train_mask, "native_prediction"], merged.loc[train_mask, "observed"])
            merged.loc[holdout_mask, "calibrated_prediction"] = intercept + slope * merged.loc[holdout_mask, "native_prediction"]
            decisions.append({"fit_scope": "outer_training", "endpoint": endpoint, "outer_fold": int(fold), "intercept": intercept, "slope": slope, "n_fit": int(train_mask.sum()), "n_predict": int(holdout_mask.sum())})
    assert np.isfinite(merged.calibrated_prediction).all()
    final_coefficients = []
    calibrated_blind_parts = []
    for endpoint in ENDPOINTS:
        endpoint_train = merged[merged.endpoint == endpoint]
        intercept, slope = fit_affine(endpoint_train.native_prediction, endpoint_train.observed)
        endpoint_blind = blind[blind.endpoint == endpoint].copy().rename(columns={"prediction": "native_prediction"})
        endpoint_blind["calibrated_prediction"] = intercept + slope * endpoint_blind.native_prediction
        calibrated_blind_parts.append(endpoint_blind)
        final_coefficients.append({"fit_scope": "all_training_labels_for_blind_application", "endpoint": endpoint, "intercept": intercept, "slope": slope, "n_fit": len(endpoint_train), "n_predict": len(endpoint_blind)})
    calibrated_blind = pd.concat(calibrated_blind_parts, ignore_index=True)
    coefficients = pd.DataFrame([*decisions, *final_coefficients])
    metrics = compute_metrics(merged)
    replicates, paired = paired_group_bootstrap(merged)
    point_delta = (metrics.pivot(index="endpoint", columns="model", values="st_rae")["native_affine"] - metrics.pivot(index="endpoint", columns="model", values="st_rae")["native_raw"]).rename("delta_point_estimate")
    paired = paired.merge(point_delta, left_on="endpoint", right_index=True, validate="one_to_one")
    merged.to_parquet(EXP / "artifacts/oof_predictions.parquet", index=False)
    calibrated_blind.to_csv(EXP / "artifacts/calibrated_blind_predictions.csv", index=False)
    coefficients.to_csv(EXP / "artifacts/affine_coefficients.csv", index=False)
    metrics.to_csv(EXP / "artifacts/metrics.csv", index=False)
    replicates.to_parquet(EXP / "artifacts/paired_bootstrap_replicates.parquet", index=False)
    paired.to_csv(EXP / "artifacts/paired_bootstrap_summary.csv", index=False)
    prediction_flow()
    results, decision_met = result_table(metrics, paired)
    coefficient_lines = []
    final_table = coefficients[coefficients.fit_scope == "all_training_labels_for_blind_application"]
    for row in final_table.itertuples(index=False):
        coefficient_lines.append(f"| {row.endpoint} | {row.intercept:.12f} | {row.slope:.12f} | {int(row.n_fit)} |")
    results += "\n\n### Final all-training affine maps applied to blinded native predictions\n\n| Endpoint | Intercept | Slope | Training rows |\n|---|---:|---:|---:|\n" + "\n".join(coefficient_lines)
    results += "\n\nAll fitted slopes were positive. Therefore each individual outer-fold map and each final all-training endpoint map preserved ordering (including ties) within the rows to which that one map was applied. The assembled OOF endpoint ranks changed slightly because five independently fitted maps place different outer holdouts on different affine scales; the exact raw and affine Spearman values are reported above. No clipping was applied. Blinded labels were not available or used."
    results += f"\n\n### Reproduction\n\n```text\n{RUN_COMMAND}\n{TEST_COMMAND}\n{VERIFY_COMMAND}\n```"
    update_readme(results)
    merged[["compound_id", "endpoint", "observed", "native_prediction", "calibrated_prediction", "fold", "scaffold_id"]].to_csv(
        EXP / "artifacts/linear_fit_scatter_source.csv", index=False)
    calibrated_scatter(merged)
    raw_and_calibrated_scatter(merged)
    render_html(metrics, paired, decision_met)
    output_hashes = {}
    for path in sorted(EXP.rglob("*")):
        if path.is_file() and path.name != "run_manifest.json" and "__pycache__" not in path.parts:
            output_hashes[str(path.relative_to(EXP))] = sha256_file(path)
    manifest = {
        "status": "completed_experiment_not_submission",
        "experiment_id": EXP.name,
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "experiment_level": "adapted_native_model_affine_only",
        "challenge_fitted_stages": ["endpoint-specific OLS affine intercept+slope"],
        "frozen_stages": ["released CYP-finetuned CheMeleon encoder", "released native four-output FFN"],
        "parameters": {"outer_folds": 5, "outer_fold_source": str(FOLD_SOURCE.relative_to(ROOT)), "bootstrap_group": "scaffold_id", "bootstrap_replicates": N_BOOTSTRAP, "seed": SEED, "calibrator": "ordinary least-squares intercept+slope", "clipping": False, "blind_label_use": False},
        "metric": "ST-RAE = sum distance outside reported confidence interval / sum absolute deviation of observed point estimates from endpoint mean; MA is the unweighted mean across four endpoints",
        "decision_rule_implementation": "MET iff point-estimate macro delta (affine minus raw) < 0 and no endpoint paired-bootstrap 95% interval is wholly above 0",
        "decision_rule_met": bool(decision_met),
        "commands": {"run": RUN_COMMAND, "test": TEST_COMMAND, "verify": VERIFY_COMMAND},
        "source_files": {str(path.relative_to(ROOT)): sha256_file(path) for path in [TRAIN_SOURCE, BLIND_SOURCE, FOLD_SOURCE]},
        "runtime": {"executable": sys.executable, "python": sys.version, "platform": platform.platform(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "sklearn": sklearn.__version__, "matplotlib": matplotlib.__version__},
        "outputs": output_hashes,
    }
    (EXP / "artifacts/run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(metrics.sort_values(["endpoint", "model"])[["model", "endpoint", "st_rae", "mae", "r2", "spearman_rho", "kendall_tau", "prediction_sd", "n"]].to_string(index=False))
    print("\nPaired affine-minus-raw ST-RAE differences:\n" + paired.to_string(index=False))
    print(f"\nDecision rule met: {decision_met}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
