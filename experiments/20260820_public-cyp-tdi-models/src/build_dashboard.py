#!/usr/bin/env python3
"""Bundle DeepMetab TDI artifacts and build the static experiment dashboard."""
from __future__ import annotations

import hashlib
import html
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

EXPERIMENT = Path(__file__).resolve().parents[1]
PROJECT = EXPERIMENT.parents[1]
SOURCE = PROJECT / "tdi_public_model_benchmark"
ARTIFACTS = EXPERIMENT / "artifacts"

COPIES = {
    "outputs/public_models_v4/metrics.csv": "metrics.csv",
    "outputs/public_models_v4/enrichment.csv": "enrichment.csv",
    "outputs/public_models_v4/bootstrap_metrics.csv": "bootstrap_metrics.csv",
    "outputs/public_models_v4/joined_predictions.csv": "joined_predictions.csv",
    "outputs/public_models_v4/run_manifest.json": "benchmark_run_manifest.json",
    "predictions/raw/deepmetab.csv": "deepmetab_predictions.csv",
    "predictions/raw/deepmetab.failures.csv": "deepmetab_failures.csv",
    "predictions/raw/deepmetab.manifest.json": "inference_manifest.json",
    "predictions/raw/cypreact.csv": "cypreact_predictions.csv",
    "predictions/raw/cypreact.failures.csv": "cypreact_failures.csv",
    "predictions/raw/cypreact.manifest.json": "cypreact_inference_manifest.json",
    "predictions/raw/cypmol_substrate.csv": "cypmol_substrate_predictions.csv",
    "predictions/raw/cypmol_substrate.manifest.json": "cypmol_substrate_manifest.json",
    "predictions/raw/cypmol_inhibitor.csv": "cypmol_inhibitor_predictions.csv",
    "predictions/raw/cypmol_inhibitor.manifest.json": "cypmol_inhibitor_manifest.json",
    "predictions/raw/cypmol.failures.csv": "cypmol_failures.csv",
    "config/model_registry.csv": "model_registry.csv",
}
SOURCE_CODE = [
    "src/tdi_public_model_benchmark/deepmetab.py",
    "src/tdi_public_model_benchmark/cypreact.py",
    "src/tdi_public_model_benchmark/cypmol.py",
    "src/tdi_public_model_benchmark/evaluation.py",
    "src/tdi_public_model_benchmark/data.py",
    "src/tdi_public_model_benchmark/cli.py",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def metric_row(metrics: pd.DataFrame, model: str, score_name: str,
               endpoint: str, metric: str) -> pd.Series:
    rows = metrics[
        (metrics.model == model) & (metrics.score_name == score_name)
        & (metrics.endpoint == endpoint) & (metrics.metric == metric)
    ]
    if len(rows) != 1:
        raise ValueError(
            f"Expected one {model}/{score_name}/{endpoint}/{metric} row; found {len(rows)}"
        )
    return rows.iloc[0]


def scale(value: float, marker: float, label: str, marker_label: str) -> str:
    value_pct = max(0.0, min(100.0, 100.0 * value))
    marker_pct = max(0.0, min(100.0, 100.0 * marker))
    return (
        '<div class="measure" role="img" '
        f'aria-label="{html.escape(label)} {value:.3f}; {html.escape(marker_label)} {marker:.3f}">'
        f'<span class="fill" style="width:{value_pct:.3f}%"></span>'
        f'<i style="left:{marker_pct:.3f}%"></i></div>'
    )


def build_html(metrics: pd.DataFrame, enrichment: pd.DataFrame, joined: pd.DataFrame) -> str:
    endpoints = ["CYP2D6", "CYP3A4"]
    series = [
        ("CYPMol", "substrate_score", "CYPMol substrate"),
        ("CYPMol", "inhibitor_score", "CYPMol inhibitor"),
        ("DeepMetab", "substrate_score", "DeepMetab substrate"),
        ("CypReact", "cypreact_public_score", "CypReact R/N"),
    ]
    prevalence = {"CYP2D6": 324 / 1497, "CYP3A4": 764 / 3584}
    summary_rows = []
    plots = []
    enrichment_rows = []
    display_names = {(model, score): display for model, score, display in series}
    for model, score_name, display in series:
        for endpoint in endpoints:
            ap = metric_row(metrics, model, score_name, endpoint, "pr_auc_average_precision")
            roc = metric_row(metrics, model, score_name, endpoint, "roc_auc")
            coverage = metric_row(metrics, model, score_name, endpoint, "coverage")
            summary_rows.append(
                f"<tr><td>{display}</td><td>{endpoint}</td><td>{int(ap.n_labeled):,}</td>"
                f"<td>{coverage.estimate:.4f}</td>"
                f"<td>{ap.estimate:.4f} <small>{ap.ci_low:.4f}–{ap.ci_high:.4f}</small></td>"
                f"<td>{roc.estimate:.4f} <small>{roc.ci_low:.4f}–{roc.ci_high:.4f}</small></td></tr>"
            )
            plots.append(
                f'<section class="endpoint"><h3>{display} · {endpoint}</h3>'
                f'<div class="plotline"><span>Average precision</span><b>{ap.estimate:.3f}</b>'
                f'{scale(float(ap.estimate), prevalence[endpoint], "average precision", "label prevalence")}'
                f'<small>vertical mark: prevalence {prevalence[endpoint]:.3f}</small></div>'
                f'<div class="plotline"><span>ROC-AUC</span><b>{roc.estimate:.3f}</b>'
                f'{scale(float(roc.estimate), 0.5, "ROC-AUC", "random ranking")}'
                '<small>vertical mark: random ranking 0.500</small></div></section>'
            )
    for row in enrichment.sort_values(["model", "endpoint", "fraction"]).itertuples(index=False):
            hits = f"{row.hits:.2f}" if abs(row.hits - round(row.hits)) > 1e-9 else f"{int(round(row.hits))}"
            display = display_names[(row.model, row.score_name)]
            enrichment_rows.append(
                f"<tr><td>{display}</td><td>{row.endpoint}</td><td>{int(round(100 * row.fraction))}%</td>"
                f"<td>{int(row.k):,}</td><td>{hits}</td>"
                f"<td>{row.positive_rate:.4f}</td><td>{row.enrichment:.3f} "
                f"<small>{row.ci_low:.3f}–{row.ci_high:.3f}</small></td></tr>"
            )

    cypreact = joined[(joined.model == "CypReact") & joined.raw_score.notna()].copy()
    call_rows = []
    for endpoint in endpoints:
        part = cypreact[cypreact.endpoint == endpoint]
        counts = part.groupby(["is_tdi", "raw_score"]).size()
        value = lambda truth, score: int(counts.get((truth, score), 0))
        call_rows.append(
            f"<tr><td>{endpoint}</td><td>{value(1, 1.0):,}</td><td>{value(1, 0.0):,}</td>"
            f"<td>{value(0, 1.0):,}</td><td>{value(0, 0.0):,}</td></tr>"
        )

    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    return f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Public CYP models versus OpenADMET TDI</title>
<style>
:root{{--ink:#20252a;--muted:#66717a;--rule:#d9dee2;--wash:#f5f6f5;--accent:#9a4d12;--paper:#fff}}
*{{box-sizing:border-box}} body{{margin:0;background:var(--paper);color:var(--ink);font:16px/1.5 ui-serif,Georgia,serif}}
main{{width:min(68rem,calc(100% - 2rem));margin:2.8rem auto 5rem}} h1,h2,h3{{line-height:1.18}} h1{{font-size:clamp(2rem,5vw,3.15rem);letter-spacing:-.03em;margin:.15rem 0 .55rem}} h2{{margin:2.5rem 0 .9rem;border-bottom:1px solid var(--rule);padding-bottom:.35rem;font-size:1.25rem}} h3{{font:700 1rem ui-sans-serif,system-ui,sans-serif;margin:0 0 .9rem}}
p{{max-width:52rem}} .eyebrow{{font:700 .76rem ui-sans-serif,system-ui,sans-serif;letter-spacing:.09em;text-transform:uppercase;color:var(--muted)}} .meta{{font:.82rem ui-monospace,monospace;color:var(--muted)}} .lede{{font-size:1.18rem;max-width:55rem}} .notice{{border-left:.25rem solid var(--accent);background:var(--wash);padding:.7rem 1rem;font: .92rem ui-sans-serif,system-ui,sans-serif;margin:1.3rem 0}}
table{{border-collapse:collapse;width:100%;font:.9rem/1.35 ui-sans-serif,system-ui,sans-serif}} th,td{{padding:.5rem .55rem;border-bottom:1px solid var(--rule);text-align:right;vertical-align:top}} th:first-child,td:first-child{{text-align:left}} th{{font-weight:650;color:#39434b}} small{{color:var(--muted);font-weight:400}} .endpoint-grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:2rem;margin:1.2rem 0}} .endpoint{{min-width:0}} .plotline{{display:grid;grid-template-columns:1fr auto;gap:.18rem .8rem;margin:.8rem 0 1.15rem;font: .87rem ui-sans-serif,system-ui,sans-serif}} .measure{{grid-column:1/-1;position:relative;height:.55rem;background:#ecefed;border-radius:1px;overflow:visible}} .fill{{display:block;height:100%;background:var(--accent)}} .measure i{{position:absolute;top:-.22rem;height:.99rem;border-left:1px solid #30373d}} .plotline small{{grid-column:1/-1}} ol.procedure{{max-width:54rem;padding-left:1.35rem}} ol.procedure li{{margin:.65rem 0;padding-left:.25rem}} .status td:nth-child(2){{text-align:left}} a{{color:#285b84}} code{{font:.88em ui-monospace,monospace;overflow-wrap:anywhere}} .links{{font:.9rem ui-sans-serif,system-ui,sans-serif}} @media(max-width:720px){{.endpoint-grid{{grid-template-columns:1fr}} table{{font-size:.8rem}} th,td{{padding:.4rem .3rem}}}}
</style></head><body><main>
<div class="eyebrow">OpenADMET CYP · experiment result</div>
<h1>Public CYP models versus TDI</h1>
<div class="meta">20260820_public-cyp-tdi-models · generated {generated}</div>
<p class="lede">Do public CYP substrate, inhibitor or reactant scores rank time-dependent inhibitors? The completed comparisons say no. CYPMol adds no useful positive TDI enrichment and is often inversely associated.</p>
<div class="notice"><strong>Experiment only — not a submission.</strong> No TDI labels were used to fit or tune the public model.</div>
<h2>Result</h2>
<table><thead><tr><th>Model</th><th>Endpoint</th><th>Labels</th><th>Coverage</th><th>Average precision <small>(95% CI)</small></th><th>ROC-AUC <small>(95% CI)</small></th></tr></thead><tbody>{''.join(summary_rows)}</tbody></table>
<div class="endpoint-grid">{''.join(plots)}</div>
<p>Average precision should be compared with label prevalence: 0.216 for CYP2D6 and 0.213 for CYP3A4. None of the completed scores supplies useful positive TDI discrimination. CYPMol CYP3A4 inhibitor has a slightly positive global ROC-AUC but depletes positives in its highest-scoring tail.</p>
<h2>Top-ranked enrichment</h2>
<table><thead><tr><th>Model</th><th>Endpoint</th><th>Top fraction</th><th>Molecules</th><th>Expected TDI positives</th><th>Positive rate</th><th>Enrichment <small>(95% CI)</small></th></tr></thead><tbody>{''.join(enrichment_rows)}</tbody></table>
<p>Expected counts are fractional when a top-k boundary cuts a score tie. This prevents arbitrary molecule-ID ordering from creating apparent enrichment for CypReact's binary R/N score.</p>
<h2>CypReact call saturation</h2>
<table><thead><tr><th>Endpoint</th><th>TDI+ called R</th><th>TDI+ called N</th><th>TDI− called R</th><th>TDI− called N</th></tr></thead><tbody>{''.join(call_rows)}</tbody></table>
<p>CypReact called 96.3% of scored CYP2D6-labelled molecules and 99.4% of scored CYP3A4-labelled molecules reactants. High recall therefore comes with almost no specificity.</p>
<h2>Procedure</h2>
<ol class="procedure"><li><strong>Labels.</strong> Build CYP2D6 and CYP3A4 TDI label tables from the released training data: 1,497 and 3,584 labelled molecules respectively.</li><li><strong>DeepMetab.</strong> Score the 4,822-molecule union using original supplied SMILES and five public multitask checkpoints. Map released <code>readout.*</code> head keys to ChemProp <code>ffn.*</code> and require a strict full-state load.</li><li><strong>CypReact.</strong> Convert original SMILES to an ID-labelled SDF and run the unmodified v1.2 JAR. Encode its faithful public calls as N=0/R=1.</li><li><strong>CYPMol.</strong> Apply the authors' deterministic RDKit 3D preparation, then load ten complete substrate and ten inhibitor checkpoints one at a time on the GPU. Reuse fixed endpoint protein contexts while preserving every checkpoint score, ensemble mean and variance. Two molecules failed conformer generation and remain missing.</li><li><strong>Evaluation.</strong> Report average precision, ROC-AUC, tie-aware top 5%/10%/20% enrichment and coverage without fitting or threshold tuning.</li><li><strong>Uncertainty.</strong> Use 2,000 percentile bootstrap replicates clustered by Bemis–Murcko scaffold over the complete labelled endpoint universe.</li></ol>
<h2>Model queue</h2>
<table class="status"><thead><tr><th>Model score</th><th>Status and next action</th></tr></thead><tbody><tr><td>DeepMetab substrate</td><td><strong>Completed.</strong> Retain as a negative-control external comparator.</td></tr><tr><td>CypReact reactant</td><td><strong>Completed.</strong> Public R/N output is too saturated for TDI prioritization.</td></tr><tr><td>CYPMol substrate + inhibitor</td><td><strong>Completed.</strong> Neither score family provides useful positive TDI prioritization.</td></tr><tr><td>DeepP450 substrate</td><td><strong>Deferred.</strong> The 58.8 GB Baidu bundle requires authenticated download; the adapter is prepared.</td></tr></tbody></table>
<h2>Artifacts</h2>
<p class="links"><a href="artifacts/metrics.csv">metrics</a> · <a href="artifacts/enrichment.csv">enrichment</a> · <a href="artifacts/bootstrap_metrics.csv">bootstrap replicates</a> · <a href="artifacts/joined_predictions.csv">joined predictions</a> · <a href="artifacts/deepmetab_predictions.csv">DeepMetab predictions</a> · <a href="artifacts/cypreact_predictions.csv">CypReact predictions</a> · <a href="artifacts/cypmol_substrate_predictions.csv">CYPMol substrate</a> · <a href="artifacts/cypmol_inhibitor_predictions.csv">CYPMol inhibitor</a> · <a href="artifacts/inference_manifest.json">DeepMetab manifest</a> · <a href="artifacts/cypreact_inference_manifest.json">CypReact manifest</a> · <a href="artifacts/cypmol_substrate_manifest.json">CYPMol substrate manifest</a> · <a href="artifacts/cypmol_inhibitor_manifest.json">CYPMol inhibitor manifest</a> · <a href="artifacts/benchmark_run_manifest.json">benchmark manifest</a> · <a href="artifacts/run_manifest.json">experiment manifest</a></p>
</main></body></html>'''


def main() -> int:
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    source_hashes = {}
    copied_hashes = {}
    for source_relative, target_name in COPIES.items():
        source_path = SOURCE / source_relative
        if not source_path.is_file():
            raise FileNotFoundError(source_path)
        target_path = ARTIFACTS / target_name
        shutil.copy2(source_path, target_path)
        source_hashes[f"tdi_public_model_benchmark/{source_relative}"] = sha256_file(source_path)
        copied_hashes[f"artifacts/{target_name}"] = sha256_file(target_path)
    for source_relative in SOURCE_CODE:
        path = SOURCE / source_relative
        source_hashes[f"tdi_public_model_benchmark/{source_relative}"] = sha256_file(path)

    metrics = pd.read_csv(ARTIFACTS / "metrics.csv")
    enrichment = pd.read_csv(ARTIFACTS / "enrichment.csv")
    joined = pd.read_csv(ARTIFACTS / "joined_predictions.csv")
    dashboard = build_html(metrics, enrichment, joined)
    index_path = EXPERIMENT / "index.html"
    index_path.write_text(dashboard, encoding="utf-8")

    manifest = {
        "status": "experiment_only_not_submitted",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "experiment_id": "20260820_public-cyp-tdi-models",
        "completed_models": [
            "DeepMetab substrate score", "CypReact public reactant score",
            "CYPMol substrate score", "CYPMol inhibitor score",
        ],
        "pending_models": [],
        "deferred_models": ["DeepP450 substrate"],
        "bootstrap_replicates": 2000,
        "bootstrap_seed": 20260820,
        "source_files": source_hashes,
        "outputs": {**copied_hashes, "index.html": sha256_file(index_path)},
    }
    (ARTIFACTS / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"experiment": manifest["experiment_id"], "outputs": len(manifest["outputs"])}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
