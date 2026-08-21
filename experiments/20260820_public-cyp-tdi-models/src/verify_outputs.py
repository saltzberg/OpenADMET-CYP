#!/usr/bin/env python3
"""Verify the public CYP versus TDI experiment bundle and dashboard."""
from __future__ import annotations

import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path

import numpy as np
import pandas as pd

EXPERIMENT = Path(__file__).resolve().parents[1]
PROJECT = EXPERIMENT.parents[1]
ARTIFACTS = EXPERIMENT / "artifacts"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class DashboardParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.counts = {"title": 0, "main": 0, "h1": 0, "h2": 0, "table": 0}
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag in self.counts:
            self.counts[tag] += 1
        if tag == "a":
            href = dict(attrs).get("href")
            if href:
                self.links.append(str(href))


def main() -> int:
    manifest = json.loads((ARTIFACTS / "run_manifest.json").read_text())
    assert manifest["status"] == "experiment_only_not_submitted"
    assert manifest["experiment_id"] == "20260820_public-cyp-tdi-models"
    assert manifest["bootstrap_replicates"] == 2000
    for relative, expected in manifest["source_files"].items():
        assert sha256_file(PROJECT / relative) == expected, relative
    for relative, expected in manifest["outputs"].items():
        assert sha256_file(EXPERIMENT / relative) == expected, relative

    predictions = pd.read_csv(ARTIFACTS / "deepmetab_predictions.csv")
    failures = pd.read_csv(ARTIFACTS / "deepmetab_failures.csv")
    assert len(predictions) == 4822 and predictions["Molecule_Name"].is_unique
    assert len(failures) == 0
    score_columns = [c for c in predictions if c.endswith("_substrate_score")]
    checkpoint_columns = [c for c in predictions if "_model_" in c]
    assert len(score_columns) == 9 and len(checkpoint_columns) == 45
    assert np.isfinite(predictions[score_columns + checkpoint_columns]).all().all()

    cypreact = pd.read_csv(ARTIFACTS / "cypreact_predictions.csv")
    cypreact_failures = pd.read_csv(ARTIFACTS / "cypreact_failures.csv")
    assert len(cypreact) == 4822 and cypreact["Molecule_Name"].is_unique
    assert len(cypreact_failures) == 1
    assert cypreact_failures.iloc[0]["Molecule_Name"] == "OCNT-2328840"
    cypreact_scores = ["CYP2D6_cypreact_public_score", "CYP3A4_cypreact_public_score"]
    assert cypreact[cypreact_scores].notna().all(axis=1).sum() == 4821
    for column in cypreact_scores:
        assert set(cypreact[column].dropna().unique()) <= {0.0, 1.0}

    cypmol_substrate = pd.read_csv(ARTIFACTS / "cypmol_substrate_predictions.csv")
    cypmol_inhibitor = pd.read_csv(ARTIFACTS / "cypmol_inhibitor_predictions.csv")
    cypmol_failures = pd.read_csv(ARTIFACTS / "cypmol_failures.csv")
    assert len(cypmol_substrate) == len(cypmol_inhibitor) == 4822
    assert cypmol_substrate["Molecule_Name"].is_unique
    assert cypmol_inhibitor["Molecule_Name"].is_unique
    assert set(cypmol_failures["Molecule_Name"]) == {"OCNT-0453746", "OCNT-0453782"}
    for task, frame in [("substrate", cypmol_substrate), ("inhibitor", cypmol_inhibitor)]:
        means = [f"CYP2D6_{task}_score", f"CYP3A4_{task}_score"]
        checkpoint_scores = [
            column for column in frame
            if column.startswith((f"CYP2D6_{task}_", f"CYP3A4_{task}_"))
            and not column.endswith(("_score", "_variance"))
        ]
        assert len(checkpoint_scores) == 20
        assert frame[means].notna().all(axis=1).sum() == 4820
        assert np.isfinite(frame[means].dropna()).all().all()

    metrics = pd.read_csv(ARTIFACTS / "metrics.csv")
    enrichment = pd.read_csv(ARTIFACTS / "enrichment.csv")
    bootstrap = pd.read_csv(ARTIFACTS / "bootstrap_metrics.csv")
    joined = pd.read_csv(ARTIFACTS / "joined_predictions.csv")
    assert len(metrics) == 48 and set(metrics.endpoint) == {"CYP2D6", "CYP3A4"}
    assert set(metrics.model) == {"CYPMol", "DeepMetab", "CypReact"}
    assert set(metrics.metric) == {
        "coverage", "pr_auc_average_precision", "roc_auc",
        "enrichment_top_5pct", "enrichment_top_10pct", "enrichment_top_20pct",
    }
    assert (metrics.valid_bootstrap_replicates == 2000).all()
    coverage = metrics[metrics.metric == "coverage"].set_index(
        ["model", "score_name", "endpoint"]
    )["estimate"]
    assert coverage[("DeepMetab", "substrate_score", "CYP2D6")] == 1
    assert coverage[("DeepMetab", "substrate_score", "CYP3A4")] == 1
    assert coverage[("CypReact", "cypreact_public_score", "CYP2D6")] == 1
    assert np.isclose(
        coverage[("CypReact", "cypreact_public_score", "CYP3A4")], 3583 / 3584
    )
    for score_name in ["substrate_score", "inhibitor_score"]:
        assert coverage[("CYPMol", score_name, "CYP2D6")] == 1
        assert np.isclose(coverage[("CYPMol", score_name, "CYP3A4")], 3582 / 3584)
    assert len(enrichment) == 24
    assert {"boundary_tie_size", "boundary_fraction"} <= set(enrichment.columns)
    assert len(bootstrap) == 96000 and bootstrap.replicate.nunique() == 2000
    assert len(joined) == 4 * (1497 + 3584)

    expected = {
        ("DeepMetab", "CYP2D6", "pr_auc_average_precision"): 0.21396990966433738,
        ("DeepMetab", "CYP2D6", "roc_auc"): 0.49168008588298445,
        ("DeepMetab", "CYP3A4", "pr_auc_average_precision"): 0.19739242872321103,
        ("DeepMetab", "CYP3A4", "roc_auc"): 0.4634723924102335,
        ("CypReact", "CYP2D6", "pr_auc_average_precision"): 0.21907323834154843,
        ("CypReact", "CYP2D6", "roc_auc"): 0.5076884215844147,
        ("CypReact", "CYP3A4", "pr_auc_average_precision"): 0.21392294036536807,
        ("CypReact", "CYP3A4", "roc_auc"): 0.5020610888343681,
        ("CYPMol", "CYP2D6", "pr_auc_average_precision", "substrate_score"): 0.22971253654677193,
        ("CYPMol", "CYP2D6", "roc_auc", "substrate_score"): 0.5104143643501416,
        ("CYPMol", "CYP2D6", "pr_auc_average_precision", "inhibitor_score"): 0.19620338674890286,
        ("CYPMol", "CYP2D6", "roc_auc", "inhibitor_score"): 0.4687069137907444,
        ("CYPMol", "CYP3A4", "pr_auc_average_precision", "substrate_score"): 0.18651665205498452,
        ("CYPMol", "CYP3A4", "roc_auc", "substrate_score"): 0.4772872316707479,
        ("CYPMol", "CYP3A4", "pr_auc_average_precision", "inhibitor_score"): 0.20686599265605818,
        ("CYPMol", "CYP3A4", "roc_auc", "inhibitor_score"): 0.5321349477368748,
    }
    for key, value in expected.items():
        score_name = key[3] if len(key) == 4 else (
            "substrate_score" if key[0] == "DeepMetab" else "cypreact_public_score"
        )
        observed = metrics.loc[
            (metrics.model == key[0]) & (metrics.score_name == score_name)
            & (metrics.endpoint == key[1]) & (metrics.metric == key[2]),
            "estimate",
        ].item()
        assert np.isclose(observed, value, rtol=0, atol=1e-14), key

    document = (EXPERIMENT / "index.html").read_text(encoding="utf-8")
    parser = DashboardParser()
    parser.feed(document)
    parser.close()
    assert parser.counts == {"title": 1, "main": 1, "h1": 1, "h2": 6, "table": 4}
    assert "not a submission" in document.lower()
    assert "0.2140" in document and "0.1974" in document
    assert "0.2191" in document and "0.2139" in document
    assert "0.2297" in document and "0.2069" in document
    assert "readout.*" in document and "ffn.*" in document
    assert "CypReact call saturation" in document
    assert "CYPMol substrate" in document and "CYPMol inhibitor" in document
    for link in parser.links:
        if link.startswith(("http://", "https://", "#")):
            continue
        assert (EXPERIMENT / link).is_file(), link

    print("source and bundle hashes: PASS")
    print("DeepMetab/CypReact/CYPMol prediction and failure contracts: PASS")
    print("combined metrics, tie-aware enrichment and 2,000-replicate bootstrap: PASS")
    print("dashboard structure, values and links: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
