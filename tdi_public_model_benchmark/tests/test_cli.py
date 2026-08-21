import json

import pandas as pd

from tdi_public_model_benchmark.cli import evaluate_registry
from tdi_public_model_benchmark.data import prepare_labels


def test_prepare_labels_and_registry_cli(tmp_path):
    source = tmp_path / "tdi.csv"
    pd.DataFrame(
        {
            "Molecule_Name": ["A", "B", "C", "D"],
            "SMILES": ["c1ccccc1", "Cc1ccccc1", "CC", "CCC"],
            "CYP2D6_is_TDI": [False, True, False, True],
            "CYP3A4_is_TDI": [True, False, None, None],
        }
    ).to_csv(source, index=False)
    data_dir = tmp_path / "data"
    manifest = prepare_labels(source, data_dir)
    assert manifest["endpoint_summary"]["CYP2D6"]["n_labels"] == 4
    assert manifest["endpoint_summary"]["CYP3A4"]["n_labels"] == 2

    predictions = tmp_path / "scores.csv"
    pd.DataFrame(
        {"Molecule_Name": ["A", "B", "C", "D"], "score": [0.1, 0.9, 0.2, 0.8]}
    ).to_csv(predictions, index=False)
    registry = tmp_path / "registry.csv"
    pd.DataFrame(
        [
            {
                "model": "Synthetic",
                "score_name": "fixture",
                "endpoint": "CYP2D6",
                "prediction_path": predictions.name,
                "id_column": "Molecule_Name",
                "score_column": "score",
                "score_direction": "higher",
                "status": "ready",
                "source_revision": "test-only",
                "notes": "synthetic fixture",
            }
        ]
    ).to_csv(registry, index=False)
    output = tmp_path / "output"
    run = evaluate_registry(
        data_dir / "tdi_labels.csv", registry, output,
        bootstrap_replicates=20, seed=3, confidence=0.95, overwrite=False,
    )
    assert run["comparisons"][0]["source_revision"] == "test-only"
    metrics = pd.read_csv(output / "metrics.csv")
    primary = metrics.loc[metrics["metric"] == "pr_auc_average_precision", "estimate"].item()
    assert primary == 1.0
    json.loads((output / "run_manifest.json").read_text())
