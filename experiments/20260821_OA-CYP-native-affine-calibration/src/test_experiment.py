#!/usr/bin/env python3
"""Focused unit tests for affine fitting and ST-RAE/bootstrap semantics."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_experiment as R


def main() -> int:
    intercept, slope = R.fit_affine([0, 1, 2, 3], [2, 5, 8, 11])
    assert np.isclose(intercept, 2.0)
    assert np.isclose(slope, 3.0)
    assert np.allclose(R.soft_error([1, 2, 4], [1.5, 1.5, 2.0], [2.5, 2.5, 3.0]), [0.5, 0.0, 1.0])
    assert np.isclose(R.st_rae([2, 4], [1, 5], [1.5, 3.5], [2.5, 4.5]), 0.5)
    rows = []
    for endpoint_offset, endpoint in enumerate(R.ENDPOINTS):
        for i in range(8):
            observed = float(i + endpoint_offset)
            rows.append({
                "compound_id": f"c{i}", "endpoint": endpoint, "observed": observed,
                "conf_low": observed - 0.1, "conf_high": observed + 0.1,
                "native_prediction": observed + 1.0,
                "calibrated_prediction": observed + 0.2,
                "scaffold_id": f"s{i // 2}", "fold": i % 2,
            })
    replicate, summary = R.paired_group_bootstrap(pd.DataFrame(rows), n_replicates=25, seed=7)
    assert len(replicate) == 25 * 5
    assert len(summary) == 5
    assert (summary.replicates == 25).all()
    assert (summary.ci_high < 0).all()
    print("PASS affine OLS, soft error/ST-RAE, paired grouped bootstrap")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
