from pathlib import Path

import pandas as pd
import pytest

from tdi_public_model_benchmark.cypreact import parse_cypreact_csv


def test_parse_cypreact_public_calls(tmp_path: Path):
    path = tmp_path / "output.csv"
    path.write_text(
        "Inchiky, SMILES, Title, 1A2, 2B6, 2A6, 2C8, 2C9, 2C19, 2D6, 2|E1, 3A4\n"
        "KEY,CCO,M1,null,null,null,null,null,null,R,null,N\n"
    )
    parsed = parse_cypreact_csv(path)
    assert parsed.loc[0, "Molecule_Name"] == "M1"
    assert parsed.loc[0, "CYP2D6_cypreact_call"] == "R"
    assert parsed.loc[0, "CYP2D6_cypreact_public_score"] == 1.0
    assert parsed.loc[0, "CYP3A4_cypreact_public_score"] == 0.0


def test_parse_cypreact_rejects_unknown_calls(tmp_path: Path):
    path = tmp_path / "output.csv"
    pd.DataFrame({"Title": ["M1"], "2D6": ["maybe"], "3A4": ["R"]}).to_csv(
        path, index=False
    )
    with pytest.raises(ValueError, match="Unexpected CYP2D6"):
        parse_cypreact_csv(path)
