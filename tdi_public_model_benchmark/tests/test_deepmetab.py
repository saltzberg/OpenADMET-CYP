import pytest

from tdi_public_model_benchmark.deepmetab import map_checkpoint_state_dict


def test_deepmetab_readout_keys_map_to_chemprop_ffn():
    state = {
        "encoder.encoder.0.W_i.weight": "encoder",
        "readout.1.weight": "hidden",
        "readout.4.bias": "output",
    }
    assert map_checkpoint_state_dict(state) == {
        "encoder.encoder.0.W_i.weight": "encoder",
        "ffn.1.weight": "hidden",
        "ffn.4.bias": "output",
    }


def test_deepmetab_key_mapping_rejects_collisions():
    with pytest.raises(ValueError, match="collision"):
        map_checkpoint_state_dict({"readout.1.weight": 1, "ffn.1.weight": 2})
