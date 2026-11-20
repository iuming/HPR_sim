# Author: Ming Liu (刘铭)
#         National Synchrotron Radiation Laboratory (NSRL),
#         University of Science and Technology of China (USTC)
# Contact: mliu1018@ustc.edu.cn
# SPDX-License-Identifier: MIT

"""Tests for configuration parsing / serialisation."""
import numpy as np
import pytest

from hpr_sim.config import (
    JetParams,
    MotionParams,
    NozzleParams,
    SimulationParams,
    WaterParams,
    from_dict,
    parse_config_dict,
    to_dict,
)


def test_roundtrip_all_sections():
    for cls in (WaterParams, JetParams, NozzleParams, MotionParams, SimulationParams):
        obj = cls()
        d = to_dict(obj)
        assert isinstance(d, dict) and d
        back = from_dict(cls, d)
        assert back == obj


def test_unknown_keys_ignored():
    p = from_dict(SimulationParams, {"clean_cos_thresh": 0.7, "nonsense": 1})
    assert p.clean_cos_thresh == 0.7


def test_auto_sentinel_survives():
    m = from_dict(MotionParams, {"x_start_mm": "auto", "trans_speed_mm_s": "20"})
    assert m.x_start_mm == "auto"
    assert m.trans_speed_mm_s == 20.0


def test_parse_config_dict():
    cfg = {
        "jet": {"pressure_bar": "120"},
        "sim": {"clean_cos_thresh": 0.6, "dmax_kd": 26},
        "nozzle": {"n_holes": 8},
    }
    p = parse_config_dict(cfg)
    assert p["jet"].pressure_bar == 120.0
    assert p["sim"].clean_cos_thresh == 0.6
    assert p["sim"].dmax_kd == 26.0
    assert p["nozzle"].n_holes == 8
    assert isinstance(p["water"], WaterParams)


def test_water_properties():
    w = WaterParams()
    assert w.rho_mm3 == pytest.approx(1000.0 * 1e-9)
