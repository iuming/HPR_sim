"""Tests for the physics models (droplet drag, cleaning criteria)."""
import numpy as np
import pytest

from hpr_sim.config import JetParams, SimulationParams, WaterParams
from hpr_sim.physics import (
    schiller_naumann_cd,
    morsi_alexander_cd,
    putnam_cd,
    droplet_diameter_m,
    breakup_length_m,
    velocity_table,
    jet_velocity_at_distance,
    impact_pressure,
    cleaning_metric,
    efold_distance_mm,
)


def test_schiller_naumann_limits():
    # Stokes limit: Cd -> 24/Re for Re -> 0
    assert schiller_naumann_cd(0.001) == pytest.approx(24e3, rel=0.01)
    assert schiller_naumann_cd(1) == pytest.approx(27.6, rel=0.05)
    assert schiller_naumann_cd(100) == pytest.approx(1.092, rel=0.05)
    assert schiller_naumann_cd(5000) == pytest.approx(0.44, rel=0.01)
    assert schiller_naumann_cd(np.array([1.0, 100.0])).shape == (2,)


def test_morsi_alexander_ranges_positive():
    Re = np.linspace(0.05, 3000, 500)
    cd = morsi_alexander_cd(Re)
    assert (cd > 0).all()
    # Stokes limit at very low Re
    assert morsi_alexander_cd(0.01) == pytest.approx(24 / 0.01, rel=0.02)
    # Newton plateau at very high Re
    assert morsi_alexander_cd(1e5) == pytest.approx(0.44, rel=0.02)


def test_putnam():
    assert putnam_cd(1) == pytest.approx(28.0, rel=0.02)
    assert putnam_cd(100) == pytest.approx(1.101, rel=0.02)
    assert putnam_cd(5000) == pytest.approx(0.44, rel=0.01)


def test_reference_case_parameters():
    w, jp = WaterParams(), JetParams()
    v0 = jp.exit_velocity(w)
    assert v0 == pytest.approx(99.0, rel=0.01)
    assert droplet_diameter_m(v0, jp, w) * 1e6 == pytest.approx(74.5, rel=0.02)
    assert breakup_length_m(v0, jp, w) * 1e3 == pytest.approx(65.0, rel=0.02)


def test_velocity_table_monotonic_and_core():
    w, jp = WaterParams(), JetParams(drop_cd_mode="schiller-naumann")
    v0 = jp.exit_velocity(w)
    d, v = velocity_table(v0, jp, w, dmax_mm=400.0)
    assert v[0] == approx(v0)
    assert (np.diff(v) <= 1e-9).all()          # non-increasing
    Lb = breakup_length_m(v0, jp, water=w) * 1e3
    assert v[np.searchsorted(d, Lb / 2)] == approx(v0)   # intact core
    assert v[-1] < 0.5 * v0                     # decays far away


def test_constant_cd_matches_exponential():
    w, jp = WaterParams(), JetParams(drop_cd_mode="constant", drop_cd=0.5)
    v0 = jp.exit_velocity(w)
    d, v = velocity_table(v0, jp, w)
    Lb = breakup_length_m(v0, jp, water=w)
    k = 6.04
    expected = v0 * np.exp(-k * np.maximum(d * 1e-3 - Lb, 0.0))
    assert np.allclose(v, expected, rtol=1e-3)


def test_correlations_agree():
    w, jp_sn = WaterParams(), JetParams(drop_cd_mode="schiller-naumann")
    v0 = jp_sn.exit_velocity(w)
    d, v_sn = velocity_table(v0, jp_sn, w)
    for mode in ("morsi-alexander", "putnam"):
        jp = JetParams(drop_cd_mode=mode)
        _, v = velocity_table(v0, jp, w)
        assert np.allclose(v, v_sn, atol=3.0)


def test_efold_reasonable():
    w, jp = WaterParams(), JetParams(drop_cd_mode="schiller-naumann")
    e = efold_distance_mm(jp.exit_velocity(w), jp, w)
    assert 80 < e < 200
    w2, jp2 = WaterParams(), JetParams(drop_cd_mode="constant", drop_cd=0.5)
    assert efold_distance_mm(jp2.exit_velocity(w), jp2, w2) == pytest.approx(166, rel=0.05)


def test_impact_pressure():
    assert impact_pressure(99.0, 1000.0) / 1e5 == pytest.approx(49.0, rel=0.01)


def test_cleaning_metric_criteria():
    sp = SimulationParams(clean_cos_thresh=0.5, d_min_mm=5, d_max_mm=50)
    intens = np.array([0.9, 0.3, 0.9, 0.9, 0.9])
    dist = np.array([10.0, 10.0, 60.0, 3.0, 10.0])
    impact = np.array([30.0, 30.0, 30.0, 30.0, 0.5])
    eff = cleaning_metric(intens, dist, impact, sp)
    assert eff.tolist() == [True, False, False, False, True]

    # min impact criterion
    sp2 = SimulationParams(min_impact_bar=1.0)
    assert cleaning_metric(np.array([0.9]), np.array([10.0]),
                           np.array([0.5]), sp2)[0] == False
    assert cleaning_metric(np.array([0.9]), np.array([10.0]),
                           np.array([1.5]), sp2)[0] == True

    # Guha stand-off override
    sp3 = SimulationParams(d_max_mm=200.0)
    assert cleaning_metric(np.array([0.9]), np.array([100.0]),
                           np.array([5.0]), sp3,
                           d_max_eff=20.0)[0] == False


def approx(x):
    return pytest.approx(x, abs=1e-6)
