import pytest
from ..mtool import LPModes
import numpy as np


def test_init():
    lps = LPModes(v=2.0, d=10)
    assert lps.uset["u01"] == lps.u_lm(0, 1)


def test_field_eq():
    v = 3
    d = 10
    lps = LPModes(v, d)
    l, m = 1, 1
    u11 = lps.u_lm(l, m)
    fcore = lps.field_core_eq(l, u11)
    fclad = lps.field_clad_eq(l, u11)
    assert np.isclose(fcore(d / 2, 0), fclad(d / 2, 0))

    v = 4
    d = 10
    lps = LPModes(v, d)
    l, m = 2, 1
    u21 = lps.u_lm(l, m)
    fcore = lps.field_core_eq(l, u21)
    fclad = lps.field_clad_eq(l, u21)
    assert np.isclose(fcore(d / 2, np.pi / 2), fclad(d / 2, np.pi / 2))


def test_LP_normalization():
    lps = LPModes(5, 10)
    lp11 = lps.LP(1, 1)
    lp21 = lps.LP(2, 1)
    assert np.isclose(np.sum(abs(lp11 * lp11)), 1.000)
    assert np.isclose(np.sum(abs(lp21 * lp21)), 1.000)
    assert np.isclose(np.sum(abs(lp21 * lp21)), np.sum(abs(lp11 * lp11)))


def test_LP():
    lps = LPModes(3, 10)
    with pytest.raises(ValueError):
        lps.LP(3,1)


def test_get_power():
    lps = LPModes(3, 10)
    lp01 = lps.LP(0, 1)
    lp11 = lps.LP(1, 1)
    assert np.isclose(lps.get_power(lp01), 1.0)
    assert np.isclose(lps.get_power(lp11), 1.0)
