from warnings import warn
from mode_analysis.mtool import LPModes
import numpy as np


def test_init():
    lps = LPModes(v=2.0, dia_core=10)
    assert lps.uset["u01"] == lps.u_lm(0, 1)


def test_field_eq():
    v = 3
    d = 10
    lps = LPModes(v, d)
    l = 1
    m = 1
    u = lps.u_lm(l, m)
    fcore = lps.field_core_eq(l)
    fclad = lps.field_clad_eq(l)
    assert np.isclose(fcore(u, d / 2, 0), fclad(u, d / 2, 0))

    v = 4
    d = 10
    lps = LPModes(v, d)
    l = 2
    m = 1
    u = lps.u_lm(l, m)
    fcore = lps.field_core_eq(l)
    fclad = lps.field_clad_eq(l)
    assert np.isclose(fcore(u, d / 2, np.pi / 2), fclad(u, d / 2, np.pi / 2))
