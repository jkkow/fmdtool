from warnings import warn
from mode_analysis.mtool import LPModes


def test_init():
    lps = LPModes(v=2.0, dia_core=10)
    assert lps.u_set["u01"] == lps.u_lm(0,1)

    
