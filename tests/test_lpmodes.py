from warnings import warn
import pytest
from mode_analysis.modeset import LPModes


def test_init():
    LP01 = LPModes(v=2.0)
    
