import pytest
from scipy.special import jn_zeros
from mode_analysis.usolver import WGFiber


def test_gen_eigen_eq():
    wgf = WGFiber(3)


def test_find_max_jn_zeros():
    wgf = WGFiber(v=2.0)
    assert wgf.find_max_jn_zeros(l=0) == 0
    assert wgf.find_max_jn_zeros(l=1) == 0
 
    wgf = WGFiber(v=2.6)
    assert wgf.find_max_jn_zeros(l=0) == 1
    assert wgf.find_max_jn_zeros(l=1) == 0
    
    wgf = WGFiber(v=5.5)
    assert wgf.find_max_jn_zeros(l=0) == 1
    assert wgf.find_max_jn_zeros(l=1) == 1
    assert wgf.find_max_jn_zeros(l=2) == 1

    
    



    



