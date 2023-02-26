import numpy as np
from scipy.special import jn_zeros
from mode_analysis.usolver import WGFiber



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

    
def test_get_init_points_to_solve():
    wgf = WGFiber(v=2.0)
    assert wgf.get_init_points_to_solve(l=0) != None
    assert np.size(wgf.get_init_points_to_solve(l=0)) == 1
    assert np.size(wgf.get_init_points_to_solve(l=1)) == 1 # 'None' is also counted.
    assert wgf.get_init_points_to_solve(l=1) == None
    
    wgf = WGFiber(v=2.406)
    assert wgf.get_init_points_to_solve(l=0) != None
    assert wgf.get_init_points_to_solve(l=1) != None
    assert wgf.get_init_points_to_solve(l=2) == None
    
    wgf = WGFiber(v=7.0)
    assert np.size(wgf.get_init_points_to_solve(l=0)) == 2
    assert np.size(wgf.get_init_points_to_solve(l=1)) == 2
    assert np.size(wgf.get_init_points_to_solve(l=2)) == 1
    assert np.size(wgf.get_init_points_to_solve(l=3)) == 1
    assert np.size(wgf.get_init_points_to_solve(l=4)) == 1
    assert wgf.get_init_points_to_solve(l=5) == None

    wgf = WGFiber(v=7.02)
    assert np.size(wgf.get_init_points_to_solve(l=0)) == 3
    assert np.size(wgf.get_init_points_to_solve(l=1)) == 2
    assert np.size(wgf.get_init_points_to_solve(l=2)) == 2
    assert np.size(wgf.get_init_points_to_solve(l=3)) == 1
    assert np.size(wgf.get_init_points_to_solve(l=4)) == 1
    assert wgf.get_init_points_to_solve(l=5) == None
