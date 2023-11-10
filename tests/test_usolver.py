import pytest
import numpy as np
from scipy.special import jn_zeros
from ..usolver import WGFiber


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
    assert np.size(wgf.get_init_points_to_solve(l=1)) == 1  # 'None' is also counted.
    assert wgf.get_init_points_to_solve(l=1) == None

    wgf = WGFiber(v=2.406)
    assert wgf.get_init_points_to_solve(l=0) != None
    assert wgf.get_init_points_to_solve(l=1) != None
    assert wgf.get_init_points_to_solve(l=2) == None

    wgf = WGFiber(v=7.0)
    assert any(x != None for x in wgf.get_init_points_to_solve(l=0))
    assert np.size(wgf.get_init_points_to_solve(l=0)) == 2
    assert np.size(wgf.get_init_points_to_solve(l=1)) == 2
    assert np.size(wgf.get_init_points_to_solve(l=2)) == 1
    assert np.size(wgf.get_init_points_to_solve(l=3)) == 1
    assert np.size(wgf.get_init_points_to_solve(l=4)) == 1
    assert wgf.get_init_points_to_solve(l=5) == None

    wgf = WGFiber(v=7.02)
    assert any(x != None for x in wgf.get_init_points_to_solve(l=0))
    assert np.size(wgf.get_init_points_to_solve(l=0)) == 3
    assert np.size(wgf.get_init_points_to_solve(l=1)) == 2
    assert np.size(wgf.get_init_points_to_solve(l=2)) == 2
    assert np.size(wgf.get_init_points_to_solve(l=3)) == 1
    assert np.size(wgf.get_init_points_to_solve(l=4)) == 1
    assert wgf.get_init_points_to_solve(l=5) == None


def test_get_roots_for_u():
    wgf = WGFiber(v=2.0)
    roots = wgf.get_roots_for_u(l=0)
    assert roots != None
    assert np.isclose(roots[0], 1.52818403)
    assert wgf.get_roots_for_u(l=1) == None

    wgf = WGFiber(7.02)
    roots = wgf.get_roots_for_u(l=0)
    assert np.isclose(roots[0], 2.10071035)
    assert np.isclose(roots[1], 4.77191413)
    assert np.isclose(roots[-1], wgf.get_init_points_to_solve(l=0)[-1])

    wgf = WGFiber(6.380)
    assert np.size(wgf.get_roots_for_u(l=0)) == 2
    assert np.size(wgf.get_roots_for_u(l=4)) == 1

    wgf = WGFiber(7.016)
    assert np.size(wgf.get_roots_for_u(l=0)) == 3
    assert np.size(wgf.get_roots_for_u(l=2)) == 2

    wgf = WGFiber(8.4174)
    assert np.size(wgf.get_roots_for_u(l=0)) == 3
    assert np.size(wgf.get_roots_for_u(l=1)) == 2
    assert np.size(wgf.get_roots_for_u(l=2)) == 2
    assert np.size(wgf.get_roots_for_u(l=3)) == 2
    assert np.size(wgf.get_roots_for_u(l=4)) == 1
    assert np.size(wgf.get_roots_for_u(l=5)) == 1


def test_get_cutoff_value():
    assert np.isclose(WGFiber.get_lp_cutoff(0, 1), 0.0)
    assert np.isclose(WGFiber.get_lp_cutoff(1, 1), 2.404825557)
    assert np.isclose(WGFiber.get_lp_cutoff(2, 1), 3.831705970)
    assert np.isclose(WGFiber.get_lp_cutoff(0, 2), 3.831705970)
    assert np.isclose(WGFiber.get_lp_cutoff(3, 1), 5.135622301)
    assert np.isclose(WGFiber.get_lp_cutoff(1, 2), 5.520078110)


def test_get_all_uset():
    wgf = WGFiber(2.0)
    assert wgf.uset["u01"] == wgf.get_roots_for_u(0)[0]
    with pytest.raises(KeyError):
        wgf.uset["u11"]
        wgf.uset["u21"]


def test_u_lm():
    wgf = WGFiber(v=2.0)
    assert wgf.u_lm(1, 1) == None
    with pytest.raises(ValueError):
        wgf.u_lm(-1, 1)
    with pytest.raises(ValueError):
        wgf.u_lm(1, 0)
    with pytest.raises(ValueError):
        wgf.u_lm(2.1, 1)
