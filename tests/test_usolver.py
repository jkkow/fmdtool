import pytest
from mode_analysis.usolver import WGFiber


def test_make_instance():
    wgf = WGFiber()
    wgf.V = 10
    assert wgf.V == 10
    

if __name__ == "__main__":
    wgf = WGFiber()
