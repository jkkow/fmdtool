import pytest
from mode_analysis.usolver import WGFiber


@pytest.fixture
def wgf():
    return WGFiber(5)


def test_gen_eigen_eq():
    wgf.gen_eigen_eq()


if __name__ == "__main__":
    wgf = WGFiber()
