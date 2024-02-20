import numpy as np
import pytest as pt

import paramo as pp

class TestSpecRelaToolkit(object):
    def test_v_rela_s(self) -> None:
        test_beta = 0.1
        expected_gamma = 1.005
        assert pp.gofb(test_beta) == pt.approx(expected_gamma, rel=1e-3)

    def test_v_rela_v(self) -> None:
        test_beta = np.array([0.1, 0.25, 0.5])
        expected_gamma = np.array([1.005, 1.033, 1.155])
        assert pp.gofb(test_beta) == pt.approx(expected_gamma, rel=1e-3)

    def test_lorentz_f_s(self) -> None:
        expected_beta = 0.1
        test_gamma = 1.005
        assert pp.bofg(test_gamma) == pt.approx(expected_beta, rel=1e-2)

    def test_lorentz_f_v(self) -> None:
        expected_beta = np.array([0.1, 0.25, 0.5])
        test_gamma = np.array([1.005, 1.033, 1.155])
        assert pp.bofg(test_gamma) == pt.approx(expected_beta, rel=1e-2)
