import math
import numpy as np
import pytest
import tleco


class TestBofg:
    """bofg: compute beta = v/c from Lorentz factor gamma."""

    def test_rest(self):
        # gamma=1 => beta=0 (particle at rest)
        assert tleco.bofg(1.0) == pytest.approx(0.0, abs=1e-15)

    def test_known_scalar(self):
        # gamma=2 => beta = sqrt(3)/2
        expected = math.sqrt(3.0) / 2.0
        assert tleco.bofg(2.0) == pytest.approx(expected, rel=1e-6)

    def test_ultrarelativistic(self):
        # large gamma => beta approaches 1
        beta = tleco.bofg(1e6)
        assert beta == pytest.approx(1.0, rel=1e-10)
        assert beta < 1.0

    def test_vector(self):
        gammas = [1.0, 2.0, 10.0]
        betas = tleco.bofg(gammas)
        assert betas[0] == pytest.approx(0.0, abs=1e-15)
        assert betas[1] == pytest.approx(math.sqrt(3.0) / 2.0, rel=1e-6)
        assert betas[2] == pytest.approx(math.sqrt(1.0 - 1.0 / 100.0), rel=1e-6)

    def test_inverse_of_gofb(self):
        # bofg(gofb(b)) == b for a range of beta values
        betas = [0.1, 0.25, 0.5, 0.9, 0.99]
        for b in betas:
            g = tleco.gofb(b)
            b_back = tleco.bofg(g)
            assert b_back == pytest.approx(b, rel=1e-6)


class TestGofb:
    """gofb: compute Lorentz factor gamma from beta = v/c."""

    def test_rest(self):
        # beta=0 => gamma=1
        assert tleco.gofb(0.0) == pytest.approx(1.0, rel=1e-10)

    def test_known_scalar(self):
        # beta=0.1 => gamma ~ 1.005037
        expected = 1.0 / math.sqrt(1.0 - 0.1**2)
        assert tleco.gofb(0.1) == pytest.approx(expected, rel=1e-3)

    def test_vector(self):
        betas = [0.1, 0.25, 0.5]
        gammas = tleco.gofb(betas)
        for b, g in zip(betas, gammas):
            expected = 1.0 / math.sqrt(1.0 - b**2)
            assert g == pytest.approx(expected, rel=1e-3)

    def test_inverse_of_bofg(self):
        # gofb(bofg(g)) == g for a range of gamma values
        gammas = [1.0, 2.0, 5.0, 100.0]
        for g in gammas:
            b = tleco.bofg(g)
            g_back = tleco.gofb(b)
            assert g_back == pytest.approx(g, rel=1e-6)
