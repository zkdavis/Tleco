import math
import numpy as np
import pytest
import tleco


NUMG = 50
G = list(np.logspace(1, 6, NUMG))          # Lorentz factor grid 10 .. 1e6
N0 = list(np.array(G) ** (-2.5))           # power-law particle distribution
GDOT = list(np.array(G) * 1e-8)            # small synchrotron-like cooling
D_ZERO = list(np.zeros(NUMG))              # no diffusion
D_SMALL = list(np.full(NUMG, 1e-5))        # small diffusion
Q_ZERO = list(np.zeros(NUMG))              # no injection
TESC = 1e200                               # effectively no escape
TLC = 1e10                                 # light-crossing time [s]
DT = 1e8                                   # one timestep [s]


class TestFpFindifDifu:

    def test_output_length(self):
        nout = tleco.fp_findif_difu(DT, G, N0, GDOT, D_ZERO, Q_ZERO, TESC, TLC)
        assert len(nout) == NUMG

    def test_no_nan_no_inf_zero_diffusion(self):
        # Before the bug fix this path crashed; now it should return finite values
        nout = tleco.fp_findif_difu(DT, G, N0, GDOT, D_ZERO, Q_ZERO, TESC, TLC)
        assert all(math.isfinite(v) for v in nout)

    def test_no_nan_no_inf_with_diffusion(self):
        nout = tleco.fp_findif_difu(DT, G, N0, GDOT, D_SMALL, Q_ZERO, TESC, TLC)
        assert all(math.isfinite(v) for v in nout)

    def test_non_negative(self):
        # Particle distributions must never go negative
        nout = tleco.fp_findif_difu(DT, G, N0, GDOT, D_SMALL, Q_ZERO, TESC, TLC)
        assert all(v >= 0.0 for v in nout)

    def test_cooling_shifts_distribution_down(self):
        # Under pure cooling (no injection, no escape) the bulk of the distribution
        # should shift to lower Lorentz factors.  The total number of particles
        # above the midpoint of the grid should decrease after one step.
        gmid = G[NUMG // 2]
        n_above_before = sum(n for g, n in zip(G, N0) if g > gmid)
        nout = tleco.fp_findif_difu(DT, G, N0, GDOT, D_ZERO, Q_ZERO, TESC, TLC)
        n_above_after = sum(n for g, n in zip(G, nout) if g > gmid)
        assert n_above_after < n_above_before

    def test_check_params_nan_raises(self):
        # Passing NaN in the grid should panic (raised as BaseException by PyO3) when check_params=True
        g_bad = list(G)
        g_bad[5] = float('nan')
        with pytest.raises(BaseException):
            tleco.fp_findif_difu(DT, g_bad, N0, GDOT, D_ZERO, Q_ZERO, TESC, TLC,
                                 check_params=True)

    def test_check_params_disabled_skips_validation(self):
        # With check_params=False the validation is skipped (no error even for NaN grid)
        g_bad = list(G)
        g_bad[5] = float('nan')
        # Should not raise — just produce (possibly garbage) output
        tleco.fp_findif_difu(DT, g_bad, N0, GDOT, D_ZERO, Q_ZERO, TESC, TLC,
                              check_params=False)


class TestEq59Park1995:

    def test_output_length(self):
        g = list(np.logspace(0, 3, 30))
        result = tleco.eq_59_park1995(1e-3, g)
        assert len(result) == 30

    def test_no_nan(self):
        g = list(np.logspace(0, 3, 30))
        result = tleco.eq_59_park1995(1e-3, g)
        assert all(math.isfinite(v) for v in result)

    def test_non_negative(self):
        g = list(np.logspace(0, 3, 30))
        result = tleco.eq_59_park1995(1e-3, g)
        assert all(v >= 0.0 for v in result)
