import math
import numpy as np
import pytest
import tleco


# Shared grids used across tests
NUMG = 50
NUMF = 60
G = list(np.logspace(1, 7, NUMG))           # Lorentz factor grid 10 .. 1e7
N = list(np.array(G) ** (-2.5))             # power-law electron distribution
B = 0.1                                     # magnetic field [G]
FREQS = list(np.logspace(8, 28, NUMF))      # frequency grid [Hz]
R = 5e14                                    # blob radius [cm]


class TestSynEmissivityFull:

    def test_output_shapes(self):
        jmbs, ambs = tleco.syn_emissivity_full(FREQS, G, N, B, True)
        assert len(jmbs) == NUMF
        assert len(ambs) == NUMF

    def test_no_nan_or_inf(self):
        jmbs, ambs = tleco.syn_emissivity_full(FREQS, G, N, B, True)
        assert all(math.isfinite(v) for v in jmbs)
        assert all(math.isfinite(v) for v in ambs)

    def test_non_negative(self):
        jmbs, ambs = tleco.syn_emissivity_full(FREQS, G, N, B, True)
        assert all(v >= 0.0 for v in jmbs)
        assert all(v >= 0.0 for v in ambs)

    def test_absorption_zero_when_not_requested(self):
        _, ambs = tleco.syn_emissivity_full(FREQS, G, N, B, False)
        assert all(v == 0.0 for v in ambs)

    def test_stronger_b_increases_emissivity(self):
        # Doubling B shifts the synchrotron peak and raises emissivity
        jmbs_lo, _ = tleco.syn_emissivity_full(FREQS, G, N, 0.1, False)
        jmbs_hi, _ = tleco.syn_emissivity_full(FREQS, G, N, 0.2, False)
        assert sum(jmbs_hi) > sum(jmbs_lo)

    def test_gmin_equals_one_no_crash(self):
        # Before the bug fix, gamma=1 in arma_trapzd caused division by zero
        g_with_rest = list(np.logspace(0, 7, NUMG))   # starts at gamma=1
        n_rest = list(np.array(g_with_rest) ** (-2.5))
        jmbs, ambs = tleco.syn_emissivity_full(FREQS, g_with_rest, n_rest, B, True)
        assert all(math.isfinite(v) for v in ambs)


class TestRadTransBlob:

    def test_output_length(self):
        jnu = [1e-20] * NUMF
        anu = [0.0] * NUMF
        inu = tleco.rad_trans_blob(R, jnu, anu)
        assert len(inu) == NUMF

    def test_optically_thin_limit(self):
        # tau -> 0 => opt_depth_blob -> 1 => I = R*j/3
        jnu = [1e-30] * NUMF           # very faint source
        anu = [1e-100] * NUMF          # negligible absorption
        inu = tleco.rad_trans_blob(R, jnu, anu)
        for j, i in zip(jnu, inu):
            assert i == pytest.approx(R * j / 3.0, rel=1e-6)

    def test_non_negative(self):
        jmbs, ambs = tleco.syn_emissivity_full(FREQS, G, N, B, True)
        inu = tleco.rad_trans_blob(R, jmbs, ambs)
        assert all(v >= 0.0 for v in inu)

    def test_no_nan(self):
        jmbs, ambs = tleco.syn_emissivity_full(FREQS, G, N, B, True)
        inu = tleco.rad_trans_blob(R, jmbs, ambs)
        assert all(math.isfinite(v) for v in inu)

    def test_larger_blob_increases_intensity(self):
        jmbs, ambs = tleco.syn_emissivity_full(FREQS, G, N, B, True)
        inu_small = tleco.rad_trans_blob(R, jmbs, ambs)
        inu_large = tleco.rad_trans_blob(10.0 * R, jmbs, ambs)
        assert sum(inu_large) > sum(inu_small)


class TestIcIsoPowlawFull:

    def _synchrotron_seed(self):
        jmbs, ambs = tleco.syn_emissivity_full(FREQS, G, N, B, True)
        inu = tleco.rad_trans_blob(R, jmbs, ambs)
        return inu

    def test_output_length(self):
        inu = self._synchrotron_seed()
        jic = tleco.ic_iso_powlaw_full(FREQS, inu, G, N)
        assert len(jic) == NUMF

    def test_no_nan_no_inf(self):
        inu = self._synchrotron_seed()
        jic = tleco.ic_iso_powlaw_full(FREQS, inu, G, N)
        assert all(math.isfinite(v) for v in jic)

    def test_non_negative(self):
        inu = self._synchrotron_seed()
        jic = tleco.ic_iso_powlaw_full(FREQS, inu, G, N)
        assert all(v >= 0.0 for v in jic)

    def test_zero_seed_gives_zero_ic(self):
        inu_zero = [0.0] * NUMF
        jic = tleco.ic_iso_powlaw_full(FREQS, inu_zero, G, N)
        assert all(v == 0.0 for v in jic)


class TestIcIsoMonochromeFull:

    def test_output_length(self):
        uext = 1e-3
        nuext = 1e14
        jic = tleco.ic_iso_monochrome_full(FREQS, uext, nuext, N, G)
        assert len(jic) == NUMF

    def test_no_nan_no_inf(self):
        uext = 1e-3
        nuext = 1e14
        jic = tleco.ic_iso_monochrome_full(FREQS, uext, nuext, N, G)
        assert all(math.isfinite(v) for v in jic)

    def test_non_negative(self):
        uext = 1e-3
        nuext = 1e14
        jic = tleco.ic_iso_monochrome_full(FREQS, uext, nuext, N, G)
        assert all(v >= 0.0 for v in jic)

    def test_higher_energy_density_gives_more_emission(self):
        nuext = 1e14
        jic_lo = tleco.ic_iso_monochrome_full(FREQS, 1e-4, nuext, N, G)
        jic_hi = tleco.ic_iso_monochrome_full(FREQS, 1e-3, nuext, N, G)
        assert sum(jic_hi) > sum(jic_lo)


class TestRadCoolMono:

    def test_output_length(self):
        dotg = tleco.rad_cool_mono(G, 1e14, 1e-3, False)
        assert len(dotg) == NUMG

    def test_no_nan(self):
        dotg = tleco.rad_cool_mono(G, 1e14, 1e-3, False)
        assert all(math.isfinite(v) for v in dotg)

    def test_non_negative(self):
        dotg = tleco.rad_cool_mono(G, 1e14, 1e-3, False)
        assert all(v >= 0.0 for v in dotg)

    def test_thomson_scales_as_p_squared(self):
        # In the Thomson limit (xi << 1) the code computes dotg ∝ p² = γ²-1.
        g_test = [100.0, 200.0, 400.0]
        dotg = tleco.rad_cool_mono(g_test, 1e6, 1e-3, False)
        expected_01 = (200.0**2 - 1.0) / (100.0**2 - 1.0)
        expected_12 = (400.0**2 - 1.0) / (200.0**2 - 1.0)
        assert dotg[1] / dotg[0] == pytest.approx(expected_01, rel=1e-4)
        assert dotg[2] / dotg[1] == pytest.approx(expected_12, rel=1e-4)

    def test_kn_reduces_cooling_at_high_gamma(self):
        # At high gamma*nu the KN cross-section suppresses cooling relative to Thomson
        g_high = list(np.logspace(5, 8, 20))
        dotg_thomson = tleco.rad_cool_mono(g_high, 1e14, 1e-3, False)
        dotg_kn = tleco.rad_cool_mono(g_high, 1e14, 1e-3, True)
        # KN cooling must be <= Thomson cooling at every gamma
        for t, k in zip(dotg_thomson, dotg_kn):
            assert k <= t * (1.0 + 1e-9)


class TestRadCoolPwl:

    def test_output_length(self):
        uu = list(np.array(FREQS) ** (-1.5) * 1e-20)
        dotg = tleco.rad_cool_pwl(G, FREQS, uu, False)
        assert len(dotg) == NUMG

    def test_no_nan(self):
        uu = list(np.array(FREQS) ** (-1.5) * 1e-20)
        dotg = tleco.rad_cool_pwl(G, FREQS, uu, False)
        assert all(math.isfinite(v) for v in dotg)

    def test_non_negative(self):
        uu = list(np.array(FREQS) ** (-1.5) * 1e-20)
        dotg = tleco.rad_cool_pwl(G, FREQS, uu, False)
        assert all(v >= 0.0 for v in dotg)


class TestRadTransSlab:

    def test_output_length(self):
        jnu = [1e-20] * NUMF
        anu = [0.0] * NUMF
        inu = tleco.rad_trans_slab(R, jnu, anu)
        assert len(inu) == NUMF

    def test_optically_thin_limit(self):
        # tau -> 0 => opt_depth_slab -> 1 => I = R*j
        jnu = [1e-30] * NUMF
        anu = [1e-100] * NUMF
        inu = tleco.rad_trans_slab(R, jnu, anu)
        for j, i in zip(jnu, inu):
            assert i == pytest.approx(R * j, rel=1e-6)

    def test_non_negative(self):
        jmbs, ambs = tleco.syn_emissivity_full(FREQS, G, N, B, True)
        inu = tleco.rad_trans_slab(R, jmbs, ambs)
        assert all(v >= 0.0 for v in inu)

    def test_no_nan(self):
        jmbs, ambs = tleco.syn_emissivity_full(FREQS, G, N, B, True)
        inu = tleco.rad_trans_slab(R, jmbs, ambs)
        assert all(math.isfinite(v) for v in inu)

    def test_larger_slab_increases_intensity(self):
        jmbs, ambs = tleco.syn_emissivity_full(FREQS, G, N, B, True)
        inu_small = tleco.rad_trans_slab(R, jmbs, ambs)
        inu_large = tleco.rad_trans_slab(10.0 * R, jmbs, ambs)
        assert sum(inu_large) > sum(inu_small)
