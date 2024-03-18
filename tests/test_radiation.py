from paramo.comparisons import ic_bb_test,ic_mono_test

class TestRadiation(object):
    def test_rad_trans_blob(self) -> None:
        assert 1 == 1

    def test_ic_iso_powlaw(self) -> None:
        num_f150_num_g150_error = 0.22399
        numg=150
        numf=150
        nu_bounds = [5e8, 5e22]
        err = ic_bb_test.ic_bb_get_error(num_f=numf, num_g=numg, nu_bounds=nu_bounds)
        assert err <= num_f150_num_g150_error


    def test_ic_iso_mono(self) -> None:
        num_f150_num_g150_error = 0.00166
        numg=150
        numf=150
        nu_bounds = [2e13,2.5e19]
        err = ic_mono_test.ic_mono_get_error(num_f=numf, num_g=numg, nu_bounds=nu_bounds)
        assert err <= num_f150_num_g150_error

