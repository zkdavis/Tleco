"""this file create mirror for referencing the rust wrapped files in lib.rs"""

# Example retrieving constant from constants file can be improved by using macros
def get_c() -> float: ...

# functions from distribs
def  eq_59_park1995(t: float, g: [float]) -> [float]: ...

def  fp_findif_difu(dt_in: float, gamma_bins: [float], nin: [float], gdot_in: [float],
               din: [float], qin: [float], tesc_in: float, tlc: float,
               check_params: bool) -> [float]: ...

def syn_emissivity_full(freqs: [float], gamma_bins: [float], n_distrib: [float], b_field: float, with_abs: bool) -> ([float], [float]): ...

def rad_trans_blob(blob_radius: float, j_nu: [float], a_nu: [float]) -> [float]: ...

def rad_trans_slab(blob_radius: float, j_nu: [float], a_nu: [float]) -> [float]: ...

def ic_iso_powlaw_full(freqs: [float], inu: [float], g: [float], n: [float]) -> [float]: ...

def rad_cool_pwl(gg: [float], freqs: [float], uu: [float], with_kn: bool) -> [float]: ...

def rad_cool_mono(gg: [float],  nu0: float, u0: float, with_kn: bool) -> [float]: ...

def ic_iso_monochrome_full(freqs: [float], uext: float,nuext: float, n: [float], g: [float]) -> [float]: ...
