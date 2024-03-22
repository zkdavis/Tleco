import numpy as np
def steady_state_solution(n0,g):
    def analytical_form(g):
        C0 = 3.48e-11
        tesc = 1e0 / (C0 * ((10e0) ** (4.5e0)))
        return (g ** 2) * np.exp(-2 * C0 * tesc * (g - 1))
    x = n0/(np.trapz(analytical_form(g),g))
    return x*analytical_form(g)