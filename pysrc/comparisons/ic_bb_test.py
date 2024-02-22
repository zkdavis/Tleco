import matplotlib.pyplot as plt
import paramo as para
import numpy as np
import constants as C
from ic_dependencies import radiation_test as rt

if __name__ == "__main__":
    p = 2
    n0 =1
    sig_T = 0.0001
    Temp = sig_T*C.energy_e/C.kBoltz
    u_bb = 1e0
    num_eps_s = 100
    numg = 100
    nu_s = np.logspace(8, 25, num_eps_s)
    eps_s = nu_s * C.hPlanck / (C.me * (C.cLight ** 2))
    g_array = np.logspace(1, 3, numg)
    n_array = n0*(g_array ** -p) / np.trapz(g_array ** -p, g_array)


    def ub(epsa):
        return rt.BB_photon_density(epsa,Temp)*epsa*u_bb/np.trapz(rt.BB_photon_density(epsa,Temp)*epsa,epsa)

    def ub_dermer(epsa):
        epsa = epsa * C.energy_e
        return rt.BB_photon_density(epsa, Temp) * epsa * u_bb / np.trapz(rt.BB_photon_density(epsa, Temp) * epsa, epsa)


    j_ic = rt.j_ic_iso_full_dermer(eps_s, eps_s,ub_dermer, n_array, g_array)#/C.energy_e


    j_ic_para = np.array(para.ic_iso_powlaw_full(nu_s, C.cLight * ub(nu_s*C.hPlanck), g_array, n_array))

    fig, ax = plt.subplots()

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.plot(nu_s, j_ic_para,label='PARAMO')
    ax.plot(nu_s, j_ic)
    ax.legend()
    plt.show()