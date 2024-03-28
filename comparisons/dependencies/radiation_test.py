import numpy as np
import typing
import scipy
from scipy import integrate
from tleco import constants as C, misc_func as mf


# Compton distribution function: see Dermer 6.75
def Fc(q, Gam):
    fc = (2 * q * np.log(q)) + ((1 + (2 * q)) * (1 - q)) + ((((Gam * q) ** 2) * (1 - q)) / (2 * (1 + (Gam * q))))
    return fc


def g_min_f(eps_s, eps):
    return (eps_s / 2) * (1 + np.sqrt(1 + (1 / (eps * eps_s))))


def g_max_f(eps_s, eps, g2):
    x = eps - eps_s
    f = np.zeros_like(x)
    for i in range(len(x)):
        if (x[i] > 0):
            f[i] = eps[i] * eps_s / (eps[i] - eps_s)
        else:
            f[i] = g2
    return f


##eq 6.89 dermer
##inverse compton emission from isotropic power law electrons and an isotropic photon field
def j_ic_iso(eps_s: float, eps_array, u: typing.Callable, n_array: [float], g_array: [float]):
    fcon = (3 / 4) * C.cLight * C.sigmaT * eps_s

    def g_int(Fc_array, g_min_i, g_max_i):
        # assumes g_array is logspaced
        return np.trapz(n_array[g_min_i:g_max_i] * Fc_array[g_min_i:g_max_i] * np.log(10) / g_array[g_min_i:g_max_i],
                        np.log10(g_array[g_min_i:g_max_i]))

    g_int_array = np.zeros_like(eps_array)
    g_max_array = g_max_f(eps_s, eps_array, g_array[-1])
    g_min_array = g_min_f(eps_s, eps_array)
    for i in range(len(eps_array)):
        eps_i = eps_array[i]
        Gamma_i_array = 4 * g_array * eps_i
        q_i_array = eps_s / (g_array * Gamma_i_array * (1 - (eps_s / g_array)))
        for j in range(len(q_i_array)):
            q_i_array[j] = mf.clamp(q_i_array[j], (4 * g_array[j]) ** (-1), 1)
        Fc_i_array = Fc(q=q_i_array, Gam=Gamma_i_array)
        gmin_i = np.argmin(np.abs(g_array - g_min_array[i]))
        gmax_i = np.argmin(np.abs(g_array - g_max_array[i]))
        g_int_array[i] = g_int(Fc_i_array, gmin_i, gmax_i)
    # assumes eps_array is logspaced
    return fcon * np.trapz(u(eps_array) * g_int_array * np.log(10) / eps_array, np.log10(eps_array))


def j_ic_iso_full_dermer(eps_s_array: [float], eps_array: [float], u: typing.Callable, n_array: [float],
                         g_array: [float]):
    j_ic = np.zeros_like(eps_s_array)
    for i in range(len(eps_s_array)):
        eps_s_i = eps_s_array[i]
        j_ic[i] = j_ic_iso(eps_s_i, eps_array, u, n_array, g_array)

    return j_ic


# 5.14 Dermer 2009
def BB_energy_density(eps_array, Temp):
    nu = eps_array/C.hPlanck
    fcon = 8*np.pi*C.hPlanck*(nu**3)/(C.cLight**3)
    f = 1/(np.exp(eps_array / (C.kBoltz * Temp)) - 1)
    return fcon * f


##eq 6.72 from Dermer
def j_ic_mono_pwl_electron_dermer(eps_s: float, eps_in, u0, p, n_array: [float], g_array: [float]):
    if (eps_in > eps_s):
        return 0
    fcon = 0.5 * C.cLight * C.sigmaT * u0 * np.power(eps_s / eps_in, 2)
    ll = np.max([eps_s / 2, 0.5 * np.sqrt(eps_s / eps_in)])
    gmin_i = np.argmin(np.abs(g_array - ll))
    fg = 1 / (2 * g_array * eps_in)
    for i in range(len(fg)):
        if (fg[i] > 1):
            fg[i] = 1
    eps_hat = eps_s / (4 * eps_in * np.power(g_array, 2))
    integrand = n_array * (fg - eps_hat) * np.power(g_array, -1) * np.log(10)
    f = np.trapz(integrand[gmin_i:], np.log10(g_array[gmin_i:]))

    return C.hPlanck * fcon * f / (eps_s * C.energy_e)


def j_ic_mono_pwl_electron_dermer_full(eps_s_array: [float], eps_in, u0, p, n_array: [float], g_array: [float]):
    j_ic = np.zeros_like(eps_s_array)
    for i in range(len(eps_s_array)):
        eps_s_i = eps_s_array[i]
        j_ic[i] = j_ic_mono_pwl_electron_dermer(eps_s_i, eps_in, u0, p, n_array, g_array)

    return j_ic


def j_syn_explicit(nu, B, n, g, rtol=1.48e-10, tol=1.48e-10, divmax=10):
    ##dermer 7.44
    ## calculates the average power emitted from a distribution
    ## this method will be slow and is not recommended

    def W(k, mu, x):
        # whittaker function
        f = scipy.special.hyperu(mu - k + (1 / 2), 1 + 2 * mu, x)
        f = f * np.exp(-x / 2) * (x ** (mu + (1 / 2)))
        return f

    def R(x):
        # dermer 7.45
        f = 0.5 * np.pi * x
        f = f * (W(0, 4 / 3, x) * W(0, 1 / 3, x) - W(1 / 2, 5 / 6, x) * W(-1 / 2, 5 / 6, x))
        return f

    def f(gv, i):
        gv = 10 ** gv
        nuB = C.eCharge * B / (np.pi * 2 * C.me * C.cLight)
        x = 2 * nu[i] / (3 * nuB * (gv ** 2))
        gs = np.argsort(np.abs(g - gv))
        m = (n[gs[0]] - n[gs[1]]) / (g[gs[0]] - g[gs[1]])
        nf = n[gs[0]] + m * (g[gs[0]] - gv)
        return nf * R(x) * gv * np.log(10)

    jcon = np.sqrt(3) * (C.eCharge ** 3) / (C.me * (C.cLight ** 2) * np.pi * 4)
    j = np.zeros(len(nu))
    for i in range(len(nu)):
        i1 = integrate.romberg(f, np.log10(g[0]), np.log10(g[-1]), args=([i]), rtol=rtol, tol=tol, divmax=divmax)
        j[i] = jcon * B * i1

    return j
