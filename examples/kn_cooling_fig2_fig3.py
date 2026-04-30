#!/usr/bin/env python3
"""
Replicate Figures 2 and 3 from Rueda-Becerril (2021), arXiv:2011.13797,
combined into a single two-panel plot.

Left panel  (Fig. 2) — blackbody (CMB) radiation field
Right panel (Fig. 3) — monochromatic radiation field

Three curves per panel:
  blue   — Thomson limit (no KN correction)
  orange — numerical KN  (this code)
  green  — analytic approximation
              Fig. 2: deep-KN asymptotic for blackbody (Dermer & Menon 2009, Eq. 6.54)
              Fig. 3: (1 + xi)^{-3/2} factor (Moderski et al. 2005)
"""

import numpy as np
import matplotlib.pyplot as plt
import tleco

# ── CGS constants ─────────────────────────────────────────────────────────────
SIGMAT   = 6.6524587321e-25   # cm²
CLIGHT   = 2.99792458e10      # cm/s
HPLANCK  = 6.62607015e-27     # erg s
KBOLTZ   = 1.380649e-16       # erg/K
ENERGY_E = 8.187105776e-7     # m_e c²  [erg]
H_MEC2   = 6.53745089e-21     # h / (m_e c²)  [s]

URAD_CONST = 4.0 * SIGMAT * CLIGHT / (3.0 * ENERGY_E)  # [cm² / (erg s)]


# ── Figure 3: monochromatic field ─────────────────────────────────────────────

def figure3():
    """
    Cooling for a monochromatic photon field at frequency nu0.

    nu0 is set so that the KN threshold falls at gamma ~ 1e4:
        gamma_KN = 1 / (4 * H_MEC2 * nu0)  =>  nu0 = 1 / (4 * H_MEC2 * 1e4)
    """
    gamma_kn = 1e4
    nu0 = 1.0 / (4.0 * H_MEC2 * gamma_kn)   # Hz  (~3.8e15 Hz, EUV)
    u0  = 1e-3                                # erg/cm³ (qualitative)

    g_arr = np.logspace(3, 9, 120)
    g     = list(g_arr)
    p2    = g_arr**2 - 1.0                   # γ² − 1
    xi    = 4.0 * g_arr * nu0 * H_MEC2       # ξ = 4γhν₀/(m_e c²)

    dotg_th  = np.array(tleco.rad_cool_mono(g, nu0, u0, False))
    dotg_kn  = np.array(tleco.rad_cool_mono(g, nu0, u0, True))

    # Moderski et al. (2005): dotg_Th * (1 + ξ)^{-3/2}
    # Valid from Thomson through the trans-KN transition.
    dotg_mod = URAD_CONST * u0 * p2 * (1.0 + xi)**(-1.5)

    return g_arr, dotg_th, dotg_kn, dotg_mod


# ── Figure 2: blackbody (CMB) field ───────────────────────────────────────────

def planck_u_nu(nu, T):
    """Spectral energy density u_ν of a blackbody [erg cm⁻³ Hz⁻¹]."""
    x = np.clip(HPLANCK * nu / (KBOLTZ * T), 0.0, 700.0)
    return (8.0 * np.pi * HPLANCK * nu**3 / CLIGHT**3) / np.expm1(x)


def dermer_menon_6_54(g_arr, T):
    """
    Dermer & Menon (2009) Eq. 6.54: deep-KN asymptotic cooling for a blackbody.

        -dotgamma = (π³/2) * (c σ_T / λ_C³) * Θ² * (ln(4γΘ) - 5/6 - C_E - C_l)

    Θ  = kT / (m_e c²)       dimensionless temperature
    λ_C = h / (m_e c)         Compton wavelength
    C_E = 0.5772...            Euler–Mascheroni constant
    C_l = 0.5700...            (6/π²) Σ_{k≥2} ln(k)/k²

    The bracket (ln(4γΘ) - 5/6 - C_E - C_l) is positive only for
    4γΘ > exp(5/6 + C_E + C_l) ≈ exp(1.98) ≈ 7.2, i.e. γ ≳ 1.8/(4Θ).
    Below this γ the formula is not valid and NaN is returned, which
    matplotlib simply omits — this is the "gap" described in the paper.
    """
    LAMBDA_C = HPLANCK * CLIGHT / ENERGY_E        # h / (m_e c)  [cm]
    C_E      = 0.5772156649015329                  # Euler–Mascheroni
    C_l      = 0.5700                              # (6/π²) Σ ln(k)/k²
    THETA    = KBOLTZ * T / ENERGY_E               # dimensionless temperature

    prefactor = (np.pi**3 / 2.0) * (CLIGHT * SIGMAT / LAMBDA_C**3)
    bracket   = np.log(4.0 * g_arr * THETA) - 5.0/6.0 - C_E - C_l

    # NaN outside the regime of validity so matplotlib leaves a visible gap
    return np.where(bracket > 0.0, prefactor * THETA**2 * bracket, np.nan)


def figure2():
    """
    Cooling for a CMB blackbody field (T = 2.73 K).

    The CMB peak is at ~160 GHz; the KN threshold for the CMB is
    γ_KN ~ 1/(4Θ) ~ 5e8. Eq. 6.54 first becomes positive at γ ~ 4e9,
    creating a visible gap between the Thomson and deep-KN curves — exactly
    the feature the numerical approach (rad_cool_pwl with KN) fills smoothly.
    """
    T_cmb = 2.73   # K

    # Frequency grid spanning the CMB peak (~160 GHz) with ample margin
    numf      = 400
    freqs_arr = np.logspace(8, 14, numf)
    u_nu_arr  = planck_u_nu(freqs_arr, T_cmb)

    freqs = list(freqs_arr)
    u_nu  = list(u_nu_arr)

    g_arr = np.logspace(7, 13, 100)
    g     = list(g_arr)

    dotg_th = np.array(tleco.rad_cool_pwl(g, freqs, u_nu, False))
    dotg_kn = np.array(tleco.rad_cool_pwl(g, freqs, u_nu, True))
    dotg_dm = dermer_menon_6_54(g_arr, T_cmb)

    return g_arr, dotg_th, dotg_kn, dotg_dm


# ── Plotting ──────────────────────────────────────────────────────────────────

def main():
    g2, th2, kn2, dm2 = figure2()
    g3, th3, kn3, mod3 = figure3()

    fig, (ax2, ax3) = plt.subplots(1, 2, figsize=(11, 4.5))

    # ── Left panel: Fig. 2 ────────────────────────────────────────────────────
    ax2.loglog(g2, th2,  color='C0', lw=2,   label='Thomson')
    ax2.loglog(g2, kn2,  color='C1', lw=2,   label='Numerical KN (this work)')
    ax2.loglog(g2, dm2,  color='C2', lw=1.5, ls='--',
               label=r'D&M (2009) Eq. 6.54 — deep-KN asymptote ($\xi\gg1$)')

    ax2.set_xlabel(r'$\gamma$', fontsize=12)
    ax2.set_ylabel(r'$|\dot{\gamma}|\;[\mathrm{s}^{-1}]$', fontsize=12)
    ax2.set_title(
        r'Blackbody field (CMB, $T=2.73\,\mathrm{K}$)'
        '\n'
        r'D&M Eq. 6.54 underestimates by $\sim\!40\%$'
        r' (deep-KN kernel invalid for $\xi<6.2$)',
        fontsize=10)
    ax2.set_xlim(1e7, 1e13)
    ax2.set_ylim(1e-6, 1e-1)
    ax2.legend(fontsize=8)
    ax2.grid(True, which='both', alpha=0.25)

    # ── Right panel: Fig. 3 ───────────────────────────────────────────────────
    ax3.loglog(g3, th3,  color='C0', lw=2,   label='Thomson')
    ax3.loglog(g3, kn3,  color='C1', lw=2,   label='Numerical KN (this work)')
    ax3.loglog(g3, mod3, color='C2', lw=1.5, ls='--',
               label='Moderski et al. (2005) approx.')

    ax3.set_xlabel(r'$\gamma$', fontsize=12)
    ax3.set_ylabel(r'$|\dot{\gamma}|\;[\mathrm{s}^{-1}]$', fontsize=12)
    ax3.set_title(r'Monochromatic field ($\nu_0 \approx 3.8\times10^{15}$ Hz, '
                  r'$u_0=10^{-3}$ erg cm$^{-3}$)', fontsize=11)
    ax3.set_xlim(1e3, 1e9)
    ax3.set_ylim(1e-5, 1e1)
    ax3.legend(fontsize=9)
    ax3.grid(True, which='both', alpha=0.25)

    plt.tight_layout()

    out_png = 'examples/Figs/kn_cooling_fig2_fig3.png'
    plt.savefig(out_png, dpi=150, bbox_inches='tight')
    print(f'Saved {out_png}')
    plt.show()


if __name__ == '__main__':
    main()
