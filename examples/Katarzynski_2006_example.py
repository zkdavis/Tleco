import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors
from tleco import misc_func
from tleco import conversion_funcs as cf
import tleco as para
from tleco import constants as cons
from tqdm import tqdm
from ebltable.tau_from_model import OptDepth
import sys

# Constants
with_abs, cool_withKN = True, True
num_t0, numg, numf = 150, 180, 150
fmin, fmax, gmin, gmax = 1e7, 1e30, 1e0, 1e10
R, B = 3.2e15, 0.05  # Blob size and magnetic field
R2, B2 = 1e15, 0.11  # Blob size and magnetic field
uB = (B ** 2) / (np.pi * 8)  # Magnetic field energy density
uB2 = (B2 ** 2) / (np.pi * 8)  # Magnetic field energy density
z = 0.03364  # mrk 501 redshift
C0 = 3.48e-11  # cooling constant density constant
t_inj = R / cons.cLight  # Injection time
tlc = t_inj  # Light crossing time
t_acc = 1 / (C0 * 1e4)  # Acceleration time scale
t_acc7 = tlc  # Acceleration time scale
t_esc = t_acc  # Escape time
t_acc8 = R2 / cons.cLight
t_inj8 = t_acc8
t_esc8 = t_acc8
t_times1 = np.array([0.01, 0.1, 0.5, 1, 2, 3, 5, 10, 15, 40]) * t_acc
t_times2 = np.array([0.1, 0.8, 2.5, 4.5, 5, 7, 10, 20, 40, 60, 80]) * t_acc7
t_times3 = np.array([0.1, 0.5, 1, 2, 2.5, 3, 3.5, 4, 4.7, 5.6, 7, 10]) * t_acc8
t = np.logspace(0, np.log10(t_times1[-1]), num_t0)
t = np.unique(np.concatenate((t, t_times1)))
numt = len(t)
t2 = np.logspace(0, np.log10(t_times2[-1]), num_t0)
t2 = np.unique(np.concatenate((t2, t_times2)))
numt2 = len(t2)
t3 = np.logspace(0, np.log10(t_times3[-1]), num_t0)
t3 = np.unique(np.concatenate((t3, t_times3)))
numt3 = len(t3)
max_numt = np.max([numt, numt2, numt3])


def run_katarzynski():
    # Arrays
    f = np.logspace(np.log10(fmin), np.log10(fmax), numf)
    g = np.logspace(np.log10(gmin), np.log10(gmax), numg)
    D_0 = 0.5 * np.power(g, 2) / t_acc
    D_07 = 0.5 * np.power(g, 2) / t_acc7
    D_08 = 0.5 * np.power(g, 2) / t_acc8
    D1 = 2 * D_0
    D2 = D1
    D3 = D1
    D4 = D1
    D5 = D1
    D6 = D1
    D7 = 2 * D_07
    D8 = 2 * D_08
    gdot1 = np.zeros([numt, numg])
    gdot2 = np.zeros([numt, numg])
    gdot3 = np.zeros([numt, numg])
    gdot4 = np.zeros([numt, numg])
    gdot5 = np.zeros([numt, numg])
    gdot6 = np.zeros([numt, numg])
    gdot7 = np.zeros([numt2, numg])
    gdot8 = np.zeros([numt3, numg])
    gdot1[0, :] = C0 * np.power(g, 2) - 4 * D_0 / g
    gdot2[0, :] = C0 * np.power(g, 2) - 4 * D_0 / g
    gdot3[0, :] = C0 * np.power(g, 2) - 4 * D_0 / g
    gdot4[0, :] = C0 * np.power(g, 2) - 4 * D_0 / g
    gdot5[0, :] = C0 * np.power(g, 2) - 4 * D_0 / g
    gdot6[0, :] = C0 * np.power(g, 2) - 4 * D_0 / g
    gdot7[0, :] = (4 * cons.sigmaT * uB * np.power(g, 2) / (3 * cons.me * cons.cLight)) - (4 * D_07 / g)
    gdot8[0, :] = (4 * cons.sigmaT * uB2 * np.power(g, 2) / (3 * cons.me * cons.cLight)) - (4 * D_08 / g)

    # Initial particle distribution
    p, gcut1, g2cut1 = 0, 1e0, 2e0
    p, gcut2, g2cut2 = 0, 1e6, 2e6
    p, gcut3, g2cut3 = 0, 1e1, 1e6
    p, gcut4, g2cut4 = 0, 1e0, 2e0
    p, gcut5, g2cut5 = 0, 1e6, 2e6
    p, gcut6, g2cut6 = 0, 1e1, 1e6
    p, gcut7, g2cut7 = 0, 1e0, 2e0
    p, gcut8, g2cut8 = 0, 1e0, 2e0
    n1 = np.zeros([numt, numg])
    n2 = np.zeros([numt, numg])
    n3 = np.zeros([numt, numg])
    n4 = np.zeros([numt, numg])
    n5 = np.zeros([numt, numg])
    n6 = np.zeros([numt, numg])
    n7 = np.zeros([numt2, numg])
    n8 = np.zeros([numt3, numg])
    n1[0, :] = misc_func.power_law(1, g, p, gcut1, g2cut1)
    n2[0, :] = misc_func.power_law(1, g, p, gcut2, g2cut2)
    n3[0, :] = misc_func.power_law(1, g, p, gcut3, g2cut3)
    n4[0, :] = misc_func.power_law(1, g, p, gcut4, g2cut4)
    n5[0, :] = misc_func.power_law(1, g, p, gcut5, g2cut5)
    n6[0, :] = misc_func.power_law(1, g, p, gcut6, g2cut6)
    n7[0, :] = misc_func.power_law(7, g, p, gcut7, g2cut7)
    n8[0, :] = misc_func.power_law(5e-2, g, p, gcut8, g2cut8)

    Q_inj8 = misc_func.power_law(5e-2, g, p, gcut8, g2cut8)

    j_s7 = np.zeros([numt2, numf])  # synchrotron emissivity
    ambs7 = np.zeros([numt2, numf])  # absorbtion coefficient
    j_ssc7 = np.zeros([numt2, numf])  # ssc emissivity
    I_s7 = np.zeros([numt2, numf])  # synchrotron Intensity
    I_ssc7 = np.zeros([numt2, numf])  # ssc Intensity

    j_s8 = np.zeros([numt3, numf])  # synchrotron emissivity
    ambs8 = np.zeros([numt3, numf])  # absorbtion coefficient
    j_ssc8 = np.zeros([numt3, numf])  # ssc emissivity
    I_s8 = np.zeros([numt3, numf])  # synchrotron Intensity
    I_ssc8 = np.zeros([numt3, numf])  # ssc Intensity

    progress_bar = tqdm(total=max_numt - 1, desc='Progress', file=sys.stdout, colour="#19C819")

    # Time loop
    for i in range(1, max_numt):
        if (i < numt):
            dt = t[i] - t[i - 1]
            n1[i, :] = para.fp_findif_difu(dt, g, n1[i - 1, :], gdot1[i - 1, :], D1, np.zeros_like(D1), 1e200, tlc,
                                           False)
            n2[i, :] = para.fp_findif_difu(dt, g, n2[i - 1, :], gdot2[i - 1, :], D2, np.zeros_like(D2), 1e200, tlc,
                                           False)
            n3[i, :] = para.fp_findif_difu(dt, g, n3[i - 1, :], gdot3[i - 1, :], D3, np.zeros_like(D3), 1e200, tlc,
                                           False)
            n4[i, :] = para.fp_findif_difu(dt, g, n4[i - 1, :], gdot4[i - 1, :], D4, (2.5e4) * n4[0, :] / t_esc, t_esc,
                                           tlc, False)
            n5[i, :] = para.fp_findif_difu(dt, g, n5[i - 1, :], gdot5[i - 1, :], D5, (2.5e4) * n5[0, :] / t_esc, t_esc,
                                           tlc, False)
            n6[i, :] = para.fp_findif_difu(dt, g, n6[i - 1, :], gdot6[i - 1, :], D6, (2.5e4) * n6[0, :] / t_esc, t_esc,
                                           tlc, False)

            gdot1[i, :] = gdot1[0, :]
            gdot2[i, :] = gdot2[0, :]
            gdot3[i, :] = gdot3[0, :]
            gdot4[i, :] = gdot4[0, :]
            gdot5[i, :] = gdot5[0, :]
            gdot6[i, :] = gdot6[0, :]

        if (i < numt2):
            dt2 = t2[i] - t2[i - 1]
            n7[i, :] = para.fp_findif_difu(dt2, g, n7[i - 1, :], gdot7[i - 1, :], D7, np.zeros_like(D7), 1e200, tlc,False)
            j_s7[i, :], ambs7[i, :] = para.syn_emissivity_full(f, g, n7[i, :], B, with_abs)  # ,sync and absorb
            I_s7[i, :] = para.rad_trans_slab(R, j_s7[i, :], ambs7[i, :])
            j_ssc7[i, :] = para.ic_iso_powlaw_full(f, I_s7[i, :], g, n7[i, :])
            I_ssc7[i, :] = para.rad_trans_slab(R, j_ssc7[i, :], ambs7[i, :])

            # katarzynski has a custom cooling approximation
            # dotgKN7 = np.zeros_like(g)
            # for j,gi in enumerate(g):
            #     vmax_kn = 3 * cons.energy_e / (4 * cons.hPlanck * gi)
            #     vmax = min([vmax_kn,fmax])
            #     vmax_i = np.argmin(np.abs(f-vmax))
            #     urad = (np.pi * 4 / cons.cLight)*np.trapz(I_s7[i,:vmax_i],f[:vmax_i])
            #     dotgKN7[j] = (4 * cons.sigmaT* urad * np.power(gi,2)/ (3 * cons.me * cons.cLight))

            dotgKN7 = para.rad_cool_pwl(g, f, 4 * np.pi * I_s7[i, :] / cons.cLight, cool_withKN)

            gdot7[i, :] = gdot7[0, :] + dotgKN7

        if (i < numt3):
            dt3 = t3[i] - t3[i - 1]
            n8[i, :] = para.fp_findif_difu(dt3, g, n8[i - 1, :], gdot8[i - 1, :], D8, Q_inj8, t_esc8, t_acc8, False)

            j_s8[i, :], ambs8[i, :] = para.syn_emissivity_full(f, g, n8[i, :], B2, with_abs)  # ,sync and absorb
            I_s8[i, :] = para.rad_trans_slab(R2, j_s8[i, :], ambs8[i, :])
            j_ssc8[i, :] = para.ic_iso_powlaw_full(f, I_s8[i, :], g, n8[i, :])
            I_ssc8[i, :] = para.rad_trans_slab(R2, j_ssc8[i, :], ambs8[i, :])

            # katarzynski has a custom cooling approximation
            # dotgKN8 = np.zeros_like(g)
            # for j, gi in enumerate(g):
            #     vmax_kn = 3 * cons.energy_e / (4 * cons.hPlanck * gi)
            #     vmax = min([vmax_kn, fmax])
            #     vmax_i = np.argmin(np.abs(f - vmax))
            #     urad = (np.pi * 4 / cons.cLight) * np.trapz(I_s8[i, :vmax_i], f[:vmax_i])
            #     dotgKN8[j] = (4 * cons.sigmaT * urad * np.power(gi, 2) / (3 * cons.me * cons.cLight))

            dotgKN8 = para.rad_cool_pwl(g, f, 4 * np.pi * I_s8[i, :] / cons.cLight, cool_withKN)

            gdot8[i, :] = gdot8[0, :] + dotgKN8

        progress_bar.update(1)

    return [n1, n2, n3, n4, n5, n6, n7, n8], g, [I_s7, I_s8], [I_ssc7, I_ssc8], f


def get_data(file_path):
    x, y = [], []
    with open(file_path, 'r') as f:
        for line in f.readlines():
            parts = line.split(", ")
            x.append(float(parts[0]))
            y.append(float(parts[1]))
    return np.array(x), np.array(y)


def n_plot(n, g, compare_fig=None):
    fig, ax = plt.subplots(figsize=(14, 12))
    ax.set_yscale('log')
    ax.set_xscale('log')
    if compare_fig:
        if (compare_fig == 'fig1_a'):
            t_temp = t
            t_times_temp = np.array([0.01, 0.1, 0.5, 1, 2, 3, 5, 10, 15, 40]) * t_acc
            t_acc_temp = t_acc
            ax.set_ylim(1e-12, 1e1)
            ax.set_xlim(1e0, 1e7)
        elif (compare_fig == 'fig1_b'):
            t_temp = t
            t_times_temp = np.array([0.01, 0.1, 0.5, 1, 2, 3, 5, 10, 15, 40]) * t_acc
            t_acc_temp = t_acc
            ax.set_ylim(1e-10, 1e3)
            ax.set_xlim(1e0, 1e7)
        elif (compare_fig == 'fig1_c'):
            t_temp = t
            t_times_temp = np.array([0.01, 0.1, 0.5, 1, 2, 3, 5, 10, 15, 40]) * t_acc
            t_acc_temp = t_acc
            ax.set_ylim(1e-9, 1e4)
            ax.set_xlim(1e0, 1e7)
        elif (compare_fig == 'fig2_a'):
            t_temp = t
            t_times_temp = np.array([0.01, 0.1, 0.5, 1, 2, 3, 5, 10]) * t_acc
            t_acc_temp = t_acc
            ax.set_ylim(1e-10, 1e5)
            ax.set_xlim(1e0, 1e7)
        elif (compare_fig == 'fig2_b'):
            t_temp = t
            t_times_temp = np.array([0.01, 0.1, 0.5, 1, 2, 3, 5, 10]) * t_acc
            t_acc_temp = t_acc
            ax.set_ylim(1e-8, 1e7)
            ax.set_xlim(1e0, 1e7)
        elif (compare_fig == 'fig2_c'):
            t_temp = t
            t_times_temp = np.array([0.01, 0.1, 0.5, 1, 2, 3, 5, 10]) * t_acc
            t_acc_temp = t_acc
            ax.set_ylim(1e-8, 1e7)
            ax.set_xlim(1e0, 1e7)
        elif (compare_fig == 'fig3_a'):
            t_temp = t2
            t_times_temp = np.array([0.1, 0.8, 2.5, 4.5, 5, 7, 10, 20, 40, 60, 80]) * t_acc7
            t_acc_temp = t_acc7
            ax.set_ylim(1e-12, 1e2)
            ax.set_xlim(1e0, 1e8)
        elif (compare_fig == 'fig4_a'):
            t_temp = t3
            t_times_temp = np.array([0.1, 0.5, 1, 2, 3, 4.7, 5.6, 7]) * t_acc8
            t_acc_temp = t_acc8
            ax.set_ylim(1e-12, 1e4)
            ax.set_xlim(1e0, 1e8)
        else:
            t_temp = t
            t_times_temp = np.array([0.01, 0.1, 0.5, 1, 2, 3, 5, 10]) * t_acc
            t_acc_temp = t_acc
            ax.set_ylim(1e-9, 1e4)
            ax.set_xlim(1e0, 1e7)

    scale_mult = 3
    colors = cm.rainbow(np.linspace(0, 1, len(t_times_temp)))
    bounds = np.append(0, np.sort(t_times_temp))
    cmap = matplotlib.colors.ListedColormap(colors)
    norm = matplotlib.colors.BoundaryNorm(boundaries=bounds, ncolors=len(t_times_temp))
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    plot_times = np.abs(np.subtract.outer(t_temp, t_times_temp)).argmin(0)
    j = 0
    for i in np.append(0, plot_times):
        color = cmap(j / (len(plot_times) + 1))
        j += 1
        ax.plot(g, n[i, :], color=color, linewidth=2 * scale_mult)

    cbar = plt.colorbar(sm, ax=ax, ticks=t_times_temp, pad=0)
    cbar.set_label('$t / t_{acc}$', rotation=-90, labelpad=40, fontsize=15 * scale_mult)
    cbar.set_ticks(bounds[:-1] + (bounds[1:] - bounds[:-1]) / 2)
    cbar.ax.set_yticklabels([f'{time / t_acc_temp:.1f}' for time in t_times_temp],  fontsize=10 * scale_mult)
    cbar.ax.tick_params(size=13, width=2, labelsize=10 * scale_mult)

    ax.tick_params(axis='both', which='major', size=12 * scale_mult, labelsize=12 * scale_mult)
    ax.tick_params(axis='both', which='minor', size=0 * scale_mult)

    ax.set_xlabel(r"$\gamma$", fontsize=15 * scale_mult)
    ax.set_ylabel(r"$n(\gamma)$ [cm$^{-3}$]", fontsize=15 * scale_mult)
    ax.set_title('')

    plt.tight_layout()
    plt.savefig(f"Figs/{compare_fig}.png")
    plt.savefig(f"Figs/{compare_fig}.pdf")
    # plt.show()


def nuFnu_plots(Im, nu, compare_fig):
    fig, ax = plt.subplots(figsize=(14, 12))
    ax.set_yscale('log')
    ax.set_xscale('log')
    if compare_fig == 'fig3_b':
        t_temp = t2
        t_times_temp = np.array([2.5, 3.5, 4.0, 4.5, 5, 7]) * t_acc7
        t_acc_temp = t_acc7
        d_temp = 4.32e26
        dop_temp = 21
        R_temp = R
        ax.set_xlim(1e10, 1e22)
        ax.set_ylim(1e-16, 1e-5)
    elif compare_fig == 'fig3_c':
        t_temp = t2
        t_times_temp = np.array([2.5, 3.5, 4.0, 4.5, 5, 7]) * t_acc7
        t_acc_temp = t_acc7
        R_temp = R
        d_temp = 4.32e26
        dop_temp = 21
        ax.set_xlim(1e18, 1e29)
        ax.set_ylim(1e-16, 1e-6)
    elif compare_fig == 'fig4_b':
        t_temp = t3
        t_times_temp = np.array([2, 2.5, 3, 3.5, 4, 4.7, 5.6, 7, 10]) * t_acc8
        t_acc_temp = t_acc8
        d_temp = 4.32e26
        R_temp = R2
        dop_temp = 33
        ax.set_xlim(1e10, 1e22)
        ax.set_ylim(1e-16, 1e-7)
    elif compare_fig == 'fig4_c':
        t_temp = t3
        t_times_temp = np.array([2.5, 3, 3.5, 4, 4.7, 5.6, 7, 10]) * t_acc8
        t_acc_temp = t_acc8
        d_temp = 4.32e26
        R_temp = R2
        dop_temp = 33
        ax.set_xlim(1e18, 1e29)
        ax.set_ylim(1e-16, 1e-8)

    scale_mult = 3
    colors = cm.rainbow(np.linspace(0, 1, len(t_times_temp)))
    bounds = np.append(0, np.sort(t_times_temp))
    cmap = matplotlib.colors.ListedColormap(colors)
    norm = matplotlib.colors.BoundaryNorm(boundaries=bounds, ncolors=len(t_times_temp))
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    tau_gg = OptDepth.readmodel(model='kneiske')
    ETeV = cf.erg2ev(cons.hPlanck * nu * dop_temp) / 1e12
    atten = np.exp(-1. * tau_gg.opt_depth(z, ETeV))
    plot_times = np.abs(np.subtract.outer(t_temp, t_times_temp)).argmin(0)
    j = 0
    for i in np.append(0, plot_times):
        color = cmap(j / (len(plot_times) + 1))
        j += 1
        nu_Fnu = atten * (dop_temp ** 4) * nu * (Im[i, :]) * (4 * np.pi * (R_temp ** 2)) / (4 * np.pi * (d_temp) ** 2)
        ax.plot(dop_temp * nu, nu_Fnu, color=color,  linewidth=2 * scale_mult)

    ax.set_xlabel(r"$\nu$ [Hz]", fontsize=15 * scale_mult)
    ax.set_ylabel(r"$\nu F_{\nu}$ $[\frac{erg}{s \ cm^{2}}]$", fontsize=15 * scale_mult)

    cbar = plt.colorbar(sm, ax=ax, ticks=t_times_temp, pad=0)
    cbar.set_label('$t / t_{acc}$', rotation=-90, labelpad=25, fontsize=15 * scale_mult)
    cbar.set_ticks(bounds[:-1] + (bounds[1:] - bounds[:-1]) / 2)
    cbar.ax.set_yticklabels([f'{time / t_acc_temp:.1f}' for time in t_times_temp], fontsize=10 * scale_mult)
    cbar.ax.tick_params(size=13, width=2, labelsize=10 * scale_mult)

    ax.tick_params(axis='both', which='major', size=12 * scale_mult, labelsize=12 * scale_mult)
    ax.tick_params(axis='both', which='minor', size=0 * scale_mult)

    plt.tight_layout()
    plt.savefig(f"Figs/{compare_fig}.png")
    plt.savefig(f"Figs/{compare_fig}.pdf")
    # plt.show()


if __name__ == '__main__':
    ns, g, I_syns, I_sscs, nu = run_katarzynski()
    compare_figs = ['fig1_a', 'fig1_b', 'fig1_c', 'fig2_a', 'fig2_b', 'fig2_c', 'fig3_a', 'fig4_a']
    compare_nfnu_figs = ['fig3_b', 'fig4_b', ]
    compare_nfnu_figs_ssc = ['fig3_c', 'fig4_c']
    for i, n in enumerate(ns):
        n_plot(n, g, compare_fig=compare_figs[i])
        if (i > 5):
            Is = I_syns[i - 6]
            Issc = I_sscs[i - 6]
            nuFnu_plots(Is, nu, compare_fig=compare_nfnu_figs[i - 6])
            nuFnu_plots(Issc, nu, compare_fig=compare_nfnu_figs_ssc[i - 6])
