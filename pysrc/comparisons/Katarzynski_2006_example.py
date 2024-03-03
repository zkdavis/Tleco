import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors
import misc_func
import paramo as para
import constants as cons

# Constants
with_abs, cool_withKN = True, True
num_t, numg, numf = 150, 300, 300
fmin, fmax, tmax, gmin, gmax = 1e8, 1e28, 1e18, 1e0, 1e8
R, B = 3.2e15, 0.05  # Blob size and magnetic field
uB = (B ** 2) / (np.pi * 8)  # Magnetic field energy density
C0 = 3.48e-11  # Energy density constant
t_inj = R / cons.cLight  # Injection time
tlc = t_inj  # Light crossing time
t_acc = 1 / (C0 * 1e4)  # Acceleration time scale
t_acc7 = tlc # Acceleration time scale
t_esc = t_acc  # Escape time
t_times = np.array([0.01, 0.1, 0.5, 1, 2, 3, 5, 10, 15, 20, 40, 60, 80]) * t_acc
def run_katarzynski():

    # Arrays
    t = np.logspace(0, np.log10(t_acc * 30), num_t)
    t = np.unique(np.concatenate((t, t_times)))
    numt = len(t)

    f = np.logspace(np.log10(fmin), np.log10(fmax), numf)
    g = np.logspace(np.log10(gmin), np.log10(gmax), numg)
    D_0 = 0.5 * np.power(g, 2) / t_acc
    D_07 = 0.5 * np.power(g, 2) / t_acc7
    D1 = 2 * D_0
    D2 = D1
    D3 = D1
    D4 = D1
    D5 = D1
    D6 = D1
    D7 = 2 * D_07
    gdot1 = np.zeros([numt, numg])
    gdot2 = np.zeros([numt, numg])
    gdot3 = np.zeros([numt, numg])
    gdot4 = np.zeros([numt, numg])
    gdot5 = np.zeros([numt, numg])
    gdot6 = np.zeros([numt, numg])
    gdot7 = np.zeros([numt, numg])
    gdot1[0, :] = C0 * np.power(g, 2) - 4 * D_0 / g
    gdot2[0, :] = C0 * np.power(g, 2) - 4 * D_0 / g
    gdot3[0, :] = C0 * np.power(g, 2) - 4 * D_0 / g
    gdot4[0, :] = C0 * np.power(g, 2) - 4 * D_0 / g
    gdot5[0, :] = C0 * np.power(g, 2) - 4 * D_0 / g
    gdot6[0, :] = C0 * np.power(g, 2) - 4 * D_0 / g
    gdot7[0, :] = (4 * cons.sigmaT* uB * np.power(g,2)/ (3 * cons.me * cons.cLight)) - (4 * D_07 / g)

    # Initial particle distribution
    p, gcut1, g2cut1 = 0, 1e0, 2e0
    p, gcut2, g2cut2 = 0, 1e6, 2e6
    p, gcut3, g2cut3 = 0, 1e1, 1e6
    p, gcut4, g2cut4 = 0, 1e0, 2e0
    p, gcut5, g2cut5 = 0, 1e6, 2e6
    p, gcut6, g2cut6 = 0, 1e1, 1e6
    p, gcut7, g2cut7 = 0, 1e0, 2e0
    n1 = np.zeros([numt, numg])
    n2 = np.zeros([numt, numg])
    n3 = np.zeros([numt, numg])
    n4 = np.zeros([numt, numg])
    n5 = np.zeros([numt, numg])
    n6 = np.zeros([numt, numg])
    n7 = np.zeros([numt, numg])
    n1[0, :] = misc_func.power_law(1, g, p, gcut1, g2cut1)
    n2[0, :] = misc_func.power_law(1, g, p, gcut2, g2cut2)
    n3[0, :] = misc_func.power_law(1, g, p, gcut3, g2cut3)
    n4[0, :] = misc_func.power_law(1, g, p, gcut4, g2cut4)
    n5[0, :] = misc_func.power_law(1, g, p, gcut5, g2cut5)
    n6[0, :] = misc_func.power_law(1, g, p, gcut6, g2cut6)
    n7[0, :] = misc_func.power_law(7, g, p, gcut7, g2cut7)

    j_s7 = np.zeros([numt, numf])  # synchrotron emissivity
    ambs7 = np.zeros([numt, numf])  # absorbtion coefficient
    j_ssc7 = np.zeros([numt, numf])  # ssc emissivity
    I_s7 = np.zeros([numt, numf])  # synchrotron Intensity
    I_ssc7 = np.zeros([numt, numf])  # ssc Intensity

    # Time loop
    for i in range(1, len(t)):
        dt = t[i] - t[i - 1]
        n1[i, :] = para.fp_findif_difu(dt, g, n1[i - 1, :], gdot1[i - 1, :], D1, np.zeros_like(D1), 1e200, tlc, False)
        n2[i, :] = para.fp_findif_difu(dt, g, n2[i - 1, :], gdot2[i - 1, :], D2, np.zeros_like(D2), 1e200, tlc, False)
        n3[i, :] = para.fp_findif_difu(dt, g, n3[i - 1, :], gdot3[i - 1, :], D3, np.zeros_like(D3), 1e200, tlc, False)
        n4[i, :] = para.fp_findif_difu(dt, g, n4[i - 1, :], gdot4[i - 1, :], D4, (4.7e4)*n4[0, :]/t_esc, t_esc, tlc, False)
        n5[i, :] = para.fp_findif_difu(dt, g, n5[i - 1, :], gdot5[i - 1, :], D5, (4.7e4)*n5[0, :]/t_esc, t_esc, tlc, False)
        n6[i, :] = para.fp_findif_difu(dt, g, n6[i - 1, :], gdot6[i - 1, :], D6, (4.7e4)*n6[0, :]/t_esc, t_esc, tlc, False)
        n7[i, :] = para.fp_findif_difu(dt, g, n7[i - 1, :], gdot7[i - 1, :], D7, np.zeros_like(D7), 1e200, tlc, False)

        j_s7[i, :], ambs7[i, :] = para.syn_emissivity_full(f, g, n7[i, :], B, with_abs)  # ,sync and absorb
        I_s7[i, :] = para.rad_trans_blob(R, j_s7[i, :], ambs7[i, :])
        j_ssc7[i, :] = para.ic_iso_powlaw_full(f, I_s7[i, :], g, n7[i, :])
        I_ssc7[i, :] = para.rad_trans_blob(R, j_ssc7[i, :], ambs7[i, :])

        dotgKN7 = para.rad_cool_pwl(g, f, 4 * np.pi * I_s7[i, :] / cons.cLight, cool_withKN)

        gdot1[i, :] = gdot1[0, :]
        gdot2[i, :] = gdot2[0, :]
        gdot3[i, :] = gdot3[0, :]
        gdot4[i, :] = gdot4[0, :]
        gdot5[i, :] = gdot5[0, :]
        gdot6[i, :] = gdot6[0, :]
        gdot7[i, :] = gdot7[0, :] + dotgKN7

    return [n1,n2,n3,n4,n5,n6,n7],g,[I_s7],[I_ssc7],f,t,t_times


def get_data(file_path):
    x, y = [], []
    with open(file_path, 'r') as f:
        for line in f.readlines():
            parts = line.split(", ")
            x.append(float(parts[0]))
            y.append(float(parts[1]))
    return x, y


def n_plot(n, g, t, plt_compare=False,compare_fig=None):
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')

    if compare_fig:
        if plt_compare:
            data_x, data_y = get_data(f"compare_data/katarzynski_digitized_plots/{compare_fig}")
            ax.plot(data_x, data_y,lw=3, marker='d', color='black', label='Comparison Data')
        if (compare_fig == 'fig1_a'):
            ax.set_ylim(1e-12, 1e1)
            ax.set_xlim(1e0, 1e7)
        elif (compare_fig == 'fig1_b'):
            ax.set_ylim(1e-10, 1e3)
            ax.set_xlim(1e0, 1e7)
        elif (compare_fig == 'fig1_c'):
            ax.set_ylim(1e-9, 1e4)
            ax.set_xlim(1e0, 1e7)
        elif (compare_fig == 'fig2_a'):
            ax.set_ylim(1e-10, 1e5)
            ax.set_xlim(1e0, 1e7)
        elif (compare_fig == 'fig2_b'):
            ax.set_ylim(1e-8, 1e7)
            ax.set_xlim(1e0, 1e7)
        elif (compare_fig == 'fig2_c'):
            ax.set_ylim(1e-8, 1e7)
            ax.set_xlim(1e0, 1e7)
        elif (compare_fig == 'fig3_a'):
            ax.set_ylim(1e-12, 1e2)
            ax.set_xlim(1e0, 1e8)
        else:
            ax.set_ylim(1e-9, 1e4)
            ax.set_xlim(1e0, 1e7)


    norm = matplotlib.colors.LogNorm(vmin=t[0], vmax=t[-1])
    cmap = cm.rainbow
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    plot_times = np.abs(np.subtract.outer(t, t_times)).argmin(0)
    for i in np.append(0,plot_times):
        ax.plot(g, n[i, :], color=cmap(norm(t[i])), lw=2, label=f'Time = {t[i]/t_acc:.2f} t_acc' if i % 20 == 0 else "")

    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Time [s]', rotation=270, labelpad=15)
    ax.set_xlabel(r"Lorentz Factor $\gamma$")
    ax.set_ylabel(r"Particle Distribution $n(\gamma)$ [cm$^{-3}$]")
    ax.set_title('Particle Distribution over Time')



    plt.tight_layout()
    plt.show()

def nuFnu_plots(Im, nu, t, compare_fig):
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')

    norm = matplotlib.colors.LogNorm(vmin=t[0], vmax=t[-1])
    cmap = cm.rainbow
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=matplotlib.colors.LogNorm(vmin=t[0], vmax=t[-1]))
    sm.set_array([])

    t_normalized = t / t_acc
    specified_times = [0.01, 0.1, 0.5, 1, 2, 3, 5, 10, 15, 20, 40, 60, 80]
    indices_to_plot = [np.argmin(np.abs(t_normalized - time)) for time in specified_times]

    max_luminosity = None
    dop = 21
    for i, time in enumerate(t_normalized):
        if i not in indices_to_plot and i not in (0, len(t) - 1) and i % 10 != 0:
            continue
        nu_Fnu = (dop**4)*nu*(Im[i,:])* (4 * np.pi * (R**2))/ (4 * np.pi * (4.32e26)**2)
        ax.plot(dop*nu, nu_Fnu, color=cmap(norm(t[i])), label=f'Time = {time:.2f} t_acc')

        # Update the maximum luminosity value
        current_max = np.max(nu_Fnu)
        max_luminosity = current_max if max_luminosity is None else max(max_luminosity, current_max)

    ax.set_xlim(nu[0], nu[-1])
    ax.set_ylim(1e-16, 1e-5) if compare_fig == 'fig3_b' else ax.set_ylim(1e-16, 2 * max_luminosity)
    ax.set_xlabel(r"$\nu$ [Hz]", fontsize=18)
    ax.set_ylabel(r"$\nu F_{\nu}$ $[\frac{erg}{s cm^{2}}]$", fontsize=18)

    colorbar_ticks = np.logspace(0, 7, 8)
    cbar = plt.colorbar(sm, ax=ax, anchor=(-0.6, 0.0), ticks=colorbar_ticks)
    cbar.ax.set_yticklabels([f"$10^{{{int(np.log10(tick))}}}$" for tick in colorbar_ticks])
    cbar.ax.set_ylabel(r"t [s]", fontsize=18)

    ax.tick_params(axis='both', which='both', labelsize=15)
    ax.tick_params(axis='both', which='major', size=10)
    ax.tick_params(axis='y', which='minor', size=5)

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    ns,g,I_syns,I_sscs,nu,t,t_times = run_katarzynski()
    compare_figs = ['fig1_a','fig1_b','fig1_c','fig2_a','fig2_b','fig2_c','fig3_a']
    compare_nfnu_figs = ['fig3_b','fig3_c']
    for i,n in enumerate(ns):
        if(i<6):
            continue
        n_plot(n, g, t,plt_compare=False, compare_fig=compare_figs[i])
        if(i>5):
            Is = I_syns[i - 6]
            Issc = I_sscs[i - 6]
            nuFnu_plots(Is, nu, t, compare_fig=compare_nfnu_figs[i - 6])
            nuFnu_plots(Issc, nu, t, compare_fig=compare_nfnu_figs[i - 5])