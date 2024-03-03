import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors
import misc_func
import paramo as para
import constants as cons

def run_katarzynski():
    # Constants
    numt, numg, numf = 200, 150, 80
    fmin, fmax, tmax, gmin, gmax = 1e8, 1e28, 1e18, 1e0, 1e8
    R, B = 3.2e15, 0.05  # Blob size and magnetic field
    uB = (B ** 2) / (np.pi * 8)  # Magnetic field energy density
    C0 = 3.48e-11  # Energy density constant
    t_acc = 1 / (C0 * 1e4)  # Acceleration time scale
    t_esc = t_acc  # Escape time
    t_inj = R / cons.cLight  # Injection time
    tlc = t_inj  # Light crossing time
    # Arrays
    t = np.logspace(0, np.log10(t_acc * 30), numt)
    t_times = np.array([0.01, 0.1, 0.5, 1, 2, 3, 5, 10, 15, 20, 40, 60, 80]) * t_acc
    t = np.unique(np.concatenate((t, t_times)))
    numt = len(t)

    f = np.logspace(np.log10(fmin), np.log10(fmax), numf)
    g = np.logspace(np.log10(gmin), np.log10(gmax), numg)
    D_0 = 0.5 * np.power(g, 2) / t_acc
    D1 = 2 * D_0
    D2 = D1
    D3 = D1
    D4 = D1
    gdot1 = np.zeros([numt, numg])
    gdot2 = np.zeros([numt, numg])
    gdot3 = np.zeros([numt, numg])
    gdot4 = np.zeros([numt, numg])
    gdot1[0, :] = C0 * np.power(g, 2) - 4 * D_0 / g
    gdot2[0, :] = C0 * np.power(g, 2) - 4 * D_0 / g
    gdot3[0, :] = C0 * np.power(g, 2) - 4 * D_0 / g
    gdot4[0, :] = C0 * np.power(g, 2) - 4 * D_0 / g

    # Initial particle distribution
    p, gcut1, g2cut1 = 0, 1e0, 2e0
    p, gcut2, g2cut2 = 0, 1e6, 2e6
    p, gcut3, g2cut3 = 0, 1e1, 1e6
    p, gcut4, g2cut4 = 0, 1e0, 2e0
    n1 = np.zeros([numt, numg])
    n2 = np.zeros([numt, numg])
    n3 = np.zeros([numt, numg])
    n4 = np.zeros([numt, numg])
    n1[0, :] = misc_func.power_law(1, g, p, gcut1, g2cut1)
    n2[0, :] = misc_func.power_law(1, g, p, gcut2, g2cut2)
    n3[0, :] = misc_func.power_law(1, g, p, gcut3, g2cut3)
    n4[0, :] = misc_func.power_law(1, g, p, gcut4, g2cut4)

    # Time loop
    for i in range(1, len(t)):
        dt = t[i] - t[i - 1]
        n1[i, :] = para.fp_findif_difu(dt, g, n1[i - 1, :], gdot1[i - 1, :], D1, np.zeros_like(D1), 1e200, tlc, False)
        n2[i, :] = para.fp_findif_difu(dt, g, n2[i - 1, :], gdot2[i - 1, :], D2, np.zeros_like(D2), 1e200, tlc, False)
        n3[i, :] = para.fp_findif_difu(dt, g, n3[i - 1, :], gdot3[i - 1, :], D3, np.zeros_like(D2), 1e200, tlc, False)
        n4[i, :] = para.fp_findif_difu(dt, g, n4[i - 1, :], gdot4[i - 1, :], D4, (4.7e4)*n4[0, :]/t_esc, t_esc, tlc, False)
        gdot1[i, :] = gdot1[0, :]
        gdot2[i, :] = gdot2[0, :]
        gdot3[i, :] = gdot3[0, :]
        gdot4[i, :] = gdot4[0, :]

    return [n1,n2,n3,n4],g,t,t_acc,t_times


def get_data(file_path):
    x, y = [], []
    with open(file_path, 'r') as f:
        for line in f.readlines():
            parts = line.split(", ")
            x.append(float(parts[0]))
            y.append(float(parts[1]))
    return x, y


def n_plot(n, g, t, t_acc, plt_compare=False,compare_fig=None):
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




if __name__ == '__main__':
    ns,g,t,t_acc,t_times = run_katarzynski()
    compare_figs = ['fig1_a','fig1_b','fig1_c','fig2_a']
    for i,n in enumerate(ns):
        n_plot(n, g, t, t_acc,plt_compare=False, compare_fig=compare_figs[i])