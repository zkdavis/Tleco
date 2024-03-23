import numpy as np
import matplotlib.pyplot as plt
from paramo import misc_func
from dependencies import fp_test
import paramo as para


def analytical_form(g):
    C0 = 3.48e-11
    tesc = 1e0 / (C0 * (10 ** (4.5)))
    return (g ** 2) * np.exp(-2 * C0 * tesc * (g - 1))


def analytical_solution(n0, g):
    integrand_values = analytical_form(g)
    integral_approx = np.trapz(integrand_values, g)
    x = n0 / integral_approx
    return x * analytical_form(g)


def get_n0(p, gmin, gmax):
    x_values = np.logspace(np.log10(gmin), np.log10(gmax), num=1000)
    y_values = x_values ** p
    return np.trapz(y_values, x_values)


def plot_setup():
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')
    return fig, ax


def add_legend(ax, location='lower center'):
    ax.legend(loc=location)
    ax.set_xlabel("$\gamma$", fontsize=18)
    ax.set_ylabel(r"n $[cm^{-3}]$", fontsize=18)
    ax.set_title(" Steady State Solutions ", y=1.0, pad=-14, loc='left')
    ax.tick_params(axis='both', which='both', labelsize=15)


def plot_data(ax, x_data, y_data, plot_type='scatter', **kwargs):
    if plot_type == 'scatter':
        ax.scatter(x_data, y_data, **kwargs)
    elif plot_type == 'plot':
        ax.plot(x_data, y_data, **kwargs)


def adjust_plot(ax, x_label, y_label, title):
    ax.set_xlabel(x_label, fontsize=18)
    ax.set_ylabel(y_label, fontsize=18)
    ax.set_title(title, y=1.0, pad=-14, loc='left')
    ax.tick_params(axis='both', which='both', labelsize=15)


def get_error_analytic(ei, eig, n0):
    ef = analytical_solution(n0, eig)
    valid_ef = ef > 0
    errors = np.where(valid_ef, ((ef - ei) / ef) ** 2, 0)
    er = np.sqrt(sum(errors) / len(eig))
    eers = np.sqrt(errors)
    return er, eers


def run_test(numg, numt, n0):
    C0 = 3.48e-11
    t_acc = 1 / (C0 * (10 ** (4.5)))
    tmax = t_acc * 300
    t = np.logspace(0, np.log10(tmax), numt)
    g = np.logspace(0, np.log10(1.5e8), numg)
    D_0 = 0.5 * np.power(g, 2) / t_acc
    D1 = 2 * D_0
    gdot1 = np.zeros([numt, numg])
    gdot1[0, :] = C0 * np.power(g, 2) - 4 * D_0 / g
    tlc = tmax
    p, gcut1, g2cut1 = 0, 1e4, 1e6
    n1 = np.zeros([numt, numg])
    n1[0, :] = misc_func.power_law(1, g, p, gcut1, g2cut1)
    n1[0, :] = n0 * n1[0, :] / np.trapz(n1[0, :], g)

    for i in range(1, numt):
        dt = t[i] - t[i - 1]
        n1[i, :] = para.fp_findif_difu(dt, g, n1[i - 1, :], gdot1[i - 1, :], D1, np.zeros_like(D1), 1e200, tlc, False)
        gdot1[i, :] = gdot1[0, :]

    return g, n1


def process_data(numt, numgs):
    xs, ys, gs, g2s, ns, nios, n0s, eers = [], [], [], [], [], [], [], []
    n0 = get_n0(0, 1e4, 1e6)
    for i, mg in enumerate(numgs):
        g, n = run_test(mg, numt, n0)
        gminci, gmaxci = np.argmin(np.abs(g - 1e1)), np.argmin(np.abs(g - 1e7))
        g2, n2 = g[gminci:gmaxci], n[-1, gminci:gmaxci]
        er, eer = get_error_analytic(n2, g2, n0)
        xs.append(mg)
        gs.append(g)
        ns.append(n)
        eers.append(eer)
        ys.append(er)
        n0s.append(n0)
        g2s.append(g2)
    return xs, ys, gs, g2s, ns, n0s, eers


def convergence_plots_analytic(numt, numgs):
    xs, ys, gs, g2s, ns, n0s, eers = process_data(numt, numgs)

    # Setup plots
    fig_error, ax_error = plot_setup()
    fig_ns, ax_ns = plot_setup()
    fig_eers, ax_eers = plot_setup()
    fig_nios, ax_nios = plot_setup()
    ax_nios.set_xscale('linear')

    xc2 = np.logspace(1.5, 3.8)
    ax_error.plot(xc2, ys[2] * (xc2 / xs[2]) ** -1, label='$\propto N^{-1}$', color='black')

    plot_data(ax_error, numgs, ys, 'scatter')
    for i in range(len(gs)):
        plot_data(ax_eers, g2s[i], eers[i], 'scatter', s=4)
        plot_data(ax_ns, gs[i], ns[i][-1, :], 'plot', label='N=' + str(numgs[i]))
    ax_ns.plot(gs[-1], analytical_solution(n0s[-1], gs[-1]), '--')

    adjust_plot(ax_error, "Number of bins (N)", "Error", "Convergence at steady state")
    adjust_plot(ax_ns, "$\gamma$", r"n $[cm^{-3}]$", "Steady State Solutions")
    adjust_plot(ax_eers, "$\gamma$", "Error", "Error Analysis")
    # adjust_plot(ax_nios, "Number of bins (N)", "Normalized Intensity", "Intensity Over Bins")

    add_legend(ax_ns)
    plt.show()


def run_time_test(numt,numg):
    t_acc = 1
    tmax = t_acc * 0.1
    t = np.logspace(np.log10(1e-4), np.log10(tmax), numt)
    g = np.logspace(0, np.log10(1.5e8), numg)
    D_0 = np.power(g, 3)
    D1 = 2 * D_0
    gdot1 = np.zeros([numt, numg])
    gdot1[0, :] = (g ** 2) - 3 * (g ** 2) - 2 * D_0 / g
    tlc = tmax * 100
    n1 = np.zeros([numt, numg])
    n1[0, :] = fp_test.eq_59_park1995(1e-4, g)
    n1[0, :] = n1[0, :] / np.trapz(n1[0, :], g)

    for i in range(1, numt):
        dt = t[i] - t[i - 1]
        n1[i, :] = para.fp_findif_difu(dt, g, n1[i - 1, :], gdot1[i - 1, :], D1, np.zeros_like(D1), 1e0, tlc, False)
        gdot1[i, :] = gdot1[0, :]

    return g, n1, t


def get_error_analytic_time(ei, ef):
    errors = np.where(ef > 0, ((ef - ei) / ef) ** 2, 0)
    eers = np.sqrt(errors)
    er = np.sqrt(np.sum(errors)/len(ei))

    return er, eers


def plot_convergence_results_time(numts, numg, t):
    results = []
    nios = []
    gmin, gmax = 1e4, 1e6
    n0 = 1  # get_n0(0, gmin, gmax)
    ts = [5e-4, 1e-3, 5e-3, 1e-2]
    for mt in numts:
        g, n, tr = run_time_test(mt, numg)
        tind = np.argmin(np.abs(tr - t))
        ei = n[tind, :]
        norm = np.trapz(fp_test.eq_59_park1995(1e-4, g)*g*np.log(10),np.log10(g))
        yt = np.array(fp_test.eq_59_park1995(t + 1e-4, g))
        norm2 = np.trapz(yt,g)
        ef = n0 * np.array(fp_test.eq_59_park1995(t + 1e-4, g)) / norm
        nio = np.trapz(ei, g) / np.trapz(n[0, :], g)
        nios.append(nio)
        minargs = np.argsort(np.abs(ei - 1e-8))
        ma1 = minargs[0]
        ma2 = minargs[1]
        for k in range(len(minargs)):
            if (np.abs(minargs[k] - ma1) > 10):
                ma2 = minargs[k]
                break
        if (ma1 > ma2):
            gmaxci = ma1
            gminci = ma2
        else:
            gmaxci = ma2
            gminci = ma1

        g = g[gminci:gmaxci]
        ei = ei[gminci:gmaxci]
        ef = ef[gminci:gmaxci]

        er, eer = get_error_analytic_time(ei, ef)
        results.append((mt, er, eer, g, ei,ef))

    plot_results(results, ts[-1], nios)


def plot_results(results, t, nios, separate_figures=True):
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']

    if separate_figures:
        figs = [plt.figure(figsize=(10, 5)) for _ in range(4)]
        axs = [fig.add_subplot(111) for fig in figs]
    else:
        fig, axs = plt.subplots(4, 1, figsize=(10, 20))
    mts = []
    ers = []
    for i, (mt, er, eer, eig, ei,ef) in enumerate(results):
        mts.append(mt)
        ers.append(er)
        axs[0].scatter(mt, er, color=colors[i % len(colors)], label=f'MT={mt}')
        axs[1].scatter(eig, eer, color=colors[i % len(colors)], s=4)
        axs[2].plot(eig, ei, label=f'Number of bins: {mt}')
        axs[3].scatter(mt, nios[i], color=colors[i % len(colors)])

    x_sol = results[-1][3]
    y_sol = results[-1][5]

    axs[3].set_xscale('log')
    axs[0].set_title(f"Convergence at t={t:.2E}")
    axs[2].plot(x_sol,  y_sol, '--',
                label='Analytic solution')
    axs[0].plot(mts, ers[2] * (np.array(mts) /mts[2]) ** -1, label='$\propto N^{-1}$', color='black')
    axs[2].set_ylim([1e-8, 2e-2])
    axs[2].set_xlim([1e1, 1e4])

    for ax in axs[:3]:
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.legend()

    # if separate_figures:
    #     for fig in figs:
    #         fig.savefig()

    plt.show()

if __name__ == '__main__':
    # garr = [25, 50, 100, 200, 400, 800, 1600, 3200]
    # # garr = [100, 500, 3500]
    # # garr = [20,50, 100,200,300,500]
    # numtarr = [600]
    # # numtarr = [200]
    #
    # convergence_plots_analytic(numtarr[0], garr)

    # numtarr = [100, 300, 900]
    numtarr = [100, 300,500]
    # numtarr = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500]
    garr = [900]
    # garr = [200]
    ts = [5e-4, 1e-3, 5e-3, 1e-2]
    t = ts[3]
    plot_convergence_results_time(numtarr, garr[0],t)
