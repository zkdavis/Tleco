import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate


def get_n0(p, gmin, gmax):
    return integrate.quad(lambda x: x ** p, gmin, gmax)[0]


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


def process_data(numt, numgs, n0):
    xs, ys, gs, ns, nios, eers = [], [], [], [], [], []
    for i, mg in enumerate(numgs):
        eissr = get_convergence_results(numt, mg)
        ei, eig = eissr.n4[:, -1], eissr.g
        nio = np.trapz(ei, eig) / np.trapz(eissr.n4[:, 0], eig)
        nios.append(nio)
        gminci, gmaxci = np.argmin(np.abs(eig - 1e0)), np.argmin(np.abs(eig - 1e7))
        eig, ei = eig[gminci:gmaxci], ei[gminci:gmaxci]
        _, eer = get_error_analytic(ei, eig, n0)
        xs.append(len(eig))
        gs.append(eig)
        ns.append(ei)
        eers.append(eer)
        ys.append(_)
    return xs, ys, gs, ns, nios, eers


def convergence_plots_analytic(numt, numgs):
    p = 0e0
    gmin, gmax = 1e4, 1e6
    n0 = get_n0(p, gmin, gmax)
    xs, ys, gs, ns, nios, eers = process_data(numt, numgs, n0)

    # Setup plots
    fig_error, ax_error = plot_setup()
    fig_ns, ax_ns = plot_setup()
    fig_eers, ax_eers = plot_setup()
    fig_nios, ax_nios = plot_setup()
    ax_nios.set_xscale('linear')  # Only change for ax_nios

    # Plot data
    plot_data(ax_error, xs, ys, 'scatter')
    for i in range(len(gs)):
        plot_data(ax_eers, gs[i], eers[i], 'scatter', s=4)
        plot_data(ax_ns, gs[i], ns[i], 'plot', label='N=' + str(xs[i]))
    ax_ns.plot(gs[-1], analytical_solution(n0, gs[-1]), '--')

    # Adjust plots
    adjust_plot(ax_error, "Number of bins (N)", "Error", "Convergence at steady state")
    adjust_plot(ax_ns, "$\gamma$", r"n $[cm^{-3}]$", "Steady State Solutions")
    adjust_plot(ax_eers, "$\gamma$", "Error", "Error Analysis")
    adjust_plot(ax_nios, "Number of bins (N)", "Normalized Intensity", "Intensity Over Bins")

    # Legend and show
    add_legend(ax_ns)
    plt.show()
