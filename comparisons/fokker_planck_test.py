import numpy as np
import matplotlib.pyplot as plt
from tleco import misc_func
from dependencies import fp_test
import tleco as tl


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
    x_values = np.logspace(np.log10(gmin), np.log10(gmax), num=10000)
    y_values = x_values ** p
    return np.trapz(y_values, x_values)


def plot_setup():
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_yscale('log')
    ax.set_xscale('log')
    return fig, ax


def add_legend(ax, location='lower center'):
    scale_mult = 2
    ax.legend(loc='best', fontsize=12 * scale_mult, title_fontsize=12)
    # ax.set_title(" Steady State Solutions ", y=1.0, pad=-14, loc='left')
    # ax.tick_params(axis='both', which='both', labelsize=15)


def plot_data(ax, x_data, y_data, plot_type='scatter', **kwargs):
    if plot_type == 'scatter':
        ax.scatter(x_data, y_data, **kwargs)
    elif plot_type == 'plot':
        ax.plot(x_data, y_data, **kwargs)


def adjust_plot(ax, x_label, y_label, title):
    scale_mult = 3
    ax.set_xlabel(x_label,fontsize=10 * scale_mult)
    ax.set_ylabel(y_label, fontsize=10 * scale_mult)
    ax.set_title(title, y=1.0, pad=-14, loc='left')
    ax.tick_params(axis='both', which='major', size=12 * scale_mult, labelsize=12 * scale_mult)
    ax.tick_params(axis='both', which='minor', size=0 * scale_mult)


def get_error_analytic(ei, eig, n0):
    ef = analytical_solution(n0, eig)
    valid_ef = ef > 0
    errors = np.where(valid_ef, ((ef - ei) / ef) ** 2, 0)
    er = np.sqrt(sum(errors)/len(ei))
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
        n1[i, :] = tl.fp_findif_difu(dt, g, n1[i - 1, :], gdot1[i - 1, :], D1, np.zeros_like(D1), 1e200, tlc, False)
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

    xc2 = np.logspace(1, 3.8)
    ax_error.plot(xc2, ys[2] * (xc2 / numgs[2]) ** -2, label='$\propto N^{-2}$',linestyle='--', color='black',linewidth=4)

    plot_data(ax_error, xs, ys, 'scatter',c='blue',s=250,marker='o')
    allowed_numgs=[25,50,100,800]
    for i in range(len(gs)):
        plot_data(ax_eers, g2s[i], eers[i], 'scatter', s=4)
        if(numgs[i] in allowed_numgs):
            plot_data(ax_ns, gs[i], ns[i][-1, :], 'plot', label='N=' + str(numgs[i]),linewidth=4)
    ax_ns.plot(gs[-1], analytical_solution(n0s[-1], gs[-1]), '--',label='Solution',color='k', linewidth=4)

    adjust_plot(ax_error, "Number of bins (N)", "Error", "")
    adjust_plot(ax_ns, "$\gamma$", r"n $[cm^{-3}]$", "")
    adjust_plot(ax_eers, "$\gamma$", "Error", "Error Analysis")
    # adjust_plot(ax_nios, "Number of bins (N)", "Normalized Intensity", "Intensity Over Bins")

    ax_ns.set_xlim([1e0,1e6])
    ax_ns.set_ylim([1e-8,1e2])

    ax_error.set_xlim([1e1, 1e3])
    ax_error.set_ylim([1e-3, 1e1])

    add_legend(ax_ns)
    add_legend(ax_error)
    fig_ns.savefig("Figs/numg_convergence_n.pdf", dpi=800, bbox_inches="tight")
    fig_ns.savefig("Figs/numg_convergence_n.png", dpi=800, bbox_inches="tight")
    fig_error.savefig("Figs/numg_convergence_rate.pdf", dpi=800, bbox_inches="tight")
    fig_error.savefig("Figs/numg_convergence_rate.png", dpi=800, bbox_inches="tight")
    # plt.show()


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
        n1[i, :] = tl.fp_findif_difu(dt, g, n1[i - 1, :], gdot1[i - 1, :], D1, np.zeros_like(D1), 1e0, tlc, False)
        gdot1[i, :] = gdot1[0, :]

    return g, n1, t


def get_error_analytic_time(ei, ef):
    errors = np.where(ef > 0, ((ef - ei) / ef) ** 2, 0)
    eers = np.sqrt(errors)
    er = np.sqrt(np.sum(errors)/len(ei))

    return er, eers


def plot_convergence_results_time(numts, numg, tt):
    results = []
    nios = []
    gmin, gmax = 1e4, 1e6
    n0 = 1  # get_n0(0, gmin, gmax)
    for mt in numts:
        g, n, tr = run_time_test(mt, numg)
        tind = np.argmin(np.abs(tr - tt))
        ei = n[tind, :]
        norm = np.trapz(fp_test.eq_59_park1995(1e-4, g)*g*np.log(10),np.log10(g))
        yt = np.array(fp_test.eq_59_park1995(tt + 1e-4, g))
        norm2 = np.trapz(yt,g)
        ef = n0 * np.array(fp_test.eq_59_park1995(tt + 1e-4, g)) / norm
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

    plot_results(results, tt, nios)


def plot_results(results, tt, nios, separate_figures=True):
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']

    if separate_figures:
        figs = [plt.figure(figsize=(8, 8)) for _ in range(4)]
        axs = [fig.add_subplot(111) for fig in figs]
    else:
        fig, axs = plt.subplots(4, 1, figsize=(10, 20))
    mts = []
    ers = []
    allowed_numts =[25,50,100,800]
    for i, (mt, er, eer, eig, ei,ef) in enumerate(results):
        mts.append(mt)
        ers.append(er)
        axs[0].scatter(mt, er, color='blue',s=250)
        axs[1].scatter(eig, eer, color=colors[i % len(colors)], s=4)
        if (mt in allowed_numts):
            axs[2].plot(eig, ei, label='N=' + str(mt), linewidth=4)
        axs[3].scatter(mt, nios[i], color=colors[i % len(colors)])

    x_sol = results[-1][3]
    y_sol = results[-1][5]
    scale_mult = 3
    axs[3].set_xscale('log')
    axs[0].set_title(f"  t={tt:.2E}", y=0.98, pad=-14, loc='left', fontsize=10 * scale_mult)
    axs[2].set_title(f"  t={tt:.2E}", y=0.98, pad=-14, loc='left', fontsize=10 * scale_mult)
    axs[0].set_xlabel('Number of bins (N)', fontsize=10 * scale_mult)
    axs[0].set_ylabel('Error', fontsize=10 * scale_mult)
    axs[2].set_xlabel("$\gamma$", fontsize=10 * scale_mult)
    axs[2].set_ylabel(r"n $[cm^{-3}]$", fontsize=10 * scale_mult)
    axs[2].plot(x_sol,  y_sol, '--',label='Solution',color='k',linewidth=4)
    axs[0].plot(np.logspace(1,3), ers[2] * (np.logspace(1,3) /mts[2]) ** -1,'--',linewidth=4, label='$\propto N^{-1}$', color='black')
    axs[2].set_ylim([1e-8, 5e-2])
    axs[2].set_xlim([1e1, 1e4])
    axs[0].set_xlim([10, 1e3])
    axs[0].set_ylim([min(ers)/2, max(ers)*2])



    for ax in axs[:3]:
        ax.tick_params(axis='both', which='major', size=12 * scale_mult, labelsize=12 * scale_mult)
        ax.tick_params(axis='both', which='minor', size=0 * scale_mult, labelsize=0 * scale_mult)
        ax.set_yscale('log')
        ax.set_xscale('log')
        scale_mult = 2
        ax.legend(loc='best', fontsize=12 * scale_mult, title_fontsize=12)

    # if separate_figures:
    #     for fig in figs:
    #         fig.savefig()

    figs[2].savefig(f"Figs/numt_{tt:.2E}_convergence_n.pdf", dpi=800, bbox_inches="tight")
    figs[2].savefig(f"Figs/numt_{tt:.2E}_convergence_n.png", dpi=800, bbox_inches="tight")
    figs[0].savefig(f"Figs/numt_{tt:.2E}_convergence_rate.pdf", dpi=800, bbox_inches="tight")
    figs[0].savefig(f"Figs/numt_{tt:.2E}_convergence_rate.png", dpi=800, bbox_inches="tight")

    # plt.show()

if __name__ == '__main__':
    # garr = [25,50, 100,200,400,800]
    # numtarr = [800]
    # convergence_plots_analytic(numtarr[0], garr)
    numtarr = [25,50, 100,200,400,800]
    garr = [800]
    ts = [5e-3, 1e-2]
    for t in ts:
        plot_convergence_results_time(numtarr, garr[0],t)
