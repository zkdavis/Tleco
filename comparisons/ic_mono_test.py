import matplotlib.pyplot as plt
import tleco as tl
import numpy as np
from tleco import constants as C
from dependencies import radiation_test as rt

def run_test(num_g,num_f):
    p = 2.0
    n0 =1
    nu_in = 1e13
    eps_in = nu_in*C.hPlanck/C.energy_e
    u0 = 1e0
    nu_s = np.logspace(8, 25, num_f)
    eps_s = nu_s * C.hPlanck / C.energy_e
    g_array = np.logspace(1, 3, num_g)
    n_array = n0*(g_array ** -p) / np.trapz(g_array ** -p, g_array)
    j_ic = rt.j_ic_mono_pwl_electron_dermer_full(eps_s, eps_in, u0, p, n_array, g_array)
    j_ic_para = np.array(tl.ic_iso_monochrome_full(nu_s, u0, nu_in, n_array, g_array))

    return nu_s,j_ic, j_ic_para


def run_convergence_test(num_g_values,num_f_values,nu_bounds):
    scale_mult = 3
    errors = np.zeros((len(num_g_values), len(num_f_values)))
    errors_spec = np.zeros((len(num_g_values), len(num_f_values)),dtype=object)
    nu_array = []
    for i, num_g in enumerate(num_g_values):
        for j, num_f in enumerate(num_f_values):
            nu_s, j_ic, j_ic_para = run_test(num_g=num_g, num_f=num_f)
            nu_min_i = np.argmin(np.abs(nu_s - nu_bounds[0]))
            nu_max_i = np.argmin(np.abs(nu_s - nu_bounds[1]))
            nu_s = nu_s[nu_min_i:nu_max_i]
            j_ic = j_ic[nu_min_i:nu_max_i]
            j_ic_para = j_ic_para[nu_min_i:nu_max_i]
            non_zero_indices = j_ic != 0
            j_ic_filtered = j_ic[non_zero_indices]
            j_ic_para_filtered = j_ic_para[non_zero_indices]
            nu_array.append(nu_s[non_zero_indices])
            er_spec = np.abs((j_ic_para_filtered - j_ic_filtered)/j_ic_filtered)
            re = np.sum(er_spec)/len(er_spec)
            errors[i, j] = re
            errors_spec[i, j] = er_spec

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(num_g_values, errors[:, -1],label='$\gamma$ convergence', marker='o', linestyle='-', color='blue')
    ax.plot(num_f_values, errors[-1, :], label='$\\nu$ convergence',marker='o', linestyle='-', color='red')
    # ax.set_title('Convergence')
    ax.set_xlabel('bin size', fontsize=10 * scale_mult)
    ax.set_ylabel('Relative Error', fontsize=10 * scale_mult)



    # ax.set_title('Plot Title', fontsize=18*scale_mult)

    ax.tick_params(axis='both', which='major', size=8 * scale_mult, labelsize=12 * scale_mult)
    ax.tick_params(axis='both', which='minor', size=8 * scale_mult, labelsize=6 * scale_mult)

    ax.legend(loc='best', fontsize=10 * scale_mult, title_fontsize=12)
    # ax.grid(True)
    # fig.savefig("error_plot.pdf", dpi=200, bbox_inches="tight")

    fig2, ax2 = plt.subplots(figsize=(16, 12))

    ax2.set_xscale("log")
    ax2.set_yscale("log")

    # Plot the error
    for i, ei in enumerate(errors_spec[-1, :]):
        ax2.plot(nu_array[i], ei, label='$\\frac{\\text{|PARAMO - Dermer|}}{\\text{|Dermer|}}$' + ' at numf={}'.format(num_f_values[i]), linewidth=2 * scale_mult, linestyle=':')

    ax2.set_xlim([2e12, 1e20])
    ax2.set_xlabel('$\\nu$ [Hz]', fontsize=15 * scale_mult)
    ax2.set_ylabel('Relative Error', fontsize=15 * scale_mult)
    # ax.set_title('Plot Title', fontsize=18*scale_mult)

    ax2.tick_params(axis='both', which='major', size=12 * scale_mult, labelsize=12 * scale_mult)
    ax2.tick_params(axis='both', which='minor', size=0 * scale_mult)

    ax2.legend(loc='best', fontsize=12 * scale_mult, title_fontsize=12)

    # fig2.savefig("error_plot.pdf", dpi=200, bbox_inches="tight")

    plt.show()


def ic_mono_get_error(num_g,num_f,nu_bounds):
    nu_s, j_ic, j_ic_para = run_test(num_g=num_g, num_f=num_f)
    nu_min_i = np.argmin(np.abs(nu_s - nu_bounds[0]))
    nu_max_i = np.argmin(np.abs(nu_s - nu_bounds[1]))
    nu_s = nu_s[nu_min_i:nu_max_i]
    j_ic = j_ic[nu_min_i:nu_max_i]
    j_ic_para = j_ic_para[nu_min_i:nu_max_i]
    non_zero_indices = j_ic != 0
    j_ic_filtered = j_ic[non_zero_indices]
    j_ic_para_filtered = j_ic_para[non_zero_indices]
    er_spec = np.abs((j_ic_para_filtered - j_ic_filtered) / j_ic_filtered)
    return np.sum(er_spec) / len(er_spec)


def plot_comparison_and_error():
    num_g = 300#150
    num_f = 300#150
    scale_mult = 4
    nu_s, j_ic, j_ic_para = run_test(num_g=num_g, num_f=num_f)

    print(np.trapz(j_ic,nu_s))

    fig, ax = plt.subplots(figsize=(16, 12))

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.plot(nu_s, j_ic_para, label='PARAMO', linewidth=3 * scale_mult, color='red')
    ax.plot(nu_s, j_ic, label='Dermer', linewidth=3 * scale_mult, linestyle='--', color='blue')

    ax.set_xlim([2e12, 1e20])
    ax.set_ylim([1e-32, 1e-27])
    ax.set_xlabel('$\\nu$ [Hz]', fontsize=15 * scale_mult)
    ax.set_ylabel('$j_\\nu$ [$\\frac{ergs}{s \ Hz \ cm^3}$]', fontsize=15 * scale_mult)
    # ax.set_title('Plot Title', fontsize=18*scale_mult)

    ax.tick_params(axis='both', which='major', size=12 * scale_mult, labelsize=12 * scale_mult)
    ax.tick_params(axis='both', which='minor', size=0 * scale_mult)

    ax.legend(loc='upper right', fontsize=12 * (scale_mult - 1), title_fontsize=12)

    fig.savefig("Figs/mono_comparison.pdf", dpi=800, bbox_inches="tight")
    fig.savefig("Figs/mono_comparison.png", dpi=800, bbox_inches="tight")

    error = np.abs((j_ic_para - j_ic) / j_ic)

    fig2, ax2 = plt.subplots(figsize=(16, 12))

    ax2.set_xscale("log")
    ax2.set_yscale("log")

    ax2.plot(nu_s, error, label='$\\frac{\\text{|PARAMO - Dermer|}}{\\text{|Dermer|}}$', linewidth=2 * scale_mult, linestyle=':', color='green')

    ax.set_xlim([2e12, 1e20])
    ax2.set_xlabel('$\\nu$ [Hz]', fontsize=15 * scale_mult)
    ax2.set_ylabel('$Relative Error', fontsize=15 * scale_mult)
    # ax.set_title('Plot Title', fontsize=18*scale_mult)

    ax2.tick_params(axis='both', which='major', size=12 * scale_mult, labelsize=12 * scale_mult)
    ax2.tick_params(axis='both', which='minor', size=0 * scale_mult)

    ax2.legend(loc='best', fontsize=12 * scale_mult, title_fontsize=12)

    # fig2.savefig("error_plot.pdf", dpi=200, bbox_inches="tight")

    # plt.show()

if __name__ == '__main__':
    plot_comparison_and_error()
    # run_convergence_test([25,50,150,300,500,1000],[25,50,150,300,500,1000],[2e13,2.5e19])