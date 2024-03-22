import matplotlib.pyplot as plt
import paramo as para
import numpy as np
from paramo import constants as C
from ic_dependencies import radiation_test as rt

def run_test(num_g,num_f):
    p = 2
    n0 =1
    sig_T = 0.0001
    Temp = sig_T*C.energy_e/C.kBoltz
    u_bb = 1e0
    nu_s = np.logspace(8, 25, num_f)
    nu_s2 = np.logspace(8, 25, num_f)
    eps_s = nu_s2 * C.hPlanck / C.energy_e
    g_array = np.logspace(1, 3, num_g)
    n_array = n0*(g_array ** -p) / np.trapz(g_array ** -p, g_array)

    def ub(epsa):
        return rt.BB_photon_density(epsa,Temp)*epsa*u_bb/np.trapz(rt.BB_photon_density(epsa,Temp)*epsa,epsa)

    def ub_dermer(epsa):
        epsa = epsa * C.energy_e
        return rt.BB_photon_density(epsa, Temp) * epsa * u_bb / np.trapz(rt.BB_photon_density(epsa, Temp) * epsa, epsa)


    j_ic = rt.j_ic_iso_full_dermer(eps_s, eps_s,ub_dermer, n_array, g_array)

    j_ic_para = np.array(para.ic_iso_powlaw_full(nu_s, C.cLight * ub(nu_s*C.hPlanck), g_array, n_array))

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

    ax2.set_xlim([5e7, 2e23])
    ax2.set_xlabel('$\\nu$ [Hz]', fontsize=15 * scale_mult)
    ax2.set_ylabel('Relative Error', fontsize=15 * scale_mult)
    # ax.set_title('Plot Title', fontsize=18*scale_mult)

    ax2.tick_params(axis='both', which='major', size=12 * scale_mult, labelsize=12 * scale_mult)
    ax2.tick_params(axis='both', which='minor', size=0 * scale_mult)

    ax2.legend(loc='best', fontsize=12 * scale_mult, title_fontsize=12)

    # fig2.savefig("error_plot.pdf", dpi=200, bbox_inches="tight")

    plt.show()


def ic_bb_get_error(num_g,num_f,nu_bounds):
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
    scale_mult = 3
    nu_s, j_ic, j_ic_para = run_test(num_g=num_g, num_f=num_f)

    fig, ax = plt.subplots(figsize=(16, 12))

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.plot(nu_s, j_ic_para, label='PARAMO', linewidth=2 * scale_mult, color='red')
    ax.plot(nu_s, j_ic, label='Dermer', linewidth=2 * scale_mult, linestyle='--', color='blue')

    ax.set_xlim([5e7, 2e23])
    ax.set_ylim([1e-26, 1e-4])
    ax.set_xlabel('$\\nu$ [Hz]', fontsize=15 * scale_mult)
    ax.set_ylabel('$j_\\nu$ [$\\frac{egs}{s \ cm^3}$]', fontsize=15 * scale_mult)
    # ax.set_title('Plot Title', fontsize=18*scale_mult)

    ax.tick_params(axis='both', which='major', size=12 * scale_mult, labelsize=12 * scale_mult)
    ax.tick_params(axis='both', which='minor', size=0 * scale_mult)

    ax.legend(loc='upper left', fontsize=12 * scale_mult, title_fontsize=12)

    # fig.savefig("test.pdf", dpi=200, bbox_inches="tight")

    error = np.abs((j_ic_para - j_ic) / j_ic)

    fig2, ax2 = plt.subplots(figsize=(16, 12))

    ax2.set_xscale("log")
    ax2.set_yscale("log")

    ax2.plot(nu_s, error, label='$\\frac{\\text{|PARAMO - Dermer|}}{\\text{|Dermer|}}$', linewidth=2 * scale_mult, linestyle=':', color='green')

    ax2.set_xlim([5e7, 2e23])
    ax2.set_xlabel('$\\nu$ [Hz]', fontsize=15 * scale_mult)
    ax2.set_ylabel('$Relative Error', fontsize=15 * scale_mult)
    # ax.set_title('Plot Title', fontsize=18*scale_mult)

    ax2.tick_params(axis='both', which='major', size=12 * scale_mult, labelsize=12 * scale_mult)
    ax2.tick_params(axis='both', which='minor', size=0 * scale_mult)

    ax2.legend(loc='best', fontsize=12 * scale_mult, title_fontsize=12)

    # fig2.savefig("error_plot.pdf", dpi=200, bbox_inches="tight")

    plt.show()

if __name__ == '__main__':
    plot_comparison_and_error()
    run_convergence_test([150,300,500],[150,300,500],[5e8,5e22])