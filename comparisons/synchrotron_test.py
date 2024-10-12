import matplotlib.pyplot as plt
import tleco as tl
import numpy as np
from dependencies import radiation_test as rt

def run_test(num_g,num_f):
    p = 2.0
    n0 =1
    nu_s = np.logspace(4, 15, num_f)
    ub=1
    B=np.sqrt(np.pi*8*ub)
    g_array = np.logspace(1, 3, num_g)
    n_array = n0*(g_array ** -p) / np.trapz(g_array ** -p, g_array)
    j_syn = rt.j_syn_explicit(nu_s,B,n_array,g_array) * 4*np.pi
    j_syn_para = np.array(tl.syn_emissivity_full(nu_s,g_array,n_array,B,False))[0] * 4*np.pi

    return nu_s,j_syn, j_syn_para


def run_convergence_test(num_g_values,num_f_values,nu_bounds):
    scale_mult = 3
    errors = np.zeros((len(num_g_values), len(num_f_values)))
    errors_spec = np.zeros((len(num_g_values), len(num_f_values)),dtype=object)
    nu_array = []
    for i, num_g in enumerate(num_g_values):
        for j, num_f in enumerate(num_f_values):
            nu_s, j_syn, j_syn_para = run_test(num_g=num_g, num_f=num_f)
            nu_min_i = np.argmin(np.abs(nu_s - nu_bounds[0]))
            nu_max_i = np.argmin(np.abs(nu_s - nu_bounds[1]))
            nu_s = nu_s[nu_min_i:nu_max_i]
            j_syn = j_syn[nu_min_i:nu_max_i]
            j_syn_para = j_syn_para[nu_min_i:nu_max_i]
            non_zero_indices = j_syn != 0
            j_syn_filtered = j_syn[non_zero_indices]
            j_syn_para_filtered = j_syn_para[non_zero_indices]
            nu_array.append(nu_s[non_zero_indices])
            er_spec = np.abs((j_syn_para_filtered - j_syn_filtered)/j_syn_filtered)
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


    ax2.set_xlabel('$\\nu$ [Hz]', fontsize=15 * scale_mult)
    ax2.set_ylabel('Relative Error', fontsize=15 * scale_mult)
    # ax.set_title('Plot Title', fontsize=18*scale_mult)

    ax2.tick_params(axis='both', which='major', size=12 * scale_mult, labelsize=12 * scale_mult)
    ax2.tick_params(axis='both', which='minor', size=0 * scale_mult)

    ax2.legend(loc='best', fontsize=12 * scale_mult, title_fontsize=12)

    # fig2.savefig("error_plot.pdf", dpi=200, bbox_inches="tight")

    plt.show()


def ic_mono_get_error(num_g,num_f,nu_bounds):
    nu_s, j_syn, j_syn_para = run_test(num_g=num_g, num_f=num_f)
    nu_min_i = np.argmin(np.abs(nu_s - nu_bounds[0]))
    nu_max_i = np.argmin(np.abs(nu_s - nu_bounds[1]))
    nu_s = nu_s[nu_min_i:nu_max_i]
    j_syn = j_syn[nu_min_i:nu_max_i]
    j_syn_para = j_syn_para[nu_min_i:nu_max_i]
    non_zero_indices = j_syn != 0
    j_syn_filtered = j_syn[non_zero_indices]
    j_syn_para_filtered = j_syn_para[non_zero_indices]
    er_spec = np.abs((j_syn_para_filtered - j_syn_filtered) / j_syn_filtered)
    return np.sum(er_spec) / len(er_spec)


def plot_comparison_and_error():
    num_g = 300#150
    num_f = 300#150
    scale_mult = 4
    nu_s, j_syn, j_syn_para = run_test(num_g=num_g, num_f=num_f)
    print(np.trapz(j_syn, nu_s))
    fig, ax = plt.subplots(figsize=(16, 12))

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.plot(nu_s, nu_s*j_syn_para, label='TLECO', linewidth=3 * scale_mult, color='red')
    ax.plot(nu_s, nu_s*j_syn, label='Dermer', linewidth=3 * scale_mult, linestyle='--', color='blue')

    ax.set_xlim([2e5, 1e15])
    ax.set_ylim([2e-16, 1e-10])
    ax.set_xlabel('$\\nu$ [Hz]', fontsize=15 * scale_mult)
    ax.set_ylabel('$\\nu \ j_\\nu$ [$\\frac{ergs}{s \ cm^3}$]', fontsize=15 * scale_mult)
    # ax.set_title('Plot Title', fontsize=18*scale_mult)

    ax.tick_params(axis='both', which='major', size=12 * scale_mult, labelsize=12 * scale_mult)
    ax.tick_params(axis='both', which='minor', size=0 * scale_mult)

    ax.legend(loc='upper left', fontsize=12 * (scale_mult - 1), title_fontsize=12)

    # fig.savefig("Figs/syn_test.pdf", dpi=800, bbox_inches="tight")
    # fig.savefig("Figs/syn_test.png", dpi=800, bbox_inches="tight")

    error = np.abs((j_syn_para - j_syn) / j_syn)

    fig2, ax2 = plt.subplots(figsize=(16, 12))

    ax2.set_xscale("log")
    ax2.set_yscale("log")

    ax2.plot(nu_s, error, label='$\\frac{\\text{|TLECO - Dermer|}}{\\text{|Dermer|}}$', linewidth=2 * scale_mult, linestyle=':', color='green')

    ax2.set_xlim([2e12, 1e20])
    ax2.set_xlabel('$\\nu$ [Hz]', fontsize=15 * scale_mult)
    ax2.set_ylabel('$Relative Error', fontsize=15 * scale_mult)
    # ax.set_title('Plot Title', fontsize=18*scale_mult)

    ax2.tick_params(axis='both', which='major', size=12 * scale_mult, labelsize=12 * scale_mult)
    ax2.tick_params(axis='both', which='minor', size=0 * scale_mult)

    ax2.legend(loc='best', fontsize=12 * scale_mult, title_fontsize=12)

    # fig2.savefig("error_plot.pdf", dpi=200, bbox_inches="tight")
    fig.savefig(f"Figs/syn_comparison.pdf", dpi=800, bbox_inches="tight")
    fig.savefig(f"Figs/syn_comparison.png", dpi=800, bbox_inches="tight")

    # plt.show()

if __name__ == '__main__':
    plot_comparison_and_error()
    # run_convergence_test([25,50,150,300,500,1000],[25,50,150,300,500,1000],[2e13,2.5e19])