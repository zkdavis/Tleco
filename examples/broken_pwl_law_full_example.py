import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors
import paramo as para
from paramo import misc_func as mf
from paramo import constants as cons
##Example injecting a broken power law into a blob that is cooled by synchrotron and ssc


# Parameters
R = 5e15
B = 0.1
n0 = 2
p1, p2 = -2.5, -3.5
gcut, g2cut = 1e2, 1e5
tlc = R / cons.cLight
tmax = tlc * 4

# Initialization
numt, numg, numf = 30, 60, 80
fmin, fmax = 1e8, 1e28
gmin, gmax = 1e0, 1e8
t = np.logspace(0, np.log10(tmax), numt)
g = np.logspace(np.log10(gmin), np.log10(gmax), numg)
f = np.logspace(np.log10(fmin), np.log10(fmax), numf)
n = np.zeros([numt, numg])
gdot = np.zeros([numt, numg])
D = np.full(numg, 1e-200)
Qinj = np.zeros(numg)
j_s = np.zeros([numt, numf])
j_ssc = np.zeros([numt, numf])
I_s = np.zeros([numt, numf])
I_ssc = np.zeros([numt, numf])
ambs = np.zeros([numt, numf])

#initial conditions
gdot[0, :] = (4 / 3) * cons.sigmaT * cons.cLight * (B**2 / (8 * np.pi)) * g**2 / (cons.me * cons.cLight**2)
n[0, :] = mf.broken_pwl(n0, g, p1, p2, gcut, g2cut)

# time loop
for i in range(1, len(t)):
    dt = t[i] - t[i-1]
    n[i, :] = para.fp_findif_difu(dt, g, n[i-1, :], gdot[i-1, :], D, Qinj, 1e200, tlc)
    j_s[i, :], ambs[i, :] = para.syn_emissivity_full(f, g, n[i, :], B, True)
    I_s[i, :] = para.rad_trans_blob(R, j_s[i, :], ambs[i, :])
    j_ssc[i, :] = para.ic_iso_powlaw_full(f, I_s[i, :], g, n[i, :])
    I_ssc[i, :] = para.rad_trans_blob(R, j_ssc[i, :], ambs[i, :])
    dotgKN = para.rad_cool_pwl(g, f, 4 * np.pi * I_s[i, :] / cons.cLight, True)
    gdot[i, :] = gdot[0, :] + dotgKN

def plot_generic(x, y, t, x_label, y_label, scale='log', cbar_label=r"t [s]", plot_every=1, filename=None, min_x=None, max_x=None, min_y=None, max_y=None):
    fig, ax = plt.subplots(figsize=(10,8))
    ax.set_xscale(scale)
    ax.set_yscale(scale)

    cmap = cm.rainbow
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=matplotlib.colors.LogNorm(vmin=t[0], vmax=t[-1]))
    for i in range(len(t)):
        if i % plot_every != 0 and i != len(t) - 1:
            continue
        ax.plot(x, y[i, :], color=cmap(i / len(t)))

    if min_x is not None and max_x is not None:
        ax.set_xlim(min_x, max_x)
    if min_y is not None and max_y is not None:
        ax.set_ylim(min_y, max_y)

    cbar = plt.colorbar(sm, ax=ax)
    cbar.ax.set_ylabel(cbar_label, fontsize=18)
    ax.set_xlabel(x_label, fontsize=18)
    ax.set_ylabel(y_label, fontsize=18)
    plt.tight_layout()
    if filename:
        plt.savefig(filename)
    plt.show()


plot_generic(g, n, t, r"$\gamma$", r"n $[cm^{-3}]$", plot_every=5, min_x=1e1, max_x=1e5, min_y=1e-6, max_y=5e-2)
plot_generic(f, (I_ssc + I_s) * f, t, r"$\nu$ [Hz]", r"$\nu I$ [ Hz erg $s^{-1}$ $cm^{-2}$]", plot_every=5, min_x=1e9, max_x=1e29, min_y=1e-4, max_y=1e4)



