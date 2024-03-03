import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import numpy as np

def plot_n(g,n,t,min_x=None,max_x=None,min_y=None,max_y=None):

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')




    cmap = cm.rainbow
    sm = plt.cm.ScalarMappable(cmap=cmap,  norm=matplotlib.colors.LogNorm(vmin=t[0], vmax=t[-1]))

    pls = []
    my = None

    for i in range(len(t)):
        if(i%5 !=0 and i != len(t) - 1 ):
            continue
        y = n[i,:]
        if(my == None):
            my = max(y)
        if(max(y)> my):
            my= max(y)
        pl = ax.plot(g,y,color=cmap(i))
        pls.append(pl)
    maxy = 2*my
    miny = my/1e10
    if(max_y is not None):
        maxy = max_y
    if(min_y is not None):
        miny = min_y
    ax.set_ylim(miny, maxy)
    minx = None
    maxx = None
    if(min_x is not None):
        minx = min_x
        maxx = minx*(2e6)
    if(max_x is not None):
        maxx = max_x
        if(minx is None):
            minx = maxx/(2e6)
    if(minx is not None and maxx is  not None):
        ax.set_xlim(minx, maxx)

    cbticks = []
    for i in range(8):
        cbticks.append(r"$10^{{{0}}}$".format(i))
    cbar = plt.colorbar(sm,ax=ax,anchor=(-0.6,0.0),ticks=np.logspace(0,7,8))
    cbar.ax.set_yticklabels(cbticks)
    cbar.ax.minorticks_off()
    cbar.ax.set_ylabel(r"t [s]",fontsize=18)
    cbar.ax.tick_params(
        labelsize=15
    )

    plt.tick_params(
        axis='x',
        which='minor',
        bottom = False,
        top = False,
        labelbottom = False
    )

    plt.tick_params(
        axis='both',
        which='both',
        labelsize=15
    )
    plt.tick_params(
        axis='both',
        which='major',
        size=10
    )
    plt.tick_params(
        axis='y',
        which='minor',
        size=5
    )

    ax.set_xlabel(r"$\gamma$",fontsize=18)
    ax.set_ylabel(r"n $[cm^{-3}]$",fontsize=18)

    plt.tight_layout()

    # plt.savefig(plots_folder+"n1vsg.png")
    plt.show()


def plot_n_test(g,ns,ns2,t):

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')




    cmap = cm.rainbow
    sm = plt.cm.ScalarMappable(cmap=cmap,  norm=matplotlib.colors.LogNorm(vmin=t[0], vmax=t[-1]))

    pls = []
    my = None
    color_list = []
    label_list = []
    for i in range(len(ns)):
        n = ns[i]
        n2 = ns2[i]
        y = n[-1,:]
        y2 = n2[-1,:]
        y = y/np.sum(y)
        y2 = y2/np.sum(y2)
        if(my == None):
            my = max(y)
        if(max(y)> my):
            my= max(y)

        pl = ax.plot(g,y,linewidth=5.5)
        pl2 = ax.plot(g,y2,linewidth=5.5)
        pls.append(pl)
        pls.append(pl2)

    ax.set_ylim(my/1e4, 2*my)
    ax.set_xlim(1e0, 2e5)



    plt.tick_params(
        axis='x',
        which='minor',
        bottom = False,
        top = False,
        labelbottom = False
    )

    plt.tick_params(
        axis='both',
        which='both',
        labelsize=25
    )
    plt.tick_params(
        axis='both',
        which='major',
        size=20
    )
    plt.tick_params(
        axis='y',
        which='minor',
        size=0
    )

    ax.set_xlabel(r"$\gamma$",fontsize=28)
    ax.set_ylabel(r"n $[cm^{-3}]$",fontsize=28)

    plt.tight_layout()
    ax.legend()

    # plt.savefig(plots_folder+"n1vsg.png")
    plt.show()

def plot_j(g,n,t):

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')




    cmap = cm.rainbow
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=matplotlib.colors.LogNorm(vmin=t[0], vmax=t[-1]))

    pls = []
    my = None

    for i in range(len(t)):
        if(i%15 !=0 and i != len(t) - 1 ):
            continue
        y = n[i,:]
        if(my == None):
            my = max(y)
        if(max(y)> my):
            my= max(y)
        pl = ax.plot(g,y,color=cmap(i))
        pls.append(pl)

    ax.set_ylim(my/1e12, 2*my)
    # ax.set_xlim(1e0, 2e6)

    cbticks = []
    for i in range(8):
        cbticks.append(r"$10^{{{0}}}$".format(i))
    cbar = plt.colorbar(sm, ax=ax,anchor=(-0.6,0.0),ticks=np.logspace(0,7,8))
    cbar.ax.set_yticklabels(cbticks)
    cbar.ax.minorticks_off()
    cbar.ax.set_ylabel(r"t [s]",fontsize=18)
    cbar.ax.tick_params(
        labelsize=15
    )

    plt.tick_params(
        axis='x',
        which='minor',
        bottom = False,
        top = False,
        labelbottom = False
    )

    plt.tick_params(
        axis='both',
        which='both',
        labelsize=15
    )
    plt.tick_params(
        axis='both',
        which='major',
        size=10
    )
    plt.tick_params(
        axis='y',
        which='minor',
        size=5
    )

    ax.set_xlabel(r"$\nu",fontsize=18)
    ax.set_ylabel(r"$\nu j$ [ Hz erg $s^{-1}$ $cm^{-3}$]",fontsize=18)

    plt.tight_layout()

    # plt.savefig(plots_folder+"n1vsg.png")
    plt.show()

def plot_I(g,n,t):

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')




    cmap = cm.rainbow
    sm = plt.cm.ScalarMappable(cmap=cmap,  norm=matplotlib.colors.LogNorm(vmin=t[0], vmax=t[-1]))

    pls = []
    my = None

    for i in range(len(t)):
        if(i%15 !=0 and i != len(t) - 1 ):
            continue
        y = n[i,:]
        if(my == None):
            my = max(y)
        if(max(y)> my):
            my= max(y)
        pl = ax.plot(g,y,color=cmap(i))
        pls.append(pl)

    ax.set_ylim(my/1e8, 2*my)
    # ax.set_xlim(1e0, 2e6)

    cbticks = []
    for i in range(8):
        cbticks.append(r"$10^{{{0}}}$".format(i))
    cbar = plt.colorbar(sm,ax=ax, anchor=(-0.6,0.0),ticks=np.logspace(0,7,8))
    cbar.ax.set_yticklabels(cbticks)
    cbar.ax.minorticks_off()
    cbar.ax.set_ylabel(r"t [s]",fontsize=18)
    cbar.ax.tick_params(
        labelsize=15
    )

    plt.tick_params(
        axis='x',
        which='minor',
        bottom = False,
        top = False,
        labelbottom = False
    )

    plt.tick_params(
        axis='both',
        which='both',
        labelsize=15
    )
    plt.tick_params(
        axis='both',
        which='major',
        size=10
    )
    plt.tick_params(
        axis='y',
        which='minor',
        size=5
    )

    ax.set_xlabel(r"$\nu",fontsize=18)
    ax.set_ylabel(r"$\nu 4 \pi I$ [ Hz erg $s^{-1}$ $cm^{-2}$]",fontsize=18)

    plt.tight_layout()

    # plt.savefig(plots_folder+"n1vsg.png")
    plt.show()