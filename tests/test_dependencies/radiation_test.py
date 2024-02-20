import numpy as np
import typing
import constants as C



def u_black_body(eps,sig):
    f1 = 8*np.pi*C.me*(C.cLight**2)*(eps**3)
    lc = C.hPlanck/(C.me*C.cLight)
    f = f1/((lc**3)*(np.exp(eps/sig) - 1 ))
    return f


# Compton distribution function: see Dermer 6.75
def Fc(e_s=None, g=None, e_i=None, q=None, Gam=None):
    # have to specify args
    # g is lorents factor, e_s is reduced scattered photon energy and e_i is reduced incoming photon energy.
    if q is None:
        Gam = 4 * g * e_i
        q = e_s / (g * Gam * (1 - (e_s / g)))
        if(q>1 or q<(1/(4*(g**2)))):
            return 0
        fc = 2 * q * np.log(q) + (1 + (2 * q)) * (1 - q) + (((Gam * q)**2) * (1 - q)) / (2 * (1 + (Gam * q)))

        return fc
    else:
        fc = 2 * q * np.log(q) + (1 + (2 * q)) * (1 - q) + (((Gam * q)**2) * (1 - q)) / (2 * (1 + (Gam * q)))
        return fc



##eq 6.89 dermer
##inverse compton emission from isotropic electrons and an isotropic photon field
def j_ic_iso(eps_s:float,eps_array,u:typing.Callable,n_array:[float],g_array:[float]):
    fcon = (3/4)*C.cLight*C.sigmaT*eps_s

    def g_int(Fc_array):
        #assumes g_array is logspaced
        return np.trapz(n_array*Fc_array*np.log(10)/g_array,g_array)


    g_int_array = np.zeros_like(eps_array)
    for i in range(len(eps_array)):
        eps_i = eps_array[i]
        Gamma_i_array = 4 * g_array * eps_i
        q_i_array = eps_s / (g_array * Gamma_i_array * (1 - (eps_s / g_array)))
        Fc_i_array = Fc(q=q_i_array,Gam=Gamma_i_array)
        g_int_array[i] = g_int(Fc_i_array)
    # assumes eps_array is logspaced
    return fcon*np.trapz(u(eps_array)*g_int_array*np.log(10)/eps_array)

def j_ic_iso_full(eps_s_array:[float],eps_array:[float],u:typing.Callable,n_array:[float],g_array:[float]):
    j_ic = np.zeros_like(eps_s_array)
    for i in range(len(eps_s_array)):
        eps_s_i = eps_s_array[i]
        j_ic[i] = j_ic_iso(eps_s_i,eps_array,u,n_array,g_array)

    return j_ic

def j_ic_iso_BB(eps_s_array:[float],eps_array:[float],n_array:[float],g_array:[float],sig_T:float):
    def ub(eps):
        return u_black_body(eps,sig_T)
    j_ic = j_ic_iso_full(eps_s_array,eps_array,ub,n_array,g_array)
    return j_ic

if(__name__=="__main__"):
    import matplotlib.pyplot as plt
    p = -1
    sig_T = 1
    eps_s = np.logspace(5,30)*C.hPlanck/(C.me*(C.cLight**2))
    eps = eps_s
    g_array = np.logspace(1,5)
    n_array = (g_array**p)/np.trapz(g_array**p,g_array)

    j_ic = j_ic_iso_BB(eps_s,eps,n_array,g_array,sig_T)

    fig,ax = plt.subplots()

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.plot(eps_s,j_ic)

    plt.show()