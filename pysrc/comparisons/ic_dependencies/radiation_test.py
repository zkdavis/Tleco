import numpy as np
import typing
from pysrc import constants as C



def clamp(n, min, max):
    if (n < min):
        return min
    elif (n > max):
        return max
    else:
        return n

def clamp_arr(n, min, max):
    for i in range(len(n)):
        if (n[i] < min):
            n[i] = min
        else:
            n[i] = max

    return n




# Compton distribution function: see Dermer 6.75
def Fc(q, Gam):
        fc = (2 * q * np.log(q)) + ((1 + (2 * q)) * (1 - q)) + ((((Gam * q)**2) * (1 - q)) / (2 * (1 + (Gam * q))))
        return fc

def g_min_f(eps_s,eps):
    return (eps_s/2)*(1 + np.sqrt(1 + (1/(eps*eps_s))))

def g_max_f(eps_s,eps,g2):
    x = eps-eps_s
    f = np.zeros_like(x)
    for i in range(len(x)):
        if(x[i]>0):
            f[i]= eps[i]*eps_s/(eps[i] - eps_s)
        else:
           f[i] = g2
    return f

##eq 6.89 dermer
##inverse compton emission from isotropic power law electrons and an isotropic photon field
def j_ic_iso(eps_s:float,eps_array,u:typing.Callable,n_array:[float],g_array:[float]):
    fcon = (3/4)*C.cLight*C.sigmaT*eps_s

    def g_int(Fc_array,g_min_i,g_max_i):
        #assumes g_array is logspaced
        return np.trapz(n_array[g_min_i:g_max_i]*Fc_array[g_min_i:g_max_i]*np.log(10)/g_array[g_min_i:g_max_i],np.log10(g_array[g_min_i:g_max_i]))


    g_int_array = np.zeros_like(eps_array)
    g_max_array = g_max_f(eps_s,eps_array,g_array[-1])
    g_min_array = g_min_f(eps_s,eps_array)
    for i in range(len(eps_array)):
        eps_i = eps_array[i]
        Gamma_i_array = 4 * g_array * eps_i
        q_i_array = eps_s / (g_array * Gamma_i_array * (1 - (eps_s / g_array)))
        for j in range(len(q_i_array)):
            q_i_array[j] = clamp(q_i_array[j],(4*g_array[j])**(-1),1)
        Fc_i_array = Fc(q=q_i_array,Gam=Gamma_i_array)
        gmin_i = np.argmin(np.abs(g_array - g_min_array[i]))
        gmax_i = np.argmin(np.abs(g_array - g_max_array[i]))
        g_int_array[i] = g_int(Fc_i_array,gmin_i,gmax_i)
    # assumes eps_array is logspaced
    return fcon*np.trapz(u(eps_array)*g_int_array*np.log(10)/eps_array,np.log10(eps_array))

def j_ic_iso_full_dermer(eps_s_array:[float],eps_array:[float],u:typing.Callable,n_array:[float],g_array:[float]):
    j_ic = np.zeros_like(eps_s_array)
    for i in range(len(eps_s_array)):
        eps_s_i = eps_s_array[i]
        j_ic[i] = j_ic_iso(eps_s_i,eps_array,u,n_array,g_array)

    return j_ic





#2.58 Blumenthal Gould 1970
def BB_photon_density(eps_array, Temp):
    fcon = ((np.pi**2)*((C.hPlanck/(2*np.pi))**3)*(C.cLight**3))**(-1)
    f = (eps_array**2)*((np.exp(eps_array/(C.kBoltz*Temp)) - 1)**(-1))
    return fcon*f

##eq 6.66 from Dermer
def j_ic_mono_pwl_electron_dermer(eps_s:float,eps_in,u0,p,n_array:[float],g_array:[float]):
    A = np.max([g_array[0],0.5*np.sqrt(eps_s/eps_in)])
    B = np.min([g_array[-1],0.5*(1/eps_in)])
    CC = np.max([g_array[0],0.5*eps_s,eps_in*0.5])
    a = (2/(p+1))*(np.power(A,-p-1) - np.power(B,-p-1))
    b = (1/(eps_in*(p+2)))*(np.power(CC,-p-2) - np.power(g_array[-1],-p-2))
    c = -(eps_s/(2*eps_in*(p+3)))*(np.power(A,-p-3) - np.power(B,-p-3) + np.power(CC,-p-3) - np.power(g_array[-1],-p-3))
    ue = C.energy_e* np.trapz(n_array*g_array)
    ke = (p-2)*ue/(C.energy_e * (np.power(g_array[0],2-p) - np.power(g_array[-1],2-p)))
    f = ke*C.cLight*C.sigmaT*0.25*np.power(eps_s/eps_in,2)*(a + b + c)/eps_s
    return f


def j_ic_mono_pwl_electron_dermer_full(eps_s_array:[float],eps_in,u0,p,n_array:[float],g_array:[float]):
    j_ic = np.zeros_like(eps_s_array)
    for i in range(len(eps_s_array)):
        eps_s_i = eps_s_array[i]
        j_ic[i] = j_ic_mono_pwl_electron_dermer(eps_s_i,eps_in,u0,p,n_array,g_array)

    return j_ic