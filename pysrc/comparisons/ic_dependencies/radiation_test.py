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

def u_black_body(eps,sig,norm):
    f1 = 8*np.pi*C.me*(C.cLight**2)*(eps**3)
    lc = C.hPlanck/(C.me*C.cLight)
    temexp = np.exp(eps/sig)
    for i in range(len(temexp)):
        if temexp[i] <= 1.0000001:
            temexp[i] = 1e200
    f = f1/((lc**3)*(temexp- 1 ))
    f = norm*f/np.trapz(f,eps)
    return f


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

# def j_ic_iso_full(eps_s_array:[float],eps_array:[float],u:typing.Callable,n_array:[float],g_array:[float]):
#     j_ic = np.zeros_like(eps_s_array)
#     for i in range(len(eps_s_array)):
#         eps_s_i = eps_s_array[i]
#         j_ic[i] = j_ic_iso(eps_s_i,eps_array,u,n_array,g_array)
#
#     return j_ic

# def j_ic_iso_BB(eps_s_array:[float],eps_array:[float],n_array:[float],g_array:[float],sig_T:float,u_bb:float):
#     def ub(eps):
#         return u_black_body(eps,sig_T,u_bb)
#     j_ic = j_ic_iso_full(eps_s_array,eps_array,ub,n_array,g_array)
#     return j_ic


def F_q(q,s,p):
    num= (2*q*np.log(q)) + (1 + q) - (2 * q**2) + (2*s*q*(1-q))
    den = ((1 + np.sqrt((s*q)/(1 + s*q)))**(p+2))*((1 + s*q)**((p+3)/2))
    return num/den

#2.58 Blumenthal Gould 1970
def BB_photon_density(eps_array, Temp):
    fcon = ((np.pi**2)*((C.hPlanck/(2*np.pi))**3)*(C.cLight**3))**(-1)
    f = (eps_array**2)*((np.exp(eps_array/(C.kBoltz*Temp)) - 1)**(-1))
    return fcon*f

#2.75 Blumenthal Gould 1970
def j_ic_iso_pwl_BB(eps_s:float,eps_array:[float],u0,Temp,p,n0,g_array):
    r0 = np.sqrt(C.sigmaT*3/(8*np.pi))
    fcon = np.pi*(r0**2)*C.cLight*n0*(2**(p+1))*(eps_s**(-(p+1)/2))*u0

    def q_int(q_array,s):
        q_array = np.flip(q_array)
        q_min_i = np.argmin(np.abs(q_array-0))
        q_max_i = np.argmin(np.abs(q_array-1))
        fq_array = F_q(q_array,s,p)
        integrand = (q**((p-1)/2))*fq_array
        return np.trapz(integrand[q_min_i:q_max_i],q_array[q_min_i:q_max_i])

    q_int_array = np.zeros_like(eps_array)
    N_bb = BB_photon_density(eps_array,Temp)

    for i in range(len(eps_array)):
        eps = eps_array[i]
        s = eps * eps_s / (C.energy_e ** 2)
        Gamma_e = 4 * eps * g_array / C.energy_e
        E1 = eps_s/(g_array*C.energy_e)
        q = E1 / ( Gamma_e * (1 - E1))
        q_int_array[i] = q_int(q,s)

    return fcon*np.trapz((eps_array**((p-1)/2))*N_bb*q_int_array,eps_array)



def j_ic_iso_full(eps_s_array:[float],eps_array:[float],u0:float,temp:float, p ,n0,g_array:[float]):
    j_ic = np.zeros_like(eps_s_array)
    for i in range(len(eps_s_array)):
        eps_s_i = eps_s_array[i]
        j_ic[i] = j_ic_iso_pwl_BB(eps_s_i,eps_array,u0,temp, p,n0,g_array)

    return j_ic