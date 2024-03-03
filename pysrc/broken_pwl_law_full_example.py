import paramo as para
import numpy as np
import paramo.plots_code as pc
import misc_func as mf

##Example injectecting a broken power law into a blob that is cooled by synchrotron and ssc

cLight = 2.99792458e10
mp = 1.67262e-24
me = 9.1093897e-28
sigmaT = 6.6524e-25
###constants
numt = 30
numg = 60
numf = 80
fmax = 1e28
fmin =1e8
tmax = 1e18
gmin = 1e1
gmax = 1e8
with_abs = True
cool_withKN = False

R = 5e14 # size of blob
sigma = 1e-6 #blobs magnetization
B = 0.1 #blobs magnetic field
n0 = B**2 / (np.pi * 4 * sigma * mp * cLight**2)#particle density
t_esc = 1e200#R/cLight #escape time
t_inj = R / cLight #injection time
tlc = R / cLight #light crossing time
p1 = -2.5 #pwl indices
p2 =p1 -1
gcut = 1e2
g2cut = 1e5

tmax = tlc * 4

###arrays
n = np.zeros([numt, numg]) #particle distribution dn/d\gamma
gdot = np.zeros([numt, numg]) #fp cooling term
j_s = np.zeros([numt, numf]) #synchrotron emissivity
j_ssc = np.zeros([numt, numf]) #ssc emissivity
j_eic = np.zeros([numt, numf]) #ssc emissivity
I_s = np.zeros([numt, numf])#synchrotron Intensity
I_ssc = np.zeros([numt, numf])#ssc Intensity
ambs = np.zeros([numt, numf]) #absorbtion coefficient
t = np.logspace(0, np.log10(tmax), numt) # logspaced time array where t_0 = 1
f = np.logspace(np.log10(fmin), np.log10(fmax), numf) # logspaced time array where t_0 = 1
g = np.logspace(np.log10(gmin), np.log10(gmax), numg) #logspaced lorentz factor array
D = np.full(numg, 1e-200) # diffusion array
gdot[0,:] = np.full(numg, 1e-200) #cooling array
Qinj = np.zeros(numg) #injection distribution


##define fp terms
# Qinj = broken_pwl(n0, g, p1, p2, gcut, g2cut) / t_inj
gdot[0,:] = (4 / 3) * sigmaT * cLight * (B**2 / (8 * np.pi)) * g**2 / (me * cLight**2)
# D = 0.5*gdot[0,:]
n[0,:] = mf.broken_pwl(n0, g, p1, p2, gcut, g2cut)
uext = 1e-3
nuout = 1e14
###time loop
for i in range(1,len(t)):
    dt = t[i] - t[i-1]
    n[i,:] = para.fp_findif_difu(dt, g, n[i-1,:], gdot[i-1,:], D, Qinj, t_esc, tlc)
    j_s[i,:],ambs[i,:] = para.syn_emissivity_full(f, g,n[i,:], B, with_abs) #,sync and absorb
    I_s[i,:] = para.rad_trans_blob(R, j_s[i,:], ambs[i,:])
    j_ssc[i,:] = para.ic_iso_powlaw_full(f, I_s[i,:], g, n[i,:])
    # j_eic[i,:] = para.ic_iso_powlaw_full(f,I_s[i,:],g,n[i,:])
    I_ssc[i, :] = para.rad_trans_blob(R, j_ssc[i,:], ambs[i,:])
    dotgKN = para.rad_cool_pwl(g, f, 4 * np.pi * I_ssc[i,:] / cLight, cool_withKN)
    gdot[i,:] = gdot[0,:] + dotgKN

pc.plot_n(g, n, t)
# pc.plot_j(f,f*(j_s),t)
pc.plot_j(f, f * (j_s + j_ssc + j_eic), t)
# pc.plot_I(f,np.pi * 4* (I_s)*f,t)
pc.plot_I(f, np.pi * 4 * (I_ssc + I_s) * f, t)
# pc.plot_n(g,gdot[0:1,:],t[0:1])
# pc.plot_n(g,gdot,t)
