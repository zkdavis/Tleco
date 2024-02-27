import misc_func
import paramo as para
import numpy as np
import paramo.plots_code as pc
import constants as cons

##Example injectecting a broken power law into a blob that is cooled by synchrotron and ssc
###constants
numt = 800
numg = 80
numf = 80
fmax = 1e28
fmin =1e8
tmax = 1e18
gmin = 1e0
gmax = 1e8
with_abs = True
cool_withKN = True

R = 3.2e15 # size of blob
B = 0.05 #blobs magnetic field
uB = (B**2)/(np.pi*8) #magnetic field energy density
C0 = 3.48e-11 #4 * sigmaT * uB / (3 * mass_e * cLight)
t_acc = 1 / (C0 * (1e4)) #acceleration time scale
t_esc = t_acc #escape time
t_inj = R / cons.cLight #injection time
tlc = R / cons.cLight #light crossing time
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
D_0 = 0.5 * np.power(g,2)/t_acc
D = 2 * D_0 # diffusion array
gdot[0,:] = C0 * np.power(g,2) - 4*D_0/g - g/t_acc #cooling array
Qinj = np.zeros(numg) #injection distribution


p1 = -2.5 #pwl indices
gcut = 1e0
g2cut = 2e0
n[0,:] = misc_func.power_law(1, g, p1, gcut, g2cut) #initial distribution


###time loop
for i in range(1,len(t)):
    dt = t[i] - t[i-1]
    n[i,:] = para.fp_findif_difu(dt, g, n[i-1,:], gdot[i-1,:], D, Qinj, 1e200, tlc)
    print(i)
    j_s[i,:],ambs[i,:] = para.syn_emissivity_full(f, g,n[i,:], B, with_abs) #,sync and absorb
    #I_s[i,:] = para.rad_trans_blob(R, j_s[i,:], ambs[i,:])
    #j_ssc[i,:] = para.ic_iso_powlaw_full(f, I_s[i,:], g, n[i,:])
    # j_eic[i,:] = para.ic_iso_powlaw_full(f,I_s[i,:],g,n[i,:])
    #I_ssc[i, :] = para.rad_trans_blob(R, j_ssc[i,:], ambs[i,:])
   # dotgKN = para.rad_cool_pwl(g, f, 4 * np.pi * I_ssc[i,:] / cons.cLight, cool_withKN)
    gdot[i,:] = gdot[0,:]# + dotgKN
    print(i)

pc.plot_n(g, n, t)
# pc.plot_j(f,f*(j_s),t)
pc.plot_j(f, f * (j_s + j_ssc + j_eic), t)
# pc.plot_I(f,np.pi * 4* (I_s)*f,t)
pc.plot_I(f, np.pi * 4 * (I_ssc + I_s) * f, t)
# pc.plot_n(g,gdot[0:1,:],t[0:1])
# pc.plot_n(g,gdot,t)
