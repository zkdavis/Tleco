from astropy import constants as const
import numpy as np

halfpi = 0.5 * np.pi
twopi = 2.0 * np.pi
lzero = np.log(1e-200)
# Physical constants in CGS
cLight = const.c.cgs.value
G = const.G.cgs.value
eCharge = const.e.gauss.value
hPlanck = const.h.cgs.value
me = const.m_e.cgs.value
mp = const.m_p.cgs.value
sigmaT = const.sigma_T.cgs.value
sigmaSB = const.sigma_sb.cgs.value
kBoltz = const.k_B.cgs.value

nuConst = eCharge / (twopi * me * cLight)
PmbConst = twopi * np.sqrt(3) * eCharge**2 / cLight
jmbConst = np.sqrt(3) * eCharge**2 / (2. * cLight)
ambConst = np.sqrt(3) * eCharge**2 / (4. * me * cLight)
chunche_c100g100 = 2.2619939050180366385e-6
chunche_c100g20 = 2.1157699720918349273e-1

Msun = const.M_sun.cgs.value