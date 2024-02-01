import pyparamo
import numpy as np
pi = pyparamo.get_pi()
print(pi)
g = np.logspace(0,5)
betav = pyparamo.bofg(g)
betas = pyparamo.bofg(g[27])
print(betav)
print(betas)
