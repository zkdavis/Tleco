import paramo as pp
import numpy as np

pi = pp.get_pi()
print(pi)
g = np.logspace(0,np.log10(1.5e8),300)
betav = pp.bofg(g)
betas = pp.bofg(g[27])
print(betav)
print(betas)

n = pp.eq_59_park1995(5e-3, g)

import matplotlib.pyplot as plt

fig,ax = plt.subplots()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([10, 1e4])
ax.set_ylim([1e-8, 1e-1])
ax.plot(g, n)
plt.show()
