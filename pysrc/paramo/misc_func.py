import numpy as np

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


def broken_pwl(n0, g, p1, p2, gmin_cut, g2_cut):
    f = np.zeros(len(g))
    mask1 = (g > gmin_cut) & (g <= g2_cut)
    f[mask1] = g[mask1] ** p1
    if np.any(g > g2_cut):
        i0 = np.argmax(g > g2_cut) - 1
        mask2 = g > g2_cut
        f[mask2] = f[i0] * (g[mask2] / g[i0]) ** p2
    f = n0 * f / np.trapz(f,g)

    return f


def power_law(n0, g, p, g_min, g_max,normalize=False):
    """
    @func: Computes a power law function.
    @param n0: Normalization constant.
    @param g: Array of input values.
    @param p: Power law index.
    @param g_min: Minimum value of Lorentz factor for the power law application.
    @param g_max: Maximum value of Lorentz factor for the power law application.
    @param normalize (default=False): bool if true normalizes the distribution before returning
    @return f: Array of output values according to the power law.
    """
    f = np.zeros_like(g)
    bounds = (g >= g_min) & (g <= g_max)
    f[bounds] = np.power(g, p)[bounds]
    if(normalize):
        f = f/np.trapz(f,g)
    return f*n0


