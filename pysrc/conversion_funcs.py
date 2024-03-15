def ev2ergs(ev):
    ergs = ev*1.6021766339999e-12
    return ergs

def erg2ev(ergs):
    ev = ergs/(1.6021766339999e-12)
    return ev

def ev2hz(ev):
    hz = 2.417989242*(1e14)*ev
    return hz

def hz2ev(hz):
    ev = hz/(2.417989242*(1e14))
    return ev