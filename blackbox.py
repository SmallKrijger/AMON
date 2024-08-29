# Expected annual production
def aep_func(x, y, fmGROSS, WS_BB, WD_BB):
    cg = fmGROSS(x, y, ws=WS_BB, wd=WD_BB, time=True)
    eap = cg.aep().sum()
    return cg, eap