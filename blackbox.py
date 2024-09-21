"""

This script contains the function for the actual blackbox.

"""

def aep_func(x, y, fmGROSS, WS_BB, WD_BB):
    """Script to evaluate the blackbox function.

    Parameters
    ----------
    x : list
        List of x coordinates. 
    y : list
        List of y coordinates. 
    fmGROSS : All2AllIterative
        Site with an associated windrose, turbine, wake model, blockage model, superposition model and turbulence model.
    WS_BB : Dataframe
        DataFrame object without the column index for the wind speed data csv.
    WD_BB : Dataframe
        DataFrame object without the column index for the wind speed data csv.

    Returns
    -------
    cg :
        fmGROSS site with associated coordinates of the wind turbines and selected wind values.
    eap :
        Eap value.
    """

    cg = fmGROSS(x, y, ws=WS_BB, wd=WD_BB, time=True, n_cpu=None)
    eap = cg.aep().sum()
    return cg, eap