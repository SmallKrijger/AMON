from py_wake.superposition_models import SquaredSum
from py_wake.site import XRSite
from py_wake.wind_turbines import WindTurbine
from py_wake.wind_turbines.power_ct_functions import PowerCtTabular
from py_wake.turbulence_models import CrespoHernandez
from py_wake.deficit_models.deficit_model import WakeDeficitModel, BlockageDeficitModel
from py_wake.deficit_models import BastankhahGaussianDeficit
from py_wake.deficit_models.no_wake import NoWakeDeficit
from py_wake.wind_farm_models import All2AllIterative
from shapely.geometry import Polygon, MultiPolygon

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import shapefile
import xarray as xr

def wind_turbine_setting(powercurve_path, diameter, hub_height):
    power_curve = pd.read_csv(powercurve_path, sep=";", skiprows=0)
    u = power_curve.WindSpeed.values
    power = power_curve.Power.values*1e3
    ct = power_curve.Ct.values
    windturbine = WindTurbine(name='windturbine',
                        diameter=diameter,
                        hub_height=hub_height,
                        powerCtFunction=PowerCtTabular(u, power,'kW', ct))
    return windturbine

def site_setting(powercurve_path, diameter, hub_height, WS_path, WD_path):
    
    ## Creating Turbine
    Turbine = wind_turbine_setting(powercurve_path, diameter, hub_height)

    ## Reading files
    WS = pd.read_csv(WS_path, index_col=0)
    WD = pd.read_csv(WD_path, index_col=0)
    
    ## Creating dataframe for wind rose values
    list_of_tuples = [[None for i in range(36)] for j in range(41)]
    wd_values = np.array([i*10 for i in range(36)])
    ws_values = np.array([0] + [0.5+i for i in range(41)])
    df = pd.DataFrame(list_of_tuples, columns=wd_values, index=ws_values[1:])
    
    N = len(WS)
    width = 360 / len(wd_values)
    for i in range(len(wd_values)):
        wd = wd_values[i]
        for j in range(len(ws_values)-1):
            lower, upper = ws_values[j], ws_values[j+1]
            if wd == 0:
                sector = (360 - 0.5*width <= WD.values) | (WD.values < 0.5*width)
            else:
                sector = (wd - 0.5*width <= WD.values) & (WD.values < wd + 0.5*width)
            TS_sector = WS.values[sector]
            df.iloc[j,i] = sum((lower <= TS_sector) & (TS_sector < upper)) / N

    ## Create site object
    site = XRSite( 
        ds =xr.Dataset(
            data_vars={"P":(("wd", "ws"), df.values.T), "TI":.1},
            coords={"wd":wd_values, "ws":ws_values[1:]}
            )
        )

    wd_tot_per_ws = []
    for i in range(len(wd_values)):
        wd_tot_per_ws.append(site.ds.P.values[i].sum())

    max_index = np.argmax(wd_tot_per_ws)

    ## Creating WindRose and saving it
    script_dir = os.path.dirname(__file__)
    results_dir = os.path.join(script_dir, 'data/')

    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    fig, ax = plt.subplots()
    site.plot_wd_distribution(n_wd=36, ax=ax)
    fig.set_size_inches(3,3)
    plt.savefig(results_dir + "/WindRose.png", dpi=130)
    plt.close()

    model = BastankhahGaussianDeficit(use_effective_ws=True) 
    blockage_deficitModel = [None, model][isinstance(model, BlockageDeficitModel)]
    wake_deficitModel = [NoWakeDeficit(), model][isinstance(model, WakeDeficitModel)]
    fmGROSS = All2AllIterative(site, Turbine, wake_deficitModel=wake_deficitModel, blockage_deficitModel=blockage_deficitModel, superpositionModel=SquaredSum(), turbulenceModel=CrespoHernandez())
    return fmGROSS, WS, WD, max_index, wd_tot_per_ws[max_index]

def terrain_setting(boundary_file, constraints_file, scale_factor=0.1):
    
    # Read domain boundaries, convert to shapely format
    Boundaries = shapefile.Reader(boundary_file)
    boundary_shapely = []

    for shape in Boundaries.shapes():
        coords = np.array(shape.points).T*scale_factor
        boundary_shapely.append(Polygon(coords.T))

    boundary_shapely = [MultiPolygon(boundary_shapely)]
    
    exclusion_zones_shapely = []
    if constraints_file != "na":
        # Read exclusion zones boundaries, convert to shapely format
        Constraints = shapefile.Reader(constraints_file)

        for shape in Constraints.shapes():
            coords = np.array(shape.points).T*scale_factor
            exclusion_zones_shapely.append(Polygon(coords.T))

    lb = [(Boundaries.bbox[0])*scale_factor, (Boundaries.bbox[1])*scale_factor]
    ub = [(Boundaries.bbox[2])*scale_factor, (Boundaries.bbox[3])*scale_factor]

    return lb, ub, boundary_shapely, exclusion_zones_shapely