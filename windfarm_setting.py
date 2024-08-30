"""

This script contains the functions to create and properly set the site, terrain, windrose and turbine for the 
blackbox.

This script requires multiple libraries that are written in the "requirements.txt" to be installed in your Python environnement. 
Make sure to install them properly with the right version.

"""

from py_wake.deficit_models import BastankhahGaussianDeficit
from py_wake.deficit_models.deficit_model import WakeDeficitModel, BlockageDeficitModel
from py_wake.deficit_models.no_wake import NoWakeDeficit
from py_wake.site import XRSite
from py_wake.superposition_models import SquaredSum
from py_wake.turbulence_models import CrespoHernandez
from py_wake.wind_farm_models import All2AllIterative
from py_wake.wind_turbines import WindTurbine
from py_wake.wind_turbines.power_ct_functions import PowerCtTabular
from shapely.geometry import Polygon, MultiPolygon

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import shapefile
import xarray as xr

def wind_turbine_setting(powercurve_path, diameter, hub_height):
    """Script to create a wind turbine with the right characteristics.

    Args:
        powercurve_path (str): path to the power curve data
        diameter (float): the diameter of the wind turbine
        hub_height (float): the hub height of the wind turbine

    Results:
        windturbine (WindTurbine): an object of class WindTurbine with the right characteristics
    """

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
    """Script to create the site and get dominant wind from the wind rose.

    Args:
        powercurve_path (str): path to the power curve data
        diameter (float): the diameter of the wind turbine
        hub_height (float): the hub height of the wind turbine
        WS_path (): path to the wind speed data csv
        WD_path (): path to the wind direction data csv

    Results:
        fmGROSS (All2AllIterative): site with an associated windrose, turbine, wake model, blockage model, superposition model and turbulence model
        WS (): panda object for the wind speed data csv
        WD (): panda object for the wind directiond data csv
        max_index (): index of the wind direction with the highest frequency in the wd_tot_per_ws list
        wd_tot_per_ws[max_index] (): value of the wind direction with the highest frequency
    """

    ## Creating Turbine
    Turbine = wind_turbine_setting(powercurve_path, diameter, hub_height)

    ## Reading csv files
    WS = pd.read_csv(WS_path, index_col=0)
    WD = pd.read_csv(WD_path, index_col=0)
    
    ## Creating dataframe for wind rose values
    list_of_tuples = [[None for i in range(36)] for j in range(41)]
    wd_values = np.array([i*10 for i in range(36)])
    ws_values = np.array([0] + [0.5+i for i in range(41)])
    df = pd.DataFrame(list_of_tuples, columns=wd_values, index=ws_values[1:])
    
    ## Going through csv data and sorting them for the wind rose
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

    ## Getting the predominant wind direction frequency
    wd_tot_per_ws = []
    for i in range(len(wd_values)):
        wd_tot_per_ws.append(site.ds.P.values[i].sum())

    max_index = np.argmax(wd_tot_per_ws)

    ## Creating wind rose and saving it
    script_dir = os.path.dirname(__file__)
    results_dir = os.path.join(script_dir, 'data/')

    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    fig, ax = plt.subplots()
    site.plot_wd_distribution(n_wd=36, ax=ax)
    fig.set_size_inches(3,3)
    plt.savefig(results_dir + "/WindRose.png", dpi=130)
    plt.close()

    ## Model for wake, blockage, deficit, superposition and turbulence
    model = BastankhahGaussianDeficit(use_effective_ws=True) 
    blockage_deficitModel = [None, model][isinstance(model, BlockageDeficitModel)]
    wake_deficitModel = [NoWakeDeficit(), model][isinstance(model, WakeDeficitModel)]
    fmGROSS = All2AllIterative(site, Turbine, wake_deficitModel=wake_deficitModel, blockage_deficitModel=blockage_deficitModel, superpositionModel=SquaredSum(), turbulenceModel=CrespoHernandez())
    print(model, blockage_deficitModel, wake_deficitModel)
    print(max_index, wd_tot_per_ws[max_index])
    return fmGROSS, WS, WD, max_index, wd_tot_per_ws[max_index]

def terrain_setting(boundary_file, constraints_file, scale_factor=0.1):
    """Script to create the terrain and get the lower and upper bound of its boundary.

    Args:
        boundary_file (str): path to the boundary shapefile
        constraints_file (str): path to the exclusion zones shapefile (optional, can be set to "na" if none)
        scale_factor (float): scale factor to control the size of the terrain (usually set to 1)

    Results:
        lb (list): lower bound of the boundary zone on the x and y axis
        ub (list): upper bound of the boundary zone on the x and y axis
        boundary_shapely (list): list of type shapely Multipolygon defining the boundary zone
        exclusion_zones_shapely (list): list of type shapely Polygon defining the exclusion zones if any
    """

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