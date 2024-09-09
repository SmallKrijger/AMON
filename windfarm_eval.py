"""Windfarm layout blackbox and constraints

This script allows the user to use the blackbox function returning the Expected Annual Output and the values 
of the placing and spacing constraints with its own solver. 
The instance file must be created following the format used during the tests. The user must provide a starting 
point with the format used during the tests as well. 

This script requires libraries that are written in the "requirements.txt" to be installed in your Python environnement. 
Make sure to install them properly with the right version.

"""

import ast
import blackbox as bb
import constraints as cst
import data as d
import sys
import time
import windfarm_setting as wf

def memoize(func):
    """Script to store in cache the site and terrain of the current optimization.

    Parameters
    ----------
    func :
        The function with outputs to store in cache.

    Returns
    -------
    wrapper :
        The output of `func`.
    """

    cache = {}
    # Inner wrapper function to store the data in the cache
    def wrapper(*args):
        if args in cache:
            return cache[args]
        else:
            result = func(*args)
            cache[args] = result
            return result
    return wrapper

@memoize
def windfarm_eval(param_file_path):
    """Script that creates the site and terrain for the current optimization.

    Parameters
    ----------
    param_file_path : str
        The instance parameter file path.

    Returns
    -------
    fmGROSS : All2AllIterative
        Site with an associated windrose, turbine, wake model, blockage model, superposition model and turbulence model.
    WS_BB : Dataframe
        Dataframe for the wind speed data csv.
    WD_BB : Dataframe
        Dataframe for the wind directiond data csv.
    D : float
        Diameter of the wind turbine (used for spacing).
    buildable_zone : MultiPolygon
        Zone created by removing the exclusion zones from the boundary zone when they are overlapping.
    """

    # Initializing site and boundary files
    nb_wt, D, hub_height, scale_factor, power_curve, boundary_file, exclusion_zone_file, wind_speed, wind_direction = d.read_param_file(param_file_path)
    fmGROSS, WS, WD, max_index, wd_max = wf.site_setting(power_curve, D, hub_height, wind_speed, wind_direction)
    WS_BB, WD_BB = d.read_csv_wind_data("data/wind_speed_1.csv", "data/wind_direction_1.csv") ## remove it at the end (speeding process)
    lb, ub, boundary_shapely, exclusion_zones_shapely = wf.terrain_setting(boundary_file, exclusion_zone_file, scale_factor=scale_factor)
    buildable_zone = cst.buildable_zone(boundary_shapely, exclusion_zones_shapely)
    return fmGROSS, WS_BB, WD_BB, D, buildable_zone

## BB function
def aep(param_file_path, x): 
    """Script that compute the EAP and the constraints values, return them and store them into the 'results\bb_result.txt' file with the time of execution.

    Parameters
    ----------
    param_file_path : str
        The instance parameter file path.
    x : str or list
        The set of wind turbines coordinates.

    Returns
    -------
    eap : float
        The EAP computed with the blackbox from py_wake library.
    s_d : float
        The value of the spacing constraint (separation of the wind turbines).
    sum_dist : float
        The value of the placing constraint (distance of the wind turbines from the terrain).
    """

    t0 = time.time()
    # Building site and wind rose
    fmGROSS, WS_BB, WD_BB, D, buildable_zone = windfarm_eval(param_file_path)

    # Getting the windturbines coordinates 
    if not isinstance(x, list):
        with open(x, "r") as file:
            content = file.read().splitlines() 
            X0 = ast.literal_eval(content[0])
            x_coords = X0[0::2]
            y_coords = X0[1::2]   
    else:
        x_coords = x[0::2]
        y_coords = x[1::2]

    # Calculating EAP and constraints
    try:
        cg, eap = bb.aep_func(x_coords, y_coords, fmGROSS, WS_BB, WD_BB)
        s_d = cst.spacing_constraint_min(x_coords, y_coords, D)
        sum_dist = cst.placing_constraint(x_coords, y_coords, buildable_zone)
        eap = float(eap)

    except ValueError:
        print("Unexpected eval error", sys.exc_info()[0])
        s_d = 1e6
        sum_dist = 1e6
        eap = 0

    t1 = time.time()
    t_EAP = t1 - t0

    # Storing results
    file = open('results\bb_result.txt', 'a')
    file.write(str(eap) + " " + str(s_d) + " "  + str(sum_dist) + " " + str(t_EAP) + "\n")
    file.close()

    return eap, s_d, sum_dist

if __name__ == '__main__':
    param_path = sys.argv[1]     # Get param file pat
    x0_path = sys.argv[2]        # Get initial solution (optional)
    eap, s_d, sum_dist = windfarm_eval(param_path, x0_path)
    print("EAP = " + str(eap) + " GWh,", "Spacing constraint = " + str(s_d) + " m,", "Placing constraint = " + str(sum_dist) + " m.")