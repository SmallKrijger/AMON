import time
import windfarm as wf
import constraints as cst
import PyNomad
import matplotlib.pyplot as plt
import data as d
import ast
import sys

def memoize(func):
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
def windfarm_opt(param_file_path):
    # Initializing site and boundary files
    nb_wt, D, hub_height, scale_factor, power_curve, boundary_file, exclusion_zone_file, wind_speed, wind_direction = d.read_param_file(param_file_path)
    fmGROSS, WS, WD, max_index, wd_max = wf.site_model(power_curve, D, hub_height, wind_speed, wind_direction)
    WS_BB, WD_BB = wf.read_csv("data/wind_speed_1.csv", "data/wind_direction_1.csv") ## remove it at the end (speeding process)
    lb, ub, boundary_shapely, exclusion_zones_shapely = wf.spatial_constraints(boundary_file, exclusion_zone_file, scale_factor=scale_factor)
    buildable_zone = cst.buildable_zone(boundary_shapely, exclusion_zones_shapely)
    return fmGROSS, WS_BB, WD_BB, D, buildable_zone

## BB function
def aep(param_file_path, x):   
    t0 = time.time()
    # Building site and wind rose
    fmGROSS, WS_BB, WD_BB, D, buildable_zone = windfarm_opt(param_file_path)

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
        cg, eap = wf.aep_func(x_coords, y_coords, fmGROSS, WS_BB, WD_BB)
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
    x0_path = sys.argv[2]    # Get initial solution (optional)
    eap, s_d, sum_dist = windfarm_opt(param_path, x0_path)
    print("EAP = " + str(eap) + " GWh,", "Spacing constraint = " + str(s_d) + " m,", "Placing constraint = " + str(sum_dist) + " m.")