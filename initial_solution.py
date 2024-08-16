import itertools
import sys
import time
from matplotlib import pyplot as plt
import numpy as np
import windfarm as wf
import constraints as cst
import data as d
import PyNomad

def initial_solution_constrained(nb_wt, D, buildable_zone, lb, ub):
    print("Generating initial solution...")
    x, y = initial_sol_nomad(nb_wt, buildable_zone, D, lb, ub, spacing=True)
    X0 = sum(map(list, zip(x, y)), [])
    X0 = list(itertools.chain(*zip(x, y)))
    return X0

def initial_sol_nomad(nb_wt, buildable_zone, D, lb, ub, spacing=True):

    def initial_sol(x):
        try:
            # Get coords
            dim = x.size()
            x_list = [x.get_coord(i) for i in range(dim)]
            x_coords = x_list[0::2]
            y_coords = x_list[1::2]

            # Calculate constraints
            sum_dist = cst.placing_constraint(x_coords, y_coords, buildable_zone)

            if spacing:
                s_d = cst.spacing_constraint_min(x_coords, y_coords, D)
                rawBBO = str(0) + " " + str(s_d) + " " + str(sum_dist)
                x.setBBO(rawBBO.encode("UTF-8"))
            else:
                rawBBO = str(0) + " " + str(sum_dist)
                x.setBBO(rawBBO.encode("UTF-8"))

        except:
            print("Unexpected eval error", sys.exc_info()[0])
            return 0
        return 1
    
    params = ["BB_INPUT_TYPE * R", "STOP_IF_FEASIBLE true", "DISPLAY_DEGREE 0", "DISPLAY_ALL_EVAL false", "DIRECTION_TYPE ORTHO 2N"] 

    if spacing:
        params += ["BB_OUTPUT_TYPE OBJ PB PB"]
    else:
        params += ["BB_OUTPUT_TYPE OBJ PB"]

    for i in range(0, nb_wt, 2):
        params += [f"VARIABLE_GROUP {i} {i+1}"]

    x = np.random.uniform(lb[0], ub[0], nb_wt)
    y = np.random.uniform(lb[1], ub[1], nb_wt)
    X0 = sum(map(list, zip(x, y)), [])
    X0 = list(itertools.chain(*zip(x, y)))

    ## NOMAD optimization 
    result = PyNomad.optimize(initial_sol, X0, lb*nb_wt, ub*nb_wt, params)

    obj_function_value = result['f_best']
    
    x_best = result['x_best'][0::2]
    y_best = result['x_best'][1::2]
    return x_best, y_best

# x, y = initial_sol_nomad(5000, 30)
# print(time_Nomad)

def initial_sol_test(param_file_name):
     # Initializing site and boundary files
    nb_wt, D, hub_height, scale_factor, power_curve, boundary_file, exclusion_zone_file, wind_speed, wind_direction = d.read_param_file(param_file_name)
    lb, ub, boundary_shapely, exclusion_zones_shapely = wf.spatial_constraints(boundary_file, exclusion_zone_file, scale_factor=scale_factor)
    buildable_zone = cst.buildable_zone(boundary_shapely, exclusion_zones_shapely)
    x, y = initial_sol_nomad(nb_wt, buildable_zone, D, lb, ub, spacing=True)
    X0 = sum(map(list, zip(x, y)), [])
    X0 = list(itertools.chain(*zip(x, y)))
    f = open("tests/5/x0.txt", 'w+')  # open file in write mode
    f.write(str(X0))
    f.close()   

initial_sol_test("tests/5/param.txt")
            