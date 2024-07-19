import itertools
import time
import windfarm as wf
import constraints as cst
import PyNomad
import sys  
import shapefile
import os
from shapely.geometry import Polygon, Point
from shapely import distance
import numpy as np
import matplotlib.pyplot as plt
import coords_colors as cc
import geopandas as gpd
import data as d
import initial_solution as init_sol

scale_factor = 0.1

def NOMAD_execution(param_file_name, x0=""):
    
    # Initializing site and boundary files
    nb_wt, D, boundary_file, exclusion_zone_file, wind_speed, wind_direction = d.read_param_file(param_file_name)
    params, nb_it = d.read_config_file("data/config.txt", nb_wt)
    site, fmGROSS, WS, WD, max_index, wd_max = wf.site_model(wind_speed, wind_direction, WS_column_name='0', WD_column_name='0', reduce_file_size=True)
    lb, ub, boundary_shapely, exclusion_zones_shapely = wf.spatial_constraints(boundary_file, exclusion_zone_file, scale_factor=scale_factor)
    buildable_zone = cst.buildable_zone(boundary_shapely, exclusion_zones_shapely)

    ## BB function
    def eap_nomad(x):
        try:
            # Get coords
            dim = x.size()
            x_list = [x.get_coord(i) for i in range(dim)]
            x_coords = x_list[0::2]
            y_coords = x_list[1::2]

            # Calculate constraints
            sum_dist = cst.placing_constraint(x_coords, y_coords, buildable_zone)
            s_d = cst.spacing_constraint_min(x_coords, y_coords, D)

            # Calculate EAP
            cg, eap = wf.aep_func(x_coords, y_coords, fmGROSS, WS, WD, TS_column_name='0', DS_column_name='0')
            
            # NOMAD output
            eap = -float(eap)*1000000  
            rawBBO = str(eap)  + " " + str(sum_dist) + " " + str(s_d)
            x.setBBO(rawBBO.encode("UTF-8"))
        except:
            print("Unexpected eval error", sys.exc_info()[0])
            return 0
        return 1
    
    plt.figure()
    test_number = param_file_name.split('/')[1] 
    stats_file_name = "tests_results/" + test_number + "/nomad_result_" + test_number + ".0.txt"
    nomad_stat_file = "STATS_FILE  "+ stats_file_name + "  BBE OBJ CONS_H SOL TIME"
    params += [nomad_stat_file]


    # Generating an valid initial solution
    print("--------- Beginning Test", test_number, "---------")
    print("Generating initial solution...")
    x, y = init_sol.initial_sol_nomad(nb_wt, buildable_zone, boundary_shapely, exclusion_zones_shapely, D, lb, ub, max_index, spacing=True)
    X0 = sum(map(list, zip(x, y)), [])
    X0 = list(itertools.chain(*zip(x, y)))
    plt.clf()

    ## NOMAD optimization 
    print("Launching NOMAD optimization...")
    t2 = time.time()
    result = PyNomad.optimize(eap_nomad, X0, nb_wt*lb, nb_wt*ub, params)
    t3 = time.time()
    t_Nomad = t3 - t2
    print("NOMAD execution time : ", t_Nomad, " seconds")

    test_number = param_file_name.split('/')[1] 
    np_evals, np_obj, best_eval, best_of = d.read_stat_file(nb_it, stat_file_name="tests_results/" + test_number + "/nomad_result_" + test_number + ".0.txt")
    d.draw_result_nomad(np_evals, np_obj, best_eval, best_of, nb_it, nb_wt, "tests_results/" + test_number + "/convergence_plot_" + test_number + ".png" )

    obj_function_value = -result['f_best']*10**(-6)
    print("Best objective function value : ", obj_function_value, "GWh")
    wf.plot_spatial_cstr_generation(result['x_best'][0::2], result['x_best'][1::2], "EAP", "GWh", obj_function_value, nb_wt, ub, boundary_shapely, exclusion_zones_shapely, max_index=max_index, save=True, save_name="tests_results/" + test_number + "/layout_" + test_number + ".png")
    print("--------- Ending Test", test_number, "---------")
    # return t_Nomad


def testing_process():
    for i in range(2,3):
        NOMAD_execution(param_file_name="tests/" + str(i) + "/param.txt")
# t_Nomad = plot_NOMAD()
# print(t_Nomad)

if __name__ == '__main__':
    # param_file_name = sys.argv[1]     # Get parameters file path
    testing_process()