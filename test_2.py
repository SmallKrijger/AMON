# load a time series of wd, ws and ti
# Install PyWake if needed
import py_wake

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os

# import and setup site and windTurbines
import data as d
import constraints as cst
import windfarm_setting as wf
import ast
import initial_solution as init_sol
import blackbox as bb
import sys 
import PyNomad
import time


def NOMAD_execution3(param_file_name, i):
    """Script to execute the NOMAD solver on the instances.

    Parameters
    ----------
    param_file_name : str
        Path to the instance parameter file.
    x0 : str 
        Path to the initial set of coordinates.

    Returns
    -------
    instances_results/.../convergence_plot.png :
        Convergence plot.
    instances_results/.../layout.png :
        Plot of the best layout found for the wind farm.
    instances_results/.../nomad_result.0.txt :
        Text file with the results of all evaluations of NOMAD.
    """

    # Initializing site and boundary files
    instance_number = param_file_name.split('/')[1]
    script_dir = os.path.dirname(__file__)
    results_instance_dir = os.path.join(script_dir, 'instances_results/' + instance_number + "/")

    if not os.path.isdir(results_instance_dir):
        os.makedirs(results_instance_dir)

    nb_wt, D, hub_height, scale_factor, power_curve, boundary_file, exclusion_zone_file, wind_speed, wind_direction = d.read_param_file(param_file_name)
    params, nb_it = d.read_config_file("data/config.txt", nb_wt)
    for j in range(0, nb_wt, 2):
        params.append(f"VARIABLE_GROUP {j} {j+1}")
    params.append(f"SEED {i}")
    fmGROSS, WS, WD, max_index, max_ws = wf.site_setting(power_curve, D, hub_height, wind_speed, wind_direction, results_instance_dir)
    WS_BB, WD_BB = d.read_csv_wind_data(wind_speed, wind_direction)
    lb, ub, boundary_shapely, exclusion_zones_shapely = wf.terrain_setting(boundary_file, exclusion_zone_file, scale_factor=scale_factor)
    buildable_zone = cst.buildable_zone(boundary_shapely, exclusion_zones_shapely)

    ## BB function
    def eap_nomad(x):
        """Script to compute the eap and constraint values with the NOMAD solver.

        Parameters
        ----------
        x : list
            List of (x,y) coordinates.

        Returns
        -------
        boolean :
            1 if the evaluation was successful, 0 otherwise.
        """

        try:
            # Get coords
            dim = x.size()
            x_list = [x.get_coord(i) for i in range(dim)]
            x_coords = x_list[0::2]
            y_coords = x_list[1::2]

            if cst.checking_same_coords(x_coords, y_coords):
                # Calculate constraints
                sum_dist = cst.placing_constraint(x_coords, y_coords, buildable_zone)

                # Calculate EAP
                if sum_dist == 0:
                    cg, eap, wl = bb.aep_func(x_coords, y_coords, fmGROSS, WS_BB, WD_BB)
                    s_d = cst.spacing_constraint_min(x_coords, y_coords, D)
                
                    # NOMAD output
                    eap = -float(eap)*1000000  
                    rawBBO = str(eap) + " " + str(s_d)
                    x.setBBO(rawBBO.encode("UTF-8"))
                else:
                    rawBBO = "-"
                    x.setBBO(rawBBO.encode("UTF-8"))
                
            else:
                rawBBO = "-"
                x.setBBO(rawBBO.encode("UTF-8"))

        except:
            print("Unexpected eval error", sys.exc_info()[0])
            return 0
        
        return 1
    
    plt.figure()
    
    stats_file_name = results_instance_dir + "/nomad_result_" + str(nb_wt) + ".0.txt"
    nomad_stat_file = "STATS_FILE  "+ stats_file_name + "  BBE OBJ CONS_H SOL TIME"
    params += [nomad_stat_file]

    # Generating an valid initial solution
    print("--------- Running instance", instance_number, "---------")
    with open("instances/" + instance_number + "/x0.txt", "r") as file:
        content = file.read().splitlines() 
    X0 = ast.literal_eval(content[0])

    ## NOMAD optimization 
    print("...Launching NOMAD optimization...")
    t2 = time.time()
    result = PyNomad.optimize(eap_nomad, X0, nb_wt*lb, nb_wt*ub, params)
    t3 = time.time()
    t_Nomad = t3 - t2
    print("...Ending NOMAD optimization...")

    instance_number = param_file_name.split('/')[1] 
    # np_evals, np_obj, best_eval, best_of = d.read_stat_file(nb_it, stat_file_name="instances_results/" + instance_number + "/nomad_result_" + instance_number + ".0.txt")
    # plt.clf()
    # plot_f.plot_result_nomad(np_evals, np_obj, best_eval, best_of, "instances_results/" + instance_number + "/convergence_plot_" + instance_number + ".png" )

    obj_function_value = -result['f_best']*10**(-6)
    cg, eap, wl = bb.aep_func(result['x_best'][0::2], result['x_best'][1::2], fmGROSS, WS_BB, WD_BB)
    # plot_f.plot_terrain(result['x_best'][0::2], result['x_best'][1::2], "EAP", "GWh", obj_function_value, nb_wt, ub, boundary_shapely, exclusion_zones_shapely, max_index=max_index, max_ws=max_ws, cg=cg, plot_flow_map=False, save=True, save_name="instances_results/" + instance_number + "/layout_" + instance_number + ".png")
    print("Best objective function value : ", obj_function_value, "GWh")
    print("NOMAD execution time : ", t_Nomad, " seconds")
    print("Wake loss : ", wl, " %")
    print("--------- Ending instance", instance_number, "---------")
    
    f = open("instances_results/output_instances.txt", 'a')
    f.write(str(i) + " " + str(obj_function_value) + " " + str(t_Nomad) + "\n")
    # f.write(str(nb_wt) + " " + str(0) + " " + str(0) + " " + str(0) + "\n")
    f.close() 

if __name__ == '__main__':
    for i in range(1, 30):
        NOMAD_execution3("instances/1/param.txt", i)