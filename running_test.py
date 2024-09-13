"""

This script is used to run the test to validate the installation of the user. 

This script requires multiple libraries that are written in the "requirements.txt" to be installed in your Python environnement. 
Make sure to install them properly with the right version.

"""

import ast
import blackbox as bb
import constraints as cst
import data as d
import matplotlib.pyplot as plt
import os
import plotting_functions as plot_f
import PyNomad
import sys
import time
import windfarm_setting as wf

def NOMAD_execution(param_file_name, x0=""):
    """Script to execute the NOMAD solver on the tests instance.

    Parameters
    ----------
    param_file_name : str
        Path to the instance parameter file.
    x0 : str 
        Path to the initial set of coordinates.

    Returns
    -------
    tests_results\...\convergence_plot.png :
        Convergence plot.
    tests_results\...\layout.png :
        Plot of the best layout found for the wind farm.
    tests_results\...\nomad_result.0.txt :
        Text file with the results of all evaluations of NOMAD.
    """

    # Initializing site and boundary files
    nb_wt, D, hub_height, scale_factor, power_curve, boundary_file, exclusion_zone_file, wind_speed, wind_direction = d.read_param_file(param_file_name)
    params, nb_it = d.read_config_file("data/config.txt", nb_wt)
    fmGROSS, WS, WD, max_index, wd_max = wf.site_setting(power_curve, D, hub_height, wind_speed, wind_direction)
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
                    cg, eap = bb.aep_func(x_coords, y_coords, fmGROSS, WS_BB, WD_BB)
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
    test_number = param_file_name.split('/')[1]
    script_dir = os.path.dirname(__file__)
    results_test_dir = os.path.join(script_dir, 'tests_results/' + test_number + "/")

    if not os.path.isdir(results_test_dir):
        os.makedirs(results_test_dir)
    stats_file_name = results_test_dir + "/nomad_result_" + test_number + ".0.txt"
    nomad_stat_file = "STATS_FILE  "+ stats_file_name + "  BBE OBJ CONS_H SOL TIME"
    params += [nomad_stat_file]

    # Generating an valid initial solution
    print("--------- Beginning Test", test_number, "---------")
    with open("tests/" + test_number + "/x0.txt", "r") as file:
        content = file.read().splitlines() 
    X0 = ast.literal_eval(content[0])

    ## NOMAD optimization 
    print("Launching NOMAD optimization...")
    t2 = time.time()
    result = PyNomad.optimize(eap_nomad, X0, nb_wt*lb, nb_wt*ub, params)
    t3 = time.time()
    t_Nomad = t3 - t2

    test_number = param_file_name.split('/')[1] 
    np_evals, np_obj, best_eval, best_of = d.read_stat_file(nb_it, stat_file_name="tests_results/" + test_number + "/nomad_result_" + test_number + ".0.txt")
    plt.clf()
    plot_f.plot_result_nomad(np_evals, np_obj, best_eval, best_of, "tests_results/" + test_number + "/convergence_plot_" + test_number + ".png" )

    obj_function_value = -result['f_best']*10**(-6)
    cg, eap = bb.aep_func(result['x_best'][0::2], result['x_best'][1::2], fmGROSS, WS_BB, WD_BB)
    plot_f.plot_terrain(result['x_best'][0::2], result['x_best'][1::2], "EAP", "GWh", obj_function_value, nb_wt, ub, boundary_shapely, exclusion_zones_shapely, max_index=max_index, cg=cg, plot_flow_map=False, save=True, save_name="tests_results/" + test_number + "/layout_" + test_number + ".png")
    # plot_f.plot_terrain([], [], "EAP", "GWh", 0, 0, ub, boundary_shapely, exclusion_zones_shapely, max_index=max_index, cg=cg, plot_flow_map=False, save=True, save_name="tests_results/" + test_number + "/layout_" + test_number + ".png")
    print("Best objective function value : ", obj_function_value, "GWh")
    print("NOMAD execution time : ", t_Nomad, " seconds")
    print("--------- Ending Test", test_number, "---------")
    
    f = open("tests_results/output_tests.txt", 'a')
    f.write(str(test_number) + " " + str(obj_function_value) + " " + str(t_Nomad) + "\n")
    f.close()   

def testing_process():
    """Script to launch the testing process.

        Parameters
        ----------
            
        Returns
        -------
        NOMAD_execution function results.
        """
    if os.path.exists("tests_results/output_tests.txt"):
        os.remove("tests_results/output_tests.txt")

    test_failed = []
    for i in range(5,6):
        try:
            NOMAD_execution(param_file_name="tests/" + str(i) + "/param.txt")
        except:
            test_failed.append(i)
    if test_failed != []:
        print("Test failed: ", test_failed)
    print("--------- Tests completed ---------")


if __name__ == '__main__':
    testing_process()