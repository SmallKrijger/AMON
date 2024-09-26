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

def NOMAD_execution(config_file_name, param_file_name, i, x0):
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
    params, nb_it = d.read_config_file(config_file_name, nb_wt)
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
                # if sum_dist == 0:
                cg, eap = bb.aep_func(x_coords, y_coords, fmGROSS, WS_BB, WD_BB)
                s_d = cst.spacing_constraint_min(x_coords, y_coords, D)
                # print(sum_dist, eap, s_d)
                # NOMAD output
                eap = -float(eap)*1000000  
                rawBBO = str(eap) + " " + str(sum_dist) + " " + str(s_d)
                x.setBBO(rawBBO.encode("UTF-8"))
                # else:
                    # rawBBO = "-"
                    # x.setBBO(rawBBO.encode("UTF-8"))
                
            else:
                rawBBO = "-"
                x.setBBO(rawBBO.encode("UTF-8"))

        except:
            print("Unexpected eval error", sys.exc_info()[0])
            return 0
        
        return 1
    
    plt.figure()
    # test_number = param_file_name.split('/')[1]
    script_dir = os.path.dirname(__file__)
    results_test_dir = os.path.join(script_dir, 'tests_results/' + i + "/")

    if not os.path.isdir(results_test_dir):
        os.makedirs(results_test_dir)
    stats_file_name = results_test_dir + "/nomad_result_" + i + ".0.txt"
    nomad_stat_file = "STATS_FILE  "+ stats_file_name + "  BBE OBJ CONS_H SOL TIME"
    params += [nomad_stat_file]

    # Generating an valid initial solution
    print("--------- Beginning Test", i, "---------")
    # with open("tests/" + i + "/x0.txt", "r") as file:
    with open(x0, "r") as file:
        content = file.read().splitlines() 
    X0 = ast.literal_eval(content[0])

    ## NOMAD optimization 
    print("Launching NOMAD optimization...")
    t2 = time.time()
    result = PyNomad.optimize(eap_nomad, X0, nb_wt*lb, nb_wt*ub, params)
    t3 = time.time()
    t_Nomad = t3 - t2

    # i = param_file_name.split('/')[1] 
    np_evals, np_obj, best_eval, best_of = d.read_stat_file(nb_it, stat_file_name="tests_results/" + i + "/nomad_result_" + i + ".0.txt")
    plt.clf()
    plot_f.plot_result_nomad(np_evals, np_obj, best_eval, best_of, "tests_results/" + i + "/convergence_plot_" + i + ".png" )

    obj_function_value = -result['f_best']*10**(-6)
    cg, eap = bb.aep_func(result['x_best'][0::2], result['x_best'][1::2], fmGROSS, WS_BB, WD_BB)
    plot_f.plot_terrain(result['x_best'][0::2], result['x_best'][1::2], "EAP", "GWh", obj_function_value, nb_wt, ub, boundary_shapely, exclusion_zones_shapely, max_index=max_index, cg=cg, plot_flow_map=False, save=True, save_name="tests_results/" + i + "/layout_" + i + ".png")
    print("Best objective function value : ", obj_function_value, "GWh")
    print("NOMAD execution time : ", t_Nomad, " seconds")
    print("--------- Ending Test", i, "---------")
    
    f = open("tests_results/output_tests.txt", 'a')
    f.write(str(i) + " " + str(obj_function_value) + " " + str(t_Nomad) + "\n")
    f.close()   

def NOMAD_execution_1(config_file_name, param_file_name, i, x0):
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
    params, nb_it = d.read_config_file(config_file_name, nb_wt)
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
                    # print(sum_dist, eap, s_d)
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
    # test_number = param_file_name.split('/')[1]
    script_dir = os.path.dirname(__file__)
    results_test_dir = os.path.join(script_dir, 'tests_results/' + i + "/")

    if not os.path.isdir(results_test_dir):
        os.makedirs(results_test_dir)
    stats_file_name = results_test_dir + "/nomad_result_" + i + ".0.txt"
    nomad_stat_file = "STATS_FILE  "+ stats_file_name + "  BBE OBJ CONS_H SOL TIME"
    params += [nomad_stat_file]

    # Generating an valid initial solution
    print("--------- Beginning Test", i, "---------")
    # with open("tests/" + i + "/x0.txt", "r") as file:
    with open(x0, "r") as file:
        content = file.read().splitlines() 
    X0 = ast.literal_eval(content[0])

    ## NOMAD optimization 
    print("Launching NOMAD optimization...")
    t2 = time.time()
    result = PyNomad.optimize(eap_nomad, X0, nb_wt*lb, nb_wt*ub, params)
    t3 = time.time()
    t_Nomad = t3 - t2

    # i = param_file_name.split('/')[1] 
    np_evals, np_obj, best_eval, best_of = d.read_stat_file(nb_it, stat_file_name="tests_results/" + i + "/nomad_result_" + i + ".0.txt")
    plt.clf()
    plot_f.plot_result_nomad(np_evals, np_obj, best_eval, best_of, "tests_results/" + i + "/convergence_plot_" + i + ".png" )

    obj_function_value = -result['f_best']*10**(-6)
    cg, eap = bb.aep_func(result['x_best'][0::2], result['x_best'][1::2], fmGROSS, WS_BB, WD_BB)
    plot_f.plot_terrain(result['x_best'][0::2], result['x_best'][1::2], "EAP", "GWh", obj_function_value, nb_wt, ub, boundary_shapely, exclusion_zones_shapely, max_index=max_index, cg=cg, plot_flow_map=False, save=True, save_name="tests_results/" + i + "/layout_" + i + ".png")
    print("Best objective function value : ", obj_function_value, "GWh")
    print("NOMAD execution time : ", t_Nomad, " seconds")
    print("--------- Ending Test", i, "---------")
    
    f = open("tests_results/output_tests.txt", 'a')
    f.write(str(i) + " " + str(obj_function_value) + " " + str(t_Nomad) + "\n")
    f.close()   

def NOMAD_execution_2(config_file_name, param_file_name, i, x0):
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
    params, nb_it = d.read_config_file(config_file_name, nb_wt)
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
                # sum_dist = cst.placing_constraint(x_coords, y_coords, buildable_zone)
                l_bool = wf.test_constraints_placing(boundary_shapely, exclusion_zones_shapely, x_coords, y_coords)
                sum_l_bool = sum(l_bool)
                n_wrong = nb_wt - sum_l_bool

                # Calculate EAP
                # if sum_dist == 0:
                cg, eap = bb.aep_func(x_coords, y_coords, fmGROSS, WS_BB, WD_BB)
                # s_d = cst.spacing_constraint_min(x_coords, y_coords, D)
                # print(sum_dist, eap, s_d)
                # NOMAD output
                eap = -float(eap)*1000000  
                rawBBO = str(eap) + " " + str(n_wrong)
                x.setBBO(rawBBO.encode("UTF-8"))
                # else:
                    # rawBBO = "-"
                    # x.setBBO(rawBBO.encode("UTF-8"))
                
            else:
                rawBBO = "-"
                x.setBBO(rawBBO.encode("UTF-8"))

        except:
            print("Unexpected eval error", sys.exc_info()[0])
            return 0
        
        return 1
    
    plt.figure()
    # test_number = param_file_name.split('/')[1]
    script_dir = os.path.dirname(__file__)
    results_test_dir = os.path.join(script_dir, 'tests_results/' + i + "/")

    if not os.path.isdir(results_test_dir):
        os.makedirs(results_test_dir)
    stats_file_name = results_test_dir + "/nomad_result_" + i + ".0.txt"
    nomad_stat_file = "STATS_FILE  "+ stats_file_name + "  BBE OBJ CONS_H SOL TIME"
    params += [nomad_stat_file]

    # Generating an valid initial solution
    print("--------- Beginning Test", i, "---------")
    # with open("tests/" + i + "/x0.txt", "r") as file:
    with open(x0, "r") as file:
        content = file.read().splitlines() 
    X0 = ast.literal_eval(content[0])

    ## NOMAD optimization 
    print("Launching NOMAD optimization...")
    t2 = time.time()
    result = PyNomad.optimize(eap_nomad, X0, nb_wt*lb, nb_wt*ub, params)
    t3 = time.time()
    t_Nomad = t3 - t2

    # i = param_file_name.split('/')[1] 
    np_evals, np_obj, best_eval, best_of = d.read_stat_file(nb_it, stat_file_name="tests_results/" + i + "/nomad_result_" + i + ".0.txt")
    plt.clf()
    plot_f.plot_result_nomad(np_evals, np_obj, best_eval, best_of, "tests_results/" + i + "/convergence_plot_" + i + ".png" )

    obj_function_value = -result['f_best']*10**(-6)
    cg, eap = bb.aep_func(result['x_best'][0::2], result['x_best'][1::2], fmGROSS, WS_BB, WD_BB)
    plot_f.plot_terrain(result['x_best'][0::2], result['x_best'][1::2], "EAP", "GWh", obj_function_value, nb_wt, ub, boundary_shapely, exclusion_zones_shapely, max_index=max_index, cg=cg, plot_flow_map=False, save=True, save_name="tests_results/" + i + "/layout_" + i + ".png")
    print("Best objective function value : ", obj_function_value, "GWh")
    print("NOMAD execution time : ", t_Nomad, " seconds")
    print("--------- Ending Test", i, "---------")
    
    f = open("tests_results/output_tests.txt", 'a')
    f.write(str(i) + " " + str(obj_function_value) + " " + str(t_Nomad) + "\n")
    f.close()   

def NOMAD_execution_3(config_file_name, param_file_name, i, x0):
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
    params, nb_it = d.read_config_file(config_file_name, nb_wt)
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
                # sum_dist = cst.placing_constraint(x_coords, y_coords, buildable_zone)
                l_bool = wf.test_constraints_placing(boundary_shapely, exclusion_zones_shapely, x_coords, y_coords)
                sum_l_bool = sum(l_bool)
                # n_wrong = nb_wt - sum_l_bool

                # Calculate EAP
                if sum_l_bool == nb_wt:
                    cg, eap = bb.aep_func(x_coords, y_coords, fmGROSS, WS_BB, WD_BB)
                    # s_d = cst.spacing_constraint_min(x_coords, y_coords, D)
                    # print(sum_dist, eap, s_d)
                    # NOMAD output
                    eap = -float(eap)*1000000  
                    rawBBO = str(eap)
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
    # test_number = param_file_name.split('/')[1]
    script_dir = os.path.dirname(__file__)
    results_test_dir = os.path.join(script_dir, 'tests_results/' + i + "/")

    if not os.path.isdir(results_test_dir):
        os.makedirs(results_test_dir)
    stats_file_name = results_test_dir + "/nomad_result_" + i + ".0.txt"
    nomad_stat_file = "STATS_FILE  "+ stats_file_name + "  BBE OBJ CONS_H SOL TIME"
    params += [nomad_stat_file]

    # Generating an valid initial solution
    print("--------- Beginning Test", i, "---------")
    # with open("tests/" + i + "/x0.txt", "r") as file:
    with open(x0, "r") as file:
        content = file.read().splitlines() 
    X0 = ast.literal_eval(content[0])

    ## NOMAD optimization 
    print("Launching NOMAD optimization...")
    t2 = time.time()
    result = PyNomad.optimize(eap_nomad, X0, nb_wt*lb, nb_wt*ub, params)
    t3 = time.time()
    t_Nomad = t3 - t2

    # i = param_file_name.split('/')[1] 
    np_evals, np_obj, best_eval, best_of = d.read_stat_file(nb_it, stat_file_name="tests_results/" + i + "/nomad_result_" + i + ".0.txt")
    plt.clf()
    plot_f.plot_result_nomad(np_evals, np_obj, best_eval, best_of, "tests_results/" + i + "/convergence_plot_" + i + ".png" )

    obj_function_value = -result['f_best']*10**(-6)
    cg, eap = bb.aep_func(result['x_best'][0::2], result['x_best'][1::2], fmGROSS, WS_BB, WD_BB)
    plot_f.plot_terrain(result['x_best'][0::2], result['x_best'][1::2], "EAP", "GWh", obj_function_value, nb_wt, ub, boundary_shapely, exclusion_zones_shapely, max_index=max_index, cg=cg, plot_flow_map=False, save=True, save_name="tests_results/" + i + "/layout_" + i + ".png")
    print("Best objective function value : ", obj_function_value, "GWh")
    print("NOMAD execution time : ", t_Nomad, " seconds")
    print("--------- Ending Test", i, "---------")
    
    f = open("tests_results/output_tests.txt", 'a')
    f.write(str(i) + " " + str(obj_function_value) + " " + str(t_Nomad) + "\n")
    f.close()   

def NOMAD_execution_4(config_file_name, param_file_name, i, x0):
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
    params, nb_it = d.read_config_file(config_file_name, nb_wt)
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
                # if sum_dist == 0:
                cg, eap = bb.aep_func(x_coords, y_coords, fmGROSS, WS_BB, WD_BB)
                # s_d = cst.spacing_constraint_min(x_coords, y_coords, D)
                # print(sum_dist, eap, s_d)
                # NOMAD output
                eap = -float(eap)*1000000  
                rawBBO = str(eap) + " " + str(sum_dist)
                x.setBBO(rawBBO.encode("UTF-8"))
                # else:
                #     rawBBO = "-"
                #     x.setBBO(rawBBO.encode("UTF-8"))
                
            else:
                rawBBO = "-"
                x.setBBO(rawBBO.encode("UTF-8"))

        except:
            print("Unexpected eval error", sys.exc_info()[0])
            return 0
        
        return 1
    
    plt.figure()
    # test_number = param_file_name.split('/')[1]
    script_dir = os.path.dirname(__file__)
    results_test_dir = os.path.join(script_dir, 'tests_results/' + i + "/")

    if not os.path.isdir(results_test_dir):
        os.makedirs(results_test_dir)
    stats_file_name = results_test_dir + "/nomad_result_" + i + ".0.txt"
    nomad_stat_file = "STATS_FILE  "+ stats_file_name + "  BBE OBJ CONS_H SOL TIME"
    params += [nomad_stat_file]

    # Generating an valid initial solution
    print("--------- Beginning Test", i, "---------")
    # with open("tests/" + i + "/x0.txt", "r") as file:
    with open(x0, "r") as file:
        content = file.read().splitlines() 
    X0 = ast.literal_eval(content[0])

    ## NOMAD optimization 
    print("Launching NOMAD optimization...")
    t2 = time.time()
    result = PyNomad.optimize(eap_nomad, X0, nb_wt*lb, nb_wt*ub, params)
    t3 = time.time()
    t_Nomad = t3 - t2

    # i = param_file_name.split('/')[1] 
    np_evals, np_obj, best_eval, best_of = d.read_stat_file(nb_it, stat_file_name="tests_results/" + i + "/nomad_result_" + i + ".0.txt")
    plt.clf()
    plot_f.plot_result_nomad(np_evals, np_obj, best_eval, best_of, "tests_results/" + i + "/convergence_plot_" + i + ".png" )

    obj_function_value = -result['f_best']*10**(-6)
    cg, eap = bb.aep_func(result['x_best'][0::2], result['x_best'][1::2], fmGROSS, WS_BB, WD_BB)
    plot_f.plot_terrain(result['x_best'][0::2], result['x_best'][1::2], "EAP", "GWh", obj_function_value, nb_wt, ub, boundary_shapely, exclusion_zones_shapely, max_index=max_index, cg=cg, plot_flow_map=False, save=True, save_name="tests_results/" + i + "/layout_" + i + ".png")
    print("Best objective function value : ", obj_function_value, "GWh")
    print("NOMAD execution time : ", t_Nomad, " seconds")
    print("--------- Ending Test", i, "---------")
    
    f = open("tests_results/output_tests.txt", 'a')
    f.write(str(i) + " " + str(obj_function_value) + " " + str(t_Nomad) + "\n")
    f.close()

def NOMAD_execution_5(config_file_name, param_file_name, i, x0):
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
    params, nb_it = d.read_config_file(config_file_name, nb_wt)
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
                    # s_d = cst.spacing_constraint_min(x_coords, y_coords, D)
                    # print(sum_dist, eap, s_d)
                    # NOMAD output
                    eap = -float(eap)*1000000  
                    rawBBO = str(eap)
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
    # test_number = param_file_name.split('/')[1]
    script_dir = os.path.dirname(__file__)
    results_test_dir = os.path.join(script_dir, 'tests_results/' + i + "/")

    if not os.path.isdir(results_test_dir):
        os.makedirs(results_test_dir)
    stats_file_name = results_test_dir + "/nomad_result_" + i + ".0.txt"
    nomad_stat_file = "STATS_FILE  "+ stats_file_name + "  BBE OBJ CONS_H SOL TIME"
    params += [nomad_stat_file]

    # Generating an valid initial solution
    print("--------- Beginning Test", i, "---------")
    # with open("tests/" + i + "/x0.txt", "r") as file:
    with open(x0, "r") as file:
        content = file.read().splitlines() 
    X0 = ast.literal_eval(content[0])

    ## NOMAD optimization 
    print("Launching NOMAD optimization...")
    t2 = time.time()
    result = PyNomad.optimize(eap_nomad, X0, nb_wt*lb, nb_wt*ub, params)
    t3 = time.time()
    t_Nomad = t3 - t2

    # i = param_file_name.split('/')[1] 
    np_evals, np_obj, best_eval, best_of = d.read_stat_file(nb_it, stat_file_name="tests_results/" + i + "/nomad_result_" + i + ".0.txt")
    plt.clf()
    plot_f.plot_result_nomad(np_evals, np_obj, best_eval, best_of, "tests_results/" + i + "/convergence_plot_" + i + ".png" )

    obj_function_value = -result['f_best']*10**(-6)
    cg, eap = bb.aep_func(result['x_best'][0::2], result['x_best'][1::2], fmGROSS, WS_BB, WD_BB)
    plot_f.plot_terrain(result['x_best'][0::2], result['x_best'][1::2], "EAP", "GWh", obj_function_value, nb_wt, ub, boundary_shapely, exclusion_zones_shapely, max_index=max_index, cg=cg, plot_flow_map=False, save=True, save_name="tests_results/" + i + "/layout_" + i + ".png")
    print("Best objective function value : ", obj_function_value, "GWh")
    print("NOMAD execution time : ", t_Nomad, " seconds")
    print("--------- Ending Test", i, "---------")
    
    f = open("tests_results/output_tests.txt", 'a')
    f.write(str(i) + " " + str(obj_function_value) + " " + str(t_Nomad) + "\n")
    f.close()


def testing_process():
    """Script to launch the testing process.

        Parameters
        ----------
            
        Returns
        -------
        NOMAD_execution function results.
        """
    # if os.path.exists("tests_results/output_tests.txt"):
    #     os.remove("tests_results/output_tests.txt")
    param_txt = ["tests/1/param.txt", "tests/4/param.txt"]
    x0_txt = ["tests/1/x0.txt", "tests/4/x0.txt"]
    test_failed = []
    for j in range(0,2):
        # for i in range(0,4):
        #     try:
        #         NOMAD_execution(config_file_name="data/config_" + str(i) + ".txt", param_file_name=param_txt[j], i=str(i), x0=x0_txt[j])
        #     except:
        #         test_failed.append(i)
        # for i in range(4,6):
        #     try:
        #         NOMAD_execution_1(config_file_name="data/config_" + str(i) + ".txt", param_file_name=param_txt[j], i=str(i), x0=x0_txt[j])
        #     except:
        #         test_failed.append(i)
        for i in range(6,8):
            try:
                NOMAD_execution_2(config_file_name="data/config_" + str(i) + ".txt", param_file_name=param_txt[j], i=str(i), x0=x0_txt[j])
            except:
                test_failed.append(i)
        for i in range(8,9):
            try:
                NOMAD_execution_3(config_file_name="data/config_" + str(i) + ".txt", param_file_name=param_txt[j], i=str(i), x0=x0_txt[j])
            except:
                test_failed.append(i)
        for i in range(9,11):
            try:
                NOMAD_execution_4(config_file_name="data/config_" + str(i) + ".txt", param_file_name=param_txt[j], i=str(i), x0=x0_txt[j])
            except:
                test_failed.append(i)
        for i in range(11,12):
            try:
                NOMAD_execution_5(config_file_name="data/config_" + str(i) + ".txt", param_file_name=param_txt[j], i=str(i), x0=x0_txt[j])
            except:
                test_failed.append(i)
    # for i in range(12,13):
    #     try:
    #         NOMAD_execution_5(config_file_name="data/config_" + str(i) + ".txt", param_file_name="tests/1/param.txt", i=str(i))
    #     except:
    #         test_failed.append(i)
    if test_failed != []:
        print("Test failed: ", test_failed)
    print("--------- Tests completed ---------")


if __name__ == '__main__':
    testing_process()