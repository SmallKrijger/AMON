import time
import windfarm as wf
import constraints as cst
import PyNomad
import matplotlib.pyplot as plt
import data as d
import ast
import sys
import os

def NOMAD_execution(param_file_name, x0=""):
    
    # Initializing site and boundary files
    nb_wt, D, hub_height, scale_factor, power_curve, boundary_file, exclusion_zone_file, wind_speed, wind_direction = d.read_param_file(param_file_name)
    params, nb_it = d.read_config_file("data/config.txt", nb_wt)
    fmGROSS, WS, WD, max_index, wd_max = wf.site_model(power_curve, D, hub_height, wind_speed, wind_direction)
    WS_BB, WD_BB = wf.read_csv(wind_speed, wind_direction)
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
            # print(x_coords, y_coords)
            if cst.checking_same_coords(x_coords, y_coords):
                # Calculate constraints
                l_bool = wf.test_constraints_placing(boundary_shapely, exclusion_zones_shapely, x_coords, y_coords)
                sum_l_bool = sum(l_bool)
                n_wrong = nb_wt - sum_l_bool
                # Calculate EAP
                # if sum_l_bool == nb_wt:
                    
                cg, eap = wf.aep_func(x_coords, y_coords, fmGROSS, WS_BB, WD_BB)
                # s_d = cst.spacing_constraint_min(x_coords, y_coords, D)
            
                # NOMAD output
                eap = -float(eap)*1000000  
                rawBBO = str(eap) + " " + str(n_wrong)
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
    d.draw_result_nomad(np_evals, np_obj, best_eval, best_of, nb_it, nb_wt, "tests_results/" + test_number + "/convergence_plot_" + test_number + ".png" )

    obj_function_value = -result['f_best']*10**(-6)
    cg, eap = wf.aep_func(result['x_best'][0::2], result['x_best'][1::2], fmGROSS, WS_BB, WD_BB)
    wf.plot_spatial_cstr_generation(result['x_best'][0::2], result['x_best'][1::2], "EAP", "GWh", obj_function_value, nb_wt, ub, boundary_shapely, exclusion_zones_shapely, max_index=max_index, cg=cg, plot_flow_map=False, save=True, save_name="tests_results/" + test_number + "/layout_" + test_number + ".png")
    print("Best objective function value : ", obj_function_value, "GWh")
    print("NOMAD execution time : ", t_Nomad, " seconds")
    print("--------- Ending Test", test_number, "---------")
    f = open("tests_results/output_tests.txt", 'a')  # open file in write mode
    f.write(str(test_number) + " " + str(obj_function_value) + " " + str(t_Nomad) + "\n")
    f.close()   

def testing_process():
    test_failed = []
    for i in range(1,2):
        try:
            NOMAD_execution(param_file_name="tests/" + str(i) + "/param.txt")
        except:
            test_failed.append(i)
    if test_failed != []:
        print("Test failed: ", test_failed)
    print("--------- Tests completed ---------")


if __name__ == '__main__':
    testing_process()