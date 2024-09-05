import itertools
import time

import numpy as np
import windfarm as wf
import initial_solution as init_sol
import matplotlib.pyplot as plt
import PyNomad
import sys
## Utilitaries function 
def aep_func_wt(coords, nb_wt):
    cg, eap = wf.aep_func(coords[nb_wt:], coords[:nb_wt])
    return eap

def create_list_starting_points(nb_wt):
    list_starting_coords = []
    for i in range(30):
        x, y, n_wt = wf.constrained_random(nb_wt_min=nb_wt, nb_wt_max=nb_wt+1) ## pas à jour
        list_starting_coords.append([list(x),list(y)])
    return list_starting_coords

## Function to plot difference between initial solution (feasible)
def plot_NOMAD_CV(total_budget, nb_wt, list_starting_coords, name):
    NOMAD_result = []
    X0 = []
    plt.figure()

    ## Setting parameters for NOMAD
    params = ["BB_INPUT_TYPE * R", "BB_OUTPUT_TYPE OBJ", f"MAX_BB_EVAL {total_budget}","DISPLAY_DEGREE 2", "DISPLAY_ALL_EVAL true", "DISPLAY_STATS BBE OBJ"] 
    params+= ["STATS_FILE" + " stats_test.txt" + " BBE OBJ SOL"]
    params+= ["GRANULARITY * 0.1", "DIRECTION_TYPE ORTHO N+1 QUAD"]
    for i in range(0,nb_wt,2):
        params+= [f"VARIABLE_GROUP {i} {i+1}"]
    lb = [52250, 521700]*nb_wt
    ub = [53800, 523000]*nb_wt
    
    ## Looping through every set of coords
    for i in range(30):
        x, y = list_starting_coords[i][0], list_starting_coords[i][1]
        X0 = sum(map(list, zip(x, y)), [])
        X0 = list(itertools.chain(*zip(x, y))) 

        result = PyNomad.optimize(eap_nomad, X0, lb, ub, params)

        fmt = ["{} = {}".format(n,v) for (n,v) in result.items()]
        output = "\n".join(fmt)
        print("\nNOMAD results \n" + output + " \n")

        ## Read the file
        with open("stats_test.0.txt", "r") as file:
            content = file.readlines()

        # Extract the data we have 3 things per line separated by  spaces
        # The first value is the number of evaluations
        # The second value is the objective function value
        # The other values are the coordinates of the point
        
        # Extract the data
        nb_real_eval = []
        obj_values = []

        for line in content:
            data = line.split()
            nb_real_eval.append(int(data[0]))
            obj_values.append(float(data[1])*10**(-6))

        np_evals = np.array(nb_real_eval)
        np_obj = np.array(obj_values)

        # Getting the success 
        best_obj = np_obj[0]
        best_of = [best_obj]
        best_eval = [0]
        for j,x in enumerate(np_obj):
            if x <= best_obj:
                best_of.append(x)
                best_eval.append(np_evals[j])
                best_obj = x

        best_eval.append(total_budget)
        best_of.append(best_of[-1])
        plt.step(best_eval, best_of, color= cc.colors[i], where="pre", label=str(i))

    plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0, ncol = 2)
    plt.xlabel("Number of function evaluations")
    plt.ylabel("Best Objective function value (GWh)")
    plt.title("Comparison of convergence plot over setted coordinates")
    plt.tight_layout()
    plt.savefig("results/" + str(total_budget) + "_cp_NOMAD_" + name + ".svg")
    plt.show()
    return NOMAD_result

# for i in range(5):
#     NOMAD_result = plot_NOMAD_CV(10, 15, cc.list_start_coords_15, str(i))
# NOMAD_result = plot_NOMAD_CV(10, 15, cc.list_start_coords_15, '0')
# print(NOMAD_result)


def plot_diff(nb_eolienne, nb_evals):
    plt.figure()
    
    for nb_wt in nb_eolienne:
        NOMAD_result = []
        MC_result = []
        NOMAD_time = []
        MC_time = []
        X0 = []
        # Generating an valid initial solution
        x, y, t_Nomad = init_sol.initial_sol_nomad(2000, nb_wt)
        X0 = sum(map(list, zip(x, y)), [])
        X0 = list(itertools.chain(*zip(x, y)))
        
        for total_budget in nb_evals:
            params = ["BB_INPUT_TYPE * R", "BB_OUTPUT_TYPE OBJ EB EB PB", f"MAX_BB_EVAL {total_budget}","DISPLAY_DEGREE 2", "DISPLAY_ALL_EVAL true", "DISPLAY_STATS BBE OBJ"] 
            params+= ["STATS_FILE" + " stats_test.txt" + " BBE OBJ CONS_H SOL"]
            # params+= ["GRANULARITY * 0.1", "DIRECTION_TYPE ORTHO 2N"]
            params+= ["DIRECTION_TYPE ORTHO 2N"]
            lb = [52250, 521700]*nb_wt
            ub = [53800, 523000]*nb_wt

            for i in range(0,nb_wt,2):
                params+= [f"VARIABLE_GROUP {i} {i+1}"]
            
            plt.clf()

            ## Monte Carlo optimization
            t0 = time.time()
            xopt, yopt, nopt, obj_function_opt, eap_set = wf.monte_carlo_cost_dependent(boundary_file=boundary_file, constraints_file=constraints_file, TS_path=TS_path, TS="TS_era.csv", TS_column_name=TS_column_name, DS_path=DS_path, DS="TS_era_direct.csv", DS_column_name=DS_column_name, nsimu=total_budget, obj_function="EAP", nb_wt_min=nb_wt, nb_wt_max=nb_wt+1, cost_factor=0.7, plot_generation=False, plot_flow_map=False, build_graph=False)
            t1 = time.time()
            t_MC = t1-t0
            MC_time.append(t_MC)
            MC_result.append(float(obj_function_opt))

            ## Getting the real successes
            eap_set = -np.array(eap_set)
            best_obj_MC = eap_set[0]
            best_of_MC = [best_obj_MC]
            best_eval_MC = [0]

            for i,x in enumerate(eap_set):
                if x <= best_obj_MC:
                    best_of_MC.append(x)
                    best_eval_MC.append(i)
                    best_obj_MC = x
            best_eval_MC.append(len(eap_set)) # For nice finish of the graph
            best_of_MC.append(best_of_MC[-1]) # //

            plt.clf()
            plt.scatter([i for i in range(len(eap_set))], eap_set, s=10, color='#999bde', marker='d',label="Monte-Carlo")
            plt.step(best_eval_MC, best_of_MC, 'b', where="post")

            ## NOMAD optimization 
            t2 = time.time()
            result = PyNomad.optimize(eap_nomad, X0, lb, ub, params)
            t3 = time.time()
            t_Nomad = t3 - t2
            NOMAD_time.append(t_Nomad)

            np_evals, np_obj, best_eval, best_of = d.read_stat_file(total_budget, stat_file_name="stats_test.0.txt")

            plt.scatter(np_evals, np_obj, color='#d79494', s=10, marker='v', label="NOMAD")
            plt.step(best_eval, best_of, 'r', where='post')
            plt.xlabel("Number of function evaluations")
            plt.ylabel("Best Objective function value (GWh)")
            plt.title("Convergence plot")
            plt.xlim(xmin=0)
            plt.legend()
            plt.tight_layout()
            plt.savefig("results/" + str(total_budget) + "_" + str(nb_wt) + "_cp_post.svg")

            obj_function_value = -result['f_best']*10**(-6)
            NOMAD_result.append(obj_function_value)
            print("Nomad : ", np.round(float(obj_function_value),3), " MC : ", np.round(float(obj_function_opt),3))

            wf.plot_spatial_cstr_generation(result['x_best'][0::2], result['x_best'][1::2], "EAP", "GWh", obj_function_value, nb_wt, Boundaries, boundary_shapely, exclusion_zones_shapely, max_index=max_index, save=True, save_name='NOMAD')
            
        plt.clf()
        plt.subplot(1,2,1)
        plt.plot(nb_evals, NOMAD_result,'#d79494')
        plt.scatter(nb_evals, NOMAD_result, color='r', marker='o', label="NOMAD")
        plt.plot(nb_evals, MC_result, '#999bde')
        plt.scatter(nb_evals, MC_result, color='b', marker='o', label="Monte-Carlo")
        plt.xlabel("Number of evaluations")
        plt.ylabel("Objective function value")
        plt.title('NOMAD and Monte-Carlo results')
        plt.legend()

        plt.subplot(1,2,2)
        plt.plot(nb_evals, NOMAD_time, '#d79494')
        plt.scatter(nb_evals, NOMAD_time, color='r', marker='o', label="NOMAD")
        plt.plot(nb_evals, MC_time, '#999bde')
        plt.scatter(nb_evals, MC_time, color='b', marker='o', label="Monte-Carlo")
        plt.xlabel("Number of evaluations")
        plt.ylabel("Running time (s)")
        plt.title("NOMAD and Monte-Carlo running time")
        plt.legend()
        
        plt.tight_layout()
        plt.savefig("results/" + str(nb_wt) + "_comp.png")
        plt.clf()

# plot_diff(nb_eolienne, nb_evals)

xopt, yopt, nopt, obj_function_opt, eap_set = wf.monte_carlo_cost_dependent(boundary_file="data/boundary_zone.shp", constraints_file="data/exclusion_zones.shp", TS_path="data", TS="wind_speed_1.csv", TS_column_name='0', DS_path="data", DS="wind_direction_1.csv", DS_column_name='0', nsimu=15000, obj_function="EAP", nb_wt_min=30, nb_wt_max=31, cost_factor=0.7, plot_generation=True, plot_flow_map=False, build_graph=False)



def dist_calc():
    fig, ax = plt.subplots()
    # np.set_printoptions(threshold=sys.maxsize)
    # x, y, n_wt = wf.constrained_random(Boundaries, boundary_shapely, exclusion_zones_shapely, nb_wt_min=30, nb_wt_max=30+1, D=D, placing=True)
    # print(n_wt == 90)
    x1, y1, dist_matrix, failed_points = wf.constrained_ramdom_2(D, Boundaries, boundary_shapely, exclusion_zones_shapely, scale_factor=scale_factor, ax=ax, n_wt=90)
    # print(dist_matrix)
    print(wf.compute_number_wrong_wt(boundary_shapely, exclusion_zones_shapely, D, x1, y1))
    # print(failed_points)
    boundary_filled = gpd.GeoSeries(boundary_shapely)
    exclusion_zone_filled = gpd.GeoSeries(exclusion_zones_shapely)
    boundary_filled_index = gpd.GeoSeries(boundary_shapely*len(exclusion_zones_shapely)).boundary

    boundary = boundary_filled.boundary
    exclusion_zone = exclusion_zone_filled.boundary

    ok_zone = boundary_filled
    for polygon in exclusion_zone_filled:
        ok_zone = ok_zone.difference(polygon)
    null_zone_boundaries = boundary_filled_index.intersection(exclusion_zone_filled)

    ax.set_facecolor("lightsteelblue")
    ok_zone.plot(ax=ax, color='lightgreen', alpha=0.5, zorder=1)
    boundary.plot(ax=ax, color=['darkgreen']*len(exclusion_zone_filled), linewidths=1, zorder=2)
    exclusion_zone_filled.plot(ax=ax, color=['gainsboro']*len(exclusion_zones_shapely), zorder=3)
    exclusion_zone.plot(ax=ax, color=['darkgrey']*len(exclusion_zones_shapely), hatch="///", linewidths=1, zorder=5)
    null_zone_boundaries.plot(ax=ax, color=['darkgreen']*len(exclusion_zones_shapely), linestyle='dashed', linewidths=1, zorder=4)
    ax.scatter(x1, y1, marker="o", s=40, color='red', linewidths=1, alpha=0.5, zorder=6, label='Wind Turbine')
    wf.plot_spatial_cstr_generation(x1, y1, "EAP", "GWh", 0, 80, Boundaries, boundary_shapely, exclusion_zones_shapely, max_index, ax=ax, save=True, save_name='example_spacing_placing_not_ok')
    plt.show()

    dist_matrix = wf.distance_matrix(x,y)
    print(dist_matrix)
    wf.plot_spatial_cstr_generation(x, y, "EAP", "GWh", 0, n_wt, Boundaries, boundary_shapely, exclusion_zones_shapely, max_index=max_index, save=True, save_name='TEST')
def dist_calc_2():
    x = 53571.48
    y = 521900.60
    p1 = Point(x, y)
    x1 = 53500.48
    y1 = 522000.60
    p2 = Point(x1, y1)
    x2 = 53671.48
    y2 = 522100.60
    p3 = Point(x2, y2)
    fig, ax = plt.subplots()

    boundary_filled = gpd.GeoSeries(boundary_shapely)
    exclusion_zone_filled = gpd.GeoSeries(exclusion_zones_shapely)
    boundary_filled_index = gpd.GeoSeries(boundary_shapely*len(exclusion_zones_shapely)).boundary

    boundary = boundary_filled.boundary
    exclusion_zone = exclusion_zone_filled.boundary

    ok_zone = boundary_filled
    for polygon in exclusion_zone_filled:
        ok_zone = ok_zone.difference(polygon)
    null_zone_boundaries = boundary_filled_index.intersection(exclusion_zone_filled)

    ax.set_facecolor("lightsteelblue")
    ok_zone.plot(ax=ax, color='lightgreen', alpha=0.5, zorder=1)
    boundary.plot(ax=ax, color=['darkgreen']*len(exclusion_zone_filled), linewidths=1, zorder=2)
    exclusion_zone_filled.plot(ax=ax, color=['gainsboro']*len(exclusion_zones_shapely), zorder=3)
    exclusion_zone.plot(ax=ax, color=['darkgrey']*len(exclusion_zones_shapely), hatch="///", linewidths=1, zorder=5)
    null_zone_boundaries.plot(ax=ax, color=['darkgreen']*len(exclusion_zones_shapely), linestyle='dashed', linewidths=1, zorder=4)
    ax.scatter(x, y, marker="o", s=40, color='red', linewidths=1, alpha=0.5, zorder=6, label='Wind Turbine')
    ax.scatter(x1, y1, marker="o", s=40, color='red', linewidths=1, alpha=0.5, zorder=6, label='Wind Turbine')
    ax.scatter(x2, y2, marker="o", s=40, color='red', linewidths=1, alpha=0.5, zorder=6, label='Wind Turbine')

    plt.show()

    dist = distance(p1, [p1,p2,p3])
    print(dist)

def comp_spiral():
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig, ax = plt.subplots()
    nsimu = 15
    list_time_mc_0 = []
    list_time_mc_1 = []
    list_success_r = []
    list_success_r2 = []
    list_unfeasible_0 = []
    list_unfeasible_1 = []
    nb_wt_min = 40
    nb_wt_max = 120

    for i in range(nb_wt_min, nb_wt_max, 20):
        under_list_r = []
        under_list_uf = []
        t1 = time.time()
        for j in range(0, nsimu):
            print("R1, i = ", i, "n_simu = ", j)
            x, y, n_wt = wf.constrained_random(Boundaries, boundary_shapely, exclusion_zones_shapely, nb_wt_min=i, nb_wt_max=i+1, D=D, placing_spacing=True)
            if n_wt != i:
                under_list_r.append(False)
                under_list_uf.append(len(wf.compute_number_wrong_wt(boundary_shapely, exclusion_zones_shapely, D, x, y)))
            else:
                under_list_r.append(True)
        t2 = time.time()

        list_success_r.append(under_list_r)
        list_unfeasible_0.append(under_list_uf)
        time_mc_classique = (t2 - t1)/nsimu
        list_time_mc_0.append(time_mc_classique)

    for i in range(nb_wt_min, nb_wt_max, 20):
        under_list_r2 = []
        under_list_uf_2 = []
        t3 = time.time()
        for j in range(0, nsimu):
            print("R2, i = ", i, "n_simu = ", j)
            x, y, d, failed = wf.constrained_ramdom_2(D, Boundaries, boundary_shapely, exclusion_zones_shapely, n_wt=i)
            if failed == True:
                under_list_r2.append(False)
                under_list_uf_2.append(len(wf.compute_number_wrong_wt(boundary_shapely, exclusion_zones_shapely, D, x, y)))
            else:
                under_list_r2.append(True)
        t4 = time.time()

        list_success_r2.append(under_list_r2)
        list_unfeasible_1.append(under_list_uf_2)
        time_mc_opt = (t4 - t3)/nsimu
        list_time_mc_1.append(time_mc_opt)

    ax1 = plt.subplot(1, 3, 1)
    ax1.bar([i for i in range(nb_wt_min, nb_wt_max, 20)], list_time_mc_0, width=1, label="Random layout")
    ax1.bar([i+1 for i in range(nb_wt_min, nb_wt_max, 20)], list_time_mc_1, width=1, label="Projected layout")
    plt.xlabel("Number of wind turbines")
    plt.ylabel("Running time (s)")

    ax2 = plt.subplot(1, 3, 2)
    ax2.bar([i for i in range(nb_wt_min, nb_wt_max, 20)], [sum(list_success_r[i]) for i in range(((nb_wt_max-nb_wt_min)//20))], width=1, label="Random layout")
    ax2.bar([i+1 for i in range(nb_wt_min, nb_wt_max, 20)], [sum(list_success_r2[i]) for i in range(((nb_wt_max-nb_wt_min)//20))], width=1, label="Projected layout")
    plt.xlabel("Number of wind turbines")
    plt.ylabel("Number of feasible plot over " + str(nsimu))

    ax3 = plt.subplot(1, 3, 3)
    ax3.bar([i for i in range(nb_wt_min, nb_wt_max, 20)], [(np.mean(list_unfeasible_0[i])/(nb_wt_min + 20*i))*100 for i in range(((nb_wt_max-nb_wt_min)//20))], width=1, label="Random layout")
    ax3.bar([i+1 for i in range(nb_wt_min, nb_wt_max, 20)], [(np.mean(list_unfeasible_1[i])/(nb_wt_min + 20*i))*100 for i in range(((nb_wt_max-nb_wt_min)//20))], width=1, label="Projected layout")
    plt.xlabel("Number of wind turbines")
    plt.ylabel("Percentage of badly placed wind turbine over " + str(nsimu))

    plt.tight_layout()
    plt.savefig("results/random_projected_comp.png")
    plt.show()

## refaire un fichier texte qui parcourt les itérations pour garder les bonnes itérations
## recuperer fonction distance des polygon

# fig, ax = plt.subplots()
# x = 53501.48
# y = 522200.60
# p1 = Point(x, y)
# dist = distance(p1, ok_zone)
# print(dist)
# ok_zone = gpd.GeoSeries(ok_zone)
# ax.set_facecolor("lightsteelblue")
# ax.scatter(x, y, marker="o", s=40, color='red', linewidths=1, alpha=0.5, zorder=6, label='Wind Turbine')
# ok_zone.plot(ax=ax, color='lightgreen', alpha=0.5, zorder=1)


# plt.show()