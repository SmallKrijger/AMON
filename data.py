import numpy as np
import matplotlib.pyplot as plt

def read_param_file(param_file_name):
    # Read the file
    with open(param_file_name, "r") as file:
        content = file.read().splitlines() 

    nb_wt = int(content[0].split()[1])
    diameter = int(content[1].split()[1])
    hub_height = int(content[2].split()[1])
    scale_factor = float(content[3].split()[1])
    power_curve = content[4].split()[1]
    boundary_file = content[5].split()[1]
    exclusion_zone_file = content[6].split()[1]
    wind_speed = content[7].split()[1]
    wind_direction = content[8].split()[1]
               
    return nb_wt, diameter, hub_height, scale_factor, power_curve, boundary_file, exclusion_zone_file, wind_speed, wind_direction

# nb_wt, D, boundary_file, exclusion_zone_file, wind_speed, wind_direction = read_param_file("tests/1/param.txt")
# print(nb_wt, D, boundary_file, exclusion_zone_file, wind_speed, wind_direction)

def read_config_file(config_file_name, nb_wt):
    params = []
    # Read the file
    with open(config_file_name, "r") as file:
        content = file.read().splitlines()
    nb_it = int(content[3].split()[1])
    params = content
    for i in range(0, nb_wt, 2):
        params.append(f"VARIABLE_GROUP {i} {i+1}")
    return params, nb_it
    
# params = read_config_file("data/config.txt", nb_wt)
# print(params)

def read_stat_file(total_budget, stat_file_name="nomad_result.0.txt"):
    # Read the file
    with open(stat_file_name, "r") as file:
        content = file.readlines()

    # Extract the data we have 3 things per line separated by  spaces
    # The first value is the number of evaluations
    # The second value is the objective function value
    # The third value is the constraint value
    # The other values are the coordinates of the point
        
    # Extract the data
    nb_real_eval = []
    obj_values = []

    for line in content:
        data = line.split()
        if data[2] != 'inf':
            nb_real_eval.append(int(data[0]))
            obj_values.append(float(data[1])*10**(-6))

    np_evals = np.array(nb_real_eval)
    np_obj = np.array(obj_values)

    ## Getting the real successes
    best_obj = obj_values[0]
    best_of = [best_obj]
    best_eval = [0]
    for i,x in enumerate(obj_values):
        if x <= best_obj:
            best_of.append(x)
            best_eval.append(np_evals[i])
            best_obj = x

    best_eval.append(total_budget) # For nice finish of the graph
    best_of.append(best_of[-1]) # //
    return np_evals, np_obj, best_eval, best_of

def draw_result_nomad(np_evals, np_obj, best_eval, best_of, total_budget, nb_wt, save_name):
    plt.scatter(np_evals, np_obj, color='#d79494', s=10, marker='v', label="NOMAD")
    plt.step(best_eval, best_of, 'r', where='post')
    plt.xlabel("Number of function evaluations")
    plt.ylabel("Best Objective function value (GWh)")
    plt.title("Convergence plot")
    plt.xlim(xmin=0)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_name)