"""Fonction qui permet de produire les résultats de nomad présentés dans les slides et dans le dossier final rendu"""

import PyNomad
import time
import os
import csv
import numpy as np
import matplotlib.pyplot as plt

def rosenbrock(x,p = 2):
    """Rosenbrock modified function.
    """
    x = np.array(x)
    return np.sum(np.abs(1-x[:-1])**p + 100*np.abs(x[1:] - x[:-1]**2)**p,0)



rosenbrock1 = lambda x: rosenbrock(x,1)
rosenbrock2 = lambda x: rosenbrock(x,2)


# On aurait pu faire un LHS mais le DoE est très petit par rapport à la dimension, son impact est moins important
def generate_initial_data(f,n=25,dim=4,ub = 10, lb = -10, seed=1):
    # generate training data
    list_x = np.random.uniform(lb,ub,(n,dim))
    # take the best one
    best = None
    best_val = np.inf
    for x in list_x:
        val = f(x)
        if val < best_val:
            best = x
            best_val = val
    return best
   

def rosen1_4(x):
    x_list = [x.get_coord(i) for i in range(4)]
    f = rosenbrock1(x_list)
    x.setBBO(str(f).encode("UTF-8"))
    return 1
   
def rosen2_4(x):
    x_list = [x.get_coord(i) for i in range(4)]
    f = rosenbrock2(x_list)
    x.setBBO(str(f).encode("UTF-8"))
    return 1

def rosen1_24(x):
    x_list = [x.get_coord(i) for i in range(24)]
    f = rosenbrock1(x_list)
    x.setBBO(str(f).encode("UTF-8"))
    return 1

def rosen2_24(x):
    x_list = [x.get_coord(i) for i in range(24)]
    f = rosenbrock2(x_list)
    x.setBBO(str(f).encode("UTF-8"))
    return 1

rosen1_24.__name__ = "rosenbrock1_24"
rosen2_24.__name__ = "rosenbrock2_24"
rosen1_4.__name__ = "rosenbrock1_4"
rosen2_4.__name__ = "rosenbrock2_4"

bbs  = {4: [rosen1_4,rosen2_4], 24: [rosen1_24,rosen2_24]}

os.makedirs("py_nomad", exist_ok=True)
proportion_doe = 0.1
for total_budget in [100,250]:
    for dim in [4,24] :
        for quad in [False, True]:
            for bb in bbs[dim]:
                lb = [-10]*dim
                ub = [10]*dim
                # Generate the initial points ans take the best one
                fu = rosenbrock1 if bb.__name__ == "rosenbrock1" else rosenbrock2
                x0 = generate_initial_data(fu, int(proportion_doe*total_budget), dim, lb = lb, ub = ub)

                budget_last = total_budget - int(proportion_doe*total_budget)


                params = ["BB_OUTPUT_TYPE OBJ", f"MAX_BB_EVAL {budget_last}","LOWER_BOUND * -10", "UPPER_BOUND * 10", "DISPLAY_DEGREE 2", "DISPLAY_ALL_EVAL false", "DISPLAY_STATS BBE OBJ"] 
                params+= ["STATS_FILE" + " stats_test.txt" + " BBE OBJ SOL"]
                if quad== False:
                    params+= ["QUAD_MODEL_SEARCH false"]
                else:
                    params+= ["EVAL_QUEUE_SORT DIR_LAST_SUCCESS","DIRECTION_TYPE ORTHO 2N"]

                result = PyNomad.optimize(bb, x0, lb, ub, params)

                fmt = ["{} = {}".format(n,v) for (n,v) in result.items()]
                output = "\n".join(fmt)
                print("\nNOMAD results \n" + output + " \n")

                # Plot the observed data

                # lets retrieve the stats_test.0.txt file
                

                # Read the file
                with open("stats_test.0.txt", "r") as file:
                    content = file.readlines()

                # Extract the data we have 3 things per line separated by  spaces
                # The first value is the number of evaluations
                # The second value is the objective function value
                # The other values are the coordinates of the point
                    
                # Extract the data
                nb_evals = []
                obj_values = []
                x_values = []

                for line in content:
                    data = line.split()
                    nb_evals.append(int(data[0]))
                    obj_values.append(float(data[1]))
                    x_values.append([float(x) for x in data[2:]])

                np_evals = np.array(nb_evals)
                np_obj = np.array(obj_values)
                np_x = np.array(x_values)

                # Plot the data
                folder_base = "py_nomad/"
                folder_function =folder_base+ bb.__name__ + "/"
                os.makedirs(folder_function, exist_ok=True)
                folder_budget = folder_function + f"budget_{total_budget}/"
                os.makedirs(folder_budget, exist_ok=True)
                folder_dim = folder_budget + f"dim_{dim}/"
                os.makedirs(folder_dim, exist_ok=True)
                if quad:
                    plt.clf()
                    plt.plot(np_evals, np_obj, label="Quad")
                    minimum = np.min(np_obj)
                    plt.text(0, minimum, f"min: {minimum:.2f}")
                    plt.xlabel("Number of evaluations")
                    plt.ylabel("Objective function value")
                    plt.legend()
                    plt.savefig(folder_dim+f"{bb.__name__}_quad_{quad}.png")
                    plt.clf()
                    plt.plot(np_evals, np_obj, label="Quad")
                    minimum = np.min(np_obj)
                    plt.text(0, minimum, f"min: {minimum:.2f}")
                    plt.xlabel("Number of evaluations")
                    plt.ylabel("Objective function value")  
                    plt.yscale("log")
                    plt.savefig(folder_dim+f"{bb.__name__}_log_quad_{quad}.png")
                    #save the data in a csv file
                    with open(folder_dim+f"{bb.__name__}_quad_{quad}.csv", "w") as file:
                        writer = csv.writer(file)
                        writer.writerow(["evals", "obj"])
                        for e, o in zip(np_evals, np_obj):
                            writer.writerow([e,o])
                else:
                    plt.clf()
                    plt.plot(np_evals, np_obj)
                    minimum = np.min(np_obj)
                    plt.text(0, minimum, f"min: {minimum:.2f}")
                    plt.xlabel("Number of evaluations")
                    plt.ylabel("Objective function value")
                    plt.savefig(folder_dim+f"{bb.__name__}_quad_{quad}.png")
                    plt.clf()
                    plt.plot(np_evals, np_obj)
                    minimum = np.min(np_obj)
                    plt.text(0, minimum, f"min: {minimum:.2f}")
                    plt.yscale("log")
                    plt.xlabel("Number of evaluations")
                    plt.ylabel("Objective function value")
                    plt.savefig(folder_dim+f"{bb.__name__}_log_quad_{quad}.png")
                    #save the data in a csv file
                    with open(folder_dim+f"{bb.__name__}_quad_{quad}.csv", "w") as file:
                        writer = csv.writer(file)
                        writer.writerow(["evals", "obj"])
                        for e, o in zip(np_evals, np_obj):
                            writer.writerow([e,o])