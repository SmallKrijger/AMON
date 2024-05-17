import csv
import os
import numpy as np
import matplotlib.pyplot as plt

forder_base = "final_curves/"

dims = [4, 24]
nb_runs = [100, 250]
ps = [1, 2]

def extract_data(csv_file) :
    with open(csv_file, mode='r') as file:
        reader = csv.reader(file)
        time = []
        value = []
        for i,row in enumerate(reader):
            if i == 0 or len(row) == 0:
                continue
            # print(row)
            time.append(float(row[0]))
            value.append(float(row[1]))
    return time, value



algos = ["nomad","nomad_quad", "botorch_single_output", "botorch_multi_output" ]
dims = [4, 24]
nb_runs = [100, 250]
ps = [1, 2]


results = {}
# results["botorch_multi_output"] = {}
# results["botorch_single_output"] = {}
# results["nomad_quad"] = {}
# results["nomad"] = {}

for algo in algos: # botorch_multi_output, botorch_single_output, nomad_quad, nomad
    results[algo] = {}
    for dim in dims:# 4 ou 24
        results[algo][dim] = {}
        for run in nb_runs: # 100 ou 250
            results[algo][dim][run] = {}
            for p in ps: # 1 ou 2
                results[algo][dim][run][p] = {}
                for time_or_eval_f_value in ["time", "eval","f_value"]:
                    results[algo][dim][run][p][time_or_eval_f_value] = []

# Lets retrieve the data from nomad

folder_nomad ="py_nomad/"
for quad in ["False", "True"]:
    for dim in dims:
        for run in nb_runs:
            for p in ps:
                csv_file = folder_nomad + f"rosenbrock{p}_{dim}/budget_{run}/dim_{dim}/rosenbrock{p}_{dim}_quad_{quad}.csv"
                eval, value = extract_data(csv_file)
                eval.append(run)
                value.append(value[-1])
                results["nomad_quad" if quad == "True" else "nomad"][dim][run][p]["f_value"] = value
                results["nomad_quad" if quad == "True" else "nomad"][dim][run][p]["eval"] = eval


# Lets retrieve the data from botorch single output

folder_botorch_single = "scp_folder/botorch_simple_output/"
for dim in dims:
    for run in nb_runs:
        for p in ps:
            csv_file = folder_botorch_single + f"better_imgs_ei{dim}_results/rosenbrock{p}/rosenbrock{p}_{run}_{dim}_results.csv"
            time, value = extract_data(csv_file)
            eval = list(range(1, len(time)+1))
            value = -np.maximum.accumulate(value)
            results["botorch_single_output"][dim][run][p]["f_value"] = value
            results["botorch_single_output"][dim][run][p]["eval"] = eval
            results["botorch_single_output"][dim][run][p]["time"] = time

# Lets retrieve the data from botorch multi output

folder_botorch_multi = "scp_folder/botorch_multi_output/"
for dim in dims:
    for run in [100] :
        for p in ps:
            csv_file = folder_botorch_multi + f"better_imgs_ei{dim}_results/h{p}/h{p}_{run}_{dim}_results.csv"
            time, value = extract_data(csv_file)
            eval = list(range(1, len(time)+1))
            value = -np.maximum.accumulate(value)
            results["botorch_multi_output"][dim][run][p]["f_value"] = value
            results["botorch_multi_output"][dim][run][p]["eval"] = eval
            results["botorch_multi_output"][dim][run][p]["time"] = time

forlder_botorch_multi2 = "scp_folder/better_imgs_ei4_results/"
for run in [250]:
    for p in ps:
        csv_file = forlder_botorch_multi2 + f"h{p}/h{p}_{run}_4_results.csv"
        time, value = extract_data(csv_file)
        eval = list(range(1, len(time)+1))
        value = -np.maximum.accumulate(value)
        results["botorch_multi_output"][4][run][p]["f_value"] = value
        results["botorch_multi_output"][4][run][p]["eval"] = eval
        results["botorch_multi_output"][4][run][p]["time"] = time



# Lets plot the 8 graphs

final_folder = "final_rendu/"
os.makedirs(final_folder, exist_ok=True)

def plot_from_result(dim,run,p,time_or_eval,show=False) :
    # Lets plot the 4 curves for the 4 algorithms together
    value_nomad = results["nomad"][dim][run][p]["f_value"]
    value_nomad_quad = results["nomad_quad"][dim][run][p]["f_value"]
    value_botorch_single = results["botorch_single_output"][dim][run][p]["f_value"]
    value_botorch_multi = results["botorch_multi_output"][dim][run][p]["f_value"]
    list_values = [value_nomad, value_nomad_quad, value_botorch_single, value_botorch_multi]
    time_or_eval_nomad = results["nomad"][dim][run][p][time_or_eval]
    time_or_eval_nomad_quad = results["nomad_quad"][dim][run][p][time_or_eval]
    time_or_eval_botorch_single = results["botorch_single_output"][dim][run][p][time_or_eval]
    time_or_eval_botorch_multi = results["botorch_multi_output"][dim][run][p][time_or_eval]
    list_time_or_eval = [time_or_eval_nomad, time_or_eval_nomad_quad, time_or_eval_botorch_single, time_or_eval_botorch_multi]

    plt.clf()
    for i in range(4):
        if list_time_or_eval[i] == []:
            print(f"Missing data for {algos[i]}")
            continue
        if i == 0:
            print("Nomad ", list_time_or_eval[i])
        plt.step(list_time_or_eval[i], list_values[i], label=algos[i])
        plt.xlabel('CPU time' if time_or_eval == "time" else 'Number of evaluations')
        plt.ylabel('Best observed value')
        plt.title(f"Minimum of rosenbrock{p}, dim = {dim}, {run} runs")
    plt.legend()
    plt.savefig(final_folder+f"rosenbrock{p}_dim_{dim}_run_{run}_p_{p}_{time_or_eval}.png")
    if show:
        plt.show()
    plt.yscale('log')
    folder_log = final_folder + "log/"
    os.makedirs(folder_log, exist_ok=True)
    plt.savefig(folder_log+f"rosenbrock{p}_dim_{dim}_run_{run}_p_{p}_{time_or_eval}_log.png")
    # plt.legend(algos)
        
for dim in dims:
    for run in nb_runs:
        for p in ps:
            plot_from_result(dim,run,p,"time")
            plot_from_result(dim,run,p,"eval")
# test

# plot_from_result(4,250,1,"time",show=True)




