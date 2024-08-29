import os
import sys
import windfarm_opt as wf_opt
import data as d
import numpy as np
import itertools

## Optimiseur type Monte Carlo
def optimiseur_test(instance, x0, nb_eval):
    if os.path.exists("bb_example\results\bb_result.txt"):
        os.remove("bb_example\results\bb_result.txt")
    eap_0, s_d_0, sum_dist_0 = wf_opt.aep(instance, x0)
    best_eap = eap_0
    for _ in range(nb_eval):
        print(best_eap)
        new_x = np.random.uniform(52274.51767349267, 53803.380417069864, 30)
        new_y = np.random.uniform(521785.03108614404, 522850.77981024975, 30)
        X = sum(map(list, zip(new_x, new_y)), [])
        X = list(itertools.chain(*zip(new_x, new_y)))
        eap_1, s_d_1, sum_dist_1 = wf_opt.aep(instance, X)
        if eap_1 > eap_0:
            best_eap = eap_1
            eap_0 = eap_1
    return best_eap

if __name__ == '__main__':
    param_path = sys.argv[1]     # Get param file path
    x_path = sys.argv[2]         # Get initial solution
    best_eap = optimiseur_test(param_path, x_path, 10)
    print("Best EAP = ", best_eap, "GWh")