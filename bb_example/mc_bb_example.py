import os
import sys
dir = os.path.basename(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', dir))
import windfarm_eval
import data as d
import numpy as np
import itertools

## Monte Carlo solver
def test_solver(instance, x0, nb_eval):
    eap_0, s_d_0, sum_dist_0 = windfarm_eval.windfarm_eval(instance, x0)
    best_eap = eap_0
    for _ in range(nb_eval):
        print(best_eap)
        new_x = np.random.uniform(52274.51767349267, 53803.380417069864, 30)
        new_y = np.random.uniform(521785.03108614404, 522850.77981024975, 30)
        X = sum(map(list, zip(new_x, new_y)), [])
        X = list(itertools.chain(*zip(new_x, new_y)))
        eap_1, s_d_1, sum_dist_1 = windfarm_eval.windfarm_eval(instance, X)
        if eap_1 > eap_0:
            best_eap = eap_1
            eap_0 = eap_1
    return best_eap

if __name__ == '__main__':
    best_eap = test_solver('bb_example/param.txt', 'bb_example/x0_instance.txt', 10)
    print("Best EAP = ", best_eap, "GWh")