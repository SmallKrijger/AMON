"""

This script contains the functions to generate initial solutions respecting the constraints.

This script requires multiple libraries that are written in the `requirements.txt` to be installed in your Python environnement. 
Make sure to install them properly with the right version.

"""

import itertools
import sys
import numpy as np
import windfarm_setting as wf
import constraints as cst
import data as d
import PyNomad

def initial_sol_nomad(nb_wt, buildable_zone, D, lb, ub, spacing=True):
    """Script that generates an initial solution respecting the constraints using the NOMAD solver.

        Parameters
        ----------
        nb_wt : int
            Number of wind turbines.
        buildable_zone : Shapely MultiPolygon
            Actual buildable zone. 
        D : int
            Diameter of the wind turbine.
        lb : list
            Lower bound of the boundary zone on the x and y axis.
        ub : list
            Upper bound of the boundary zone on the x and y axis.
        spacing : boolean
            If True then the spacing constraint is taken into account for the initial solution. 
        
        Returns
        -------
        x_best : list
            List of x coordinates of the initial solution.
        y_best : list
            List of y coordinates of the initial solution.
    """

    def initial_sol(x):
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

            # Calculate constraints
            sum_dist = cst.placing_constraint(x_coords, y_coords, buildable_zone)

            if spacing:
                s_d = cst.spacing_constraint_min(x_coords, y_coords, D)
                rawBBO = str(0) + " " + str(s_d) + " " + str(sum_dist)
                x.setBBO(rawBBO.encode("UTF-8"))
            else:
                rawBBO = str(0) + " " + str(sum_dist)
                x.setBBO(rawBBO.encode("UTF-8"))

        except:
            print("Unexpected eval error", sys.exc_info()[0])
            return 0
        return 1
    
    params = ["BB_INPUT_TYPE * R", "STOP_IF_FEASIBLE true", "DISPLAY_DEGREE 0", "DISPLAY_ALL_EVAL false", "DIRECTION_TYPE ORTHO 2N"] 

    if spacing:
        params += ["BB_OUTPUT_TYPE OBJ PB PB"]
    else:
        params += ["BB_OUTPUT_TYPE OBJ PB"]

    for i in range(0, nb_wt, 2):
        params += [f"VARIABLE_GROUP {i} {i+1}"]

    x = np.random.uniform(lb[0], ub[0], nb_wt)
    y = np.random.uniform(lb[1], ub[1], nb_wt)
    X0 = sum(map(list, zip(x, y)), [])
    X0 = list(itertools.chain(*zip(x, y)))

    ## NOMAD optimization 
    result = PyNomad.optimize(initial_sol, X0, lb*nb_wt, ub*nb_wt, params)
    
    x_best = result['x_best'][0::2]
    y_best = result['x_best'][1::2]
    return x_best, y_best

def initial_sol_test(param_file_name, x0_file_name):
    """Script to compute the eap and constraint values with the NOMAD solver.

        Parameters
        ----------
        param_file_name : str
            Path to the parameter text file.
        x0_file_name : str
            Path where the function will create the x0.txt with the coordinates of the initial solution.
        
        Returns
        -------
        x0_file_name\x0.txt : text file
            A text file with the coordinates of the initial solution.
    """
    
    # Initializing site and boundary files
    nb_wt, D, hub_height, scale_factor, power_curve, boundary_file, exclusion_zone_file, wind_speed, wind_direction = d.read_param_file(param_file_name)
    lb, ub, boundary_shapely, exclusion_zones_shapely = wf.terrain_setting(boundary_file, exclusion_zone_file, scale_factor=scale_factor)
    buildable_zone = cst.buildable_zone(boundary_shapely, exclusion_zones_shapely)
    x, y = initial_sol_nomad(nb_wt, buildable_zone, D, lb, ub, spacing=True)
    X0 = sum(map(list, zip(x, y)), [])
    X0 = list(itertools.chain(*zip(x, y)))
    x0_file_path = x0_file_name + "/x0.txt"
    f = open(x0_file_path, 'w+')  # open file in write mode
    f.write(str(X0))
    f.close()   
            