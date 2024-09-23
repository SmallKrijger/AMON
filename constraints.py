"""

This script contains the functions to handle and compute the multiple constraints used during the optimization process.

This script requires multiple libraries that are written in the `requirements.txt` to be installed in your Python environnement. 
Make sure to install them properly with the right version.

"""

from shapely.geometry import Point
from shapely import distance, difference, equals_exact

def creating_list_point(x_coords, y_coords):
    """Script that creates a list of Shapely Point object.

    Parameters
    ----------
    x_coords : list
        List of x coordinates.
    y_coords : list
        List of y coordinates.
    
    Returns
    -------
    list_points : list
        List of Shapely Point object. 
    """

    list_points = []
    nb_wt = len(x_coords)
    for i in range(nb_wt):
        list_points.append(Point(x_coords[i], y_coords[i]))
    return list_points

def buildable_zone(boundary_shapely, exclusion_zones_shapely):
    """Script that generates the actual buildable zone by taking out the exclusion zones from the boundary zone.

    Parameters
    ----------
    boundary_shapely : list
        List of Shapely Polygon of the boundary zone. 
    exclusion_zones_shapely : list
        List of Shapely Polygon of the exclusion zone. 
    
    Returns
    -------
    ok_zone : Shapely MultiPolygon
        Actual buildable zone. 
    """

    ok_zone = boundary_shapely
    for zone in exclusion_zones_shapely:
        ok_zone = difference(ok_zone, zone)
    return ok_zone

def dist_matrix(list_points):
    """Script that computes the matrix with the distance between every wind turbines.

    Parameters
    ----------
    list_points : list
        List of Shapely Point (representing the wind turbines)
    
    Returns
    -------
    d_matrix : list
        List of distances between every wind turbines.
    """

    nb_wt = len(list_points)
    d_matrix = [distance(list_points[i], list_points) for i in range(nb_wt)]
    for i in range(nb_wt):
        for j in range(0, i):
            if (j != i):
                d_matrix[i][j] = 0
    return d_matrix

def placing_constraint(x_coords, y_coords, ok_zone):
    """Script that computes the sum of distances of every wind turbines to the closest edge of the buildable zone MultiPolygon (the distance is set to 0
    if the wind turbine is already in the MultiPolygon).

    Parameters
    ----------
    x_coords : list
        List of x coordinates.
    y_coords : list
        List of y coordinates.
    ok_zone : Shapely MultiPolygon
        Actual buildable zone. 
    
    Returns
    -------
    sum_dist : float
        Sum of distances of every wind turbines to the closest edge of the buildable zone MultiPolygon.
    """

    list_points = creating_list_point(x_coords, y_coords)
    dist = distance(list_points, ok_zone)
    sum_dist = sum(dist)
    return sum_dist

def spacing_constraint_min(x_coords, y_coords, D):
    """Script that computes the sum of distances between wind turbines if they are too close to each other.

    Parameters
    ----------
    x_coords : list
        List of x coordinates.
    y_coords : list
        List of y coordinates.
    D : int
        Minimal distance between two wind turbines.
    
    Returns
    -------
    -s_d : float
        Sum of distances between every wind turbines if they are too close.
    """

    list_points = creating_list_point(x_coords, y_coords)
    d_matrix = dist_matrix(list_points)
    s_d = 0
    for list_d in d_matrix:
        for d in list_d:
            if (d != 0):
                s_d += min(d - 2*D, 0)
    return -s_d

def checking_same_coords(x_coords, y_coords):
    """Script that checks if two wind turbines have the same coordinates.

    Parameters
    ----------
    x_coords : list
        List of x coordinates.
    y_coords : list
        List of y coordinates.
    
    Returns
    -------
    boolean : boolean
        Boolean set to False if two wind turbines have the same coordinates, True otherwise.
    """

    list_points = creating_list_point(x_coords, y_coords)
    for i, point in enumerate(list_points):
        if i != len(list_points):
            list_bool = point.equals_exact(list_points[i+1:], 1e-9)
            s = sum(list_bool)
            if s > 0:
                return False
    return True
