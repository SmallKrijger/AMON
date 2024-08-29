from shapely.geometry import Point
from shapely import distance, difference, equals_exact

def creating_list_point(x_coords, y_coords):
    list_points = []
    nb_wt = len(x_coords)
    for i in range(nb_wt):
        list_points.append(Point(x_coords[i], y_coords[i]))
    return list_points

def buildable_zone(boundary_shapely, exclusion_zones_shapely):
    ok_zone = boundary_shapely
    for zone in exclusion_zones_shapely:
        ok_zone = difference(ok_zone, zone)
    return ok_zone

def dist_matrix(list_points):
    nb_wt = len(list_points)
    d_matrix = [distance(list_points[i], list_points) for i in range(nb_wt)]
    for i in range(nb_wt):
        for j in range(0, i):
            if (j != i):
                d_matrix[i][j] = 0
    return d_matrix

def placing_constraint(x_coords, y_coords, ok_zone):
    list_points = creating_list_point(x_coords, y_coords)
    dist = distance(list_points, ok_zone)
    sum_dist = sum(dist)
    return sum_dist

def spacing_constraint_min(x_coords, y_coords, D):
    list_points = creating_list_point(x_coords, y_coords)
    d_matrix = dist_matrix(list_points)
    s_d = 0
    for list_d in d_matrix:
        for d in list_d:
            if (d != 0):
                s_d += min(d - D, 0)
    return -s_d

def checking_same_coords(x_coords, y_coords):
    list_points = creating_list_point(x_coords, y_coords)
    for i, point in enumerate(list_points):
        if i != len(list_points):
            list_bool = point.equals_exact(list_points[i+1:], 1e-9)
            s = sum(list_bool)
            if s > 0:
                return False
    return True
