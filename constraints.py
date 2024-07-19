from shapely.geometry import Polygon, Point
from shapely import distance, difference
import windfarm as wf

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

def spacing_constraint(x_coords, y_coords, D):
    list_points = creating_list_point(x_coords, y_coords)
    d_matrix = dist_matrix(list_points)
    count_wt = 0
    for list_d in d_matrix:
        for d in list_d:
            if (d != 0) and (d <= D):
                count_wt += 1
    return count_wt

def spacing_constraint_min(x_coords, y_coords, D):
    list_points = creating_list_point(x_coords, y_coords)
    d_matrix = dist_matrix(list_points)
    s_d = 0
    for list_d in d_matrix:
        for d in list_d:
            if (d != 0):
                s_d += min(d - D, 0)
    return -s_d

def multiple_constr(x_coords, y_coords):
    list_points = creating_list_point(x_coords, y_coords)
    d_matrix = dist_matrix(list_points)
    cstr = " "
    for list_d in d_matrix:
        for d in list_d:
            if d != 0:
                cstr = cstr + " " + str(-d + 80)
    return str(cstr)

# x_coords = [53571.48, 53500.48, 53671.48]       
# y_coords = [521900.60, 522000.60, 522100.60]
# cst = multiple_constr(x_coords, y_coords)
# print(cst)
# d_matrix = spacing_constraint(x_coords, y_coords)
# d_matrix1 = wf.distance_matrix(x_coords, y_coords)
# print(d_matrix)
# print(d_matrix1)

# dist_boundary, dist_exclusion = placing_constraint(x_coords, y_coords, wf.boundary_shapely, wf.exclusion_zones_shapely)
# print(dist_boundary, dist_exclusion)
## if dist_boundary > 0 then not in bound
## if dist_exclusion > 0 then in exclusion zones
