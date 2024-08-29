from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np

def plot_terrain(x, y, obj_function, units, obj_function_value, n_wt, ub, boundary_shapely, exclusion_zones_shapely, max_index="", cg="", ax="", plot_flow_map=False, full_wind_rose=False, save=False, save_name=""):
    if ax == "":    
        fig, ax = plt.subplots()

    boundary_filled = gpd.GeoSeries(boundary_shapely)
    boundary = boundary_filled.boundary
    ok_zone = boundary_filled
    ax.set_facecolor("lightsteelblue")

    if exclusion_zones_shapely != []:
        exclusion_zone_filled = gpd.GeoSeries(exclusion_zones_shapely)
        boundary_filled_index = gpd.GeoSeries(boundary_shapely*len(exclusion_zones_shapely)).boundary
        exclusion_zone = exclusion_zone_filled.boundary
        for polygon in exclusion_zone_filled:
            ok_zone = ok_zone.difference(polygon)
            null_zone_boundaries = boundary_filled_index.intersection(exclusion_zone_filled)
        ok_zone.plot(ax=ax, color='lightgreen', alpha=0.5, zorder=1)
        exclusion_zone_filled.plot(ax=ax, color=['gainsboro']*len(exclusion_zones_shapely), zorder=3)
        exclusion_zone.plot(ax=ax, color=['darkgrey']*len(exclusion_zones_shapely), hatch="///", linewidths=1, zorder=5)
        null_zone_boundaries.plot(ax=ax, color=['darkgreen']*len(exclusion_zones_shapely), linestyle='dashed', linewidths=1, zorder=4)
        ax.scatter(x, y, marker="o", s=40, color='red', linewidths=1, alpha=0.5, zorder=6, label='Wind Turbine')

    else:
        ok_zone.plot(ax=ax, color='lightgreen', alpha=0.5, zorder=1)
        ax.scatter(x, y, marker="o", s=40, color='red', linewidths=1, alpha=0.5, zorder=3, label='Wind Turbine')
    
    if isinstance(boundary_shapely, list): 
        boundary.plot(ax=ax, color=['darkgreen']*len(boundary_shapely), linewidths=1, zorder=2)
    else:
        boundary.plot(ax=ax, color=['darkgreen'], linewidths=1, zorder=2)
    plt.title( str(obj_function + " = %s " + units + ", Nb wind turbines : %s")%(obj_function_value, n_wt))

    if full_wind_rose:
        wr_plot = inset_axes(ax,
                        width="20%", # width = 20% of parent_bbox
                        height="29%", 
                        loc=1)
        wr_plot.patch.set_edgecolor('black')  
        wr_plot.patch.set_linewidth(2) 
        im = plt.imread('data/WindRose.png')
        wr_plot.imshow(im)
        wr_plot.axis('off')
    
    else:
        u = np.cos(np.radians(max_index*10))
        v = np.sin(np.radians(max_index*10))
        ax.quiver(ub[0]-100, ub[1]-100, u, v, 40, angles=270-(max_index*10), scale=15)

    if plot_flow_map:
        cg.flow_map().plot_wake_map() 

    plt.xlabel("X")
    plt.ylabel("Y")
    plt.tight_layout()
    ax.legend(loc='lower left')

    if save:
        plt.savefig(save_name)

def plot_result_nomad(np_evals, np_obj, best_eval, best_of, total_budget, nb_wt, save_name):
    plt.scatter(np_evals, np_obj, color='#d79494', s=10, marker='v', label="NOMAD")
    plt.step(best_eval, best_of, 'r', where='post')
    plt.xlabel("Number of function evaluations")
    plt.ylabel("Best Objective function value (GWh)")
    plt.title("Convergence plot")
    plt.xlim(xmin=0)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_name)
