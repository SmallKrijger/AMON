import time
import PyNomad

#############################
# General imports and paths #
#############################
import os
from IPython.display import display, clear_output
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
from py_wake.superposition_models import SquaredSum
from py_wake.site import XRSite
from py_wake.wind_turbines import WindTurbine
from py_wake.wind_turbines.power_ct_functions import PowerCtTabular
import xarray as xr
from py_wake.site._site import UniformSite
import shapefile
from shapely.geometry import Polygon, Point, MultiPolygon
from py_wake.flow_map import XYGrid
from py_wake.turbulence_models import CrespoHernandez
from py_wake.deficit_models.deficit_model import WakeDeficitModel, BlockageDeficitModel
from py_wake.deficit_models import BastankhahGaussianDeficit
from py_wake.deficit_models.no_wake import NoWakeDeficit
from py_wake.wind_farm_models import All2AllIterative
from py_wake.utils.plotting import setup_plot
from matplotlib.colors import ListedColormap
from matplotlib import cm
import geopandas as gpd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# specify where to load data, and save results
data_path = "data"
results_path = "results"

def setting_wind_turbine(powercurve_path, diameter, hub_height):
    power_curve = pd.read_csv(powercurve_path, sep=";", skiprows=0)
    u = power_curve.WindSpeed.values
    power = power_curve.Power.values*1e3
    ct = power_curve.Ct.values
    windturbine = WindTurbine(name='windturbine',
                        diameter=diameter,
                        hub_height=hub_height,
                        powerCtFunction=PowerCtTabular(u, power,'kW', ct))
    return windturbine

def reduce_file_size(file_path, years, hours):
    file = pd.read_csv(file_path, index_col=0)
    file.index = pd.to_datetime(file.index)
    file = file[file.index.year >= years]
    file = file[file.index.hour >= hours]
    return file

###############################
# Read Wind Speed Time-Series #
###############################


def site_model(powercurve_path, diameter, hub_height, WS_path, WD_path):
    
    ## Creating Turbine
    Turbine = setting_wind_turbine(powercurve_path, diameter, hub_height)

    ## Reading files
    WS = pd.read_csv(WS_path, index_col=0)
    WD = pd.read_csv(WD_path, index_col=0)
    
    ## Creating dataframe for wind rose values
    list_of_tuples = [[None for i in range(36)] for j in range(41)]
    wd_values = np.array([i*10 for i in range(36)])
    ws_values = np.array([0] + [0.5+i for i in range(41)])
    df = pd.DataFrame(list_of_tuples, columns=wd_values, index=ws_values[1:])
    
    N = len(WS)
    width = 360 / len(wd_values)
    for i in range(len(wd_values)):
        wd = wd_values[i]
        for j in range(len(ws_values)-1):
            lower, upper = ws_values[j], ws_values[j+1]
            if wd == 0:
                sector = (360 - 0.5*width <= WD.values) | (WD.values < 0.5*width)
            else:
                sector = (wd - 0.5*width <= WD.values) & (WD.values < wd + 0.5*width)
            TS_sector = WS.values[sector]
            df.iloc[j,i] = sum((lower <= TS_sector) & (TS_sector < upper)) / N

    ## Create site object
    site = XRSite( 
        ds =xr.Dataset(
            data_vars={"P":(("wd", "ws"), df.values.T), "TI":.1},
            coords={"wd":wd_values, "ws":ws_values[1:]}
            )
        )

    wd_tot_per_ws = []
    for i in range(len(wd_values)):
        wd_tot_per_ws.append(site.ds.P.values[i].sum())

    max_index = np.argmax(wd_tot_per_ws)

    ## Creating WindRose and saving it
    script_dir = os.path.dirname(__file__)
    results_dir = os.path.join(script_dir, 'results/')

    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    fig, ax = plt.subplots()
    site.plot_wd_distribution(n_wd=36, ax=ax)
    fig.set_size_inches(3,3)
    plt.savefig("results/WindRose.png", dpi=130)
    plt.close()

    model = BastankhahGaussianDeficit(use_effective_ws=True) 
    blockage_deficitModel = [None, model][isinstance(model, BlockageDeficitModel)]
    wake_deficitModel = [NoWakeDeficit(), model][isinstance(model, WakeDeficitModel)]
    fmGROSS = All2AllIterative(site, Turbine, wake_deficitModel=wake_deficitModel, blockage_deficitModel=blockage_deficitModel, superpositionModel=SquaredSum(), turbulenceModel=CrespoHernandez())
    return fmGROSS, WS, WD, max_index, wd_tot_per_ws[max_index]

# site, fmGROSS, TS, DS = site_model(plot_distribution=False)

#######################
# SPATIAL CONSTRAINTS #
#######################

def spatial_constraints(boundary_file, constraints_file, scale_factor=0.1):
    
    # Read domain boundaries, convert to shapely format
    Boundaries = shapefile.Reader(boundary_file)
    boundary_shapely = []

    for shape in Boundaries.shapes():
        coords = np.array(shape.points).T*scale_factor
        boundary_shapely.append(Polygon(coords.T))

    boundary_shapely = [MultiPolygon(boundary_shapely)]
    
    exclusion_zones_shapely = []
    if constraints_file != "na":
        # Read exclusion zones boundaries, convert to shapely format
        Constraints = shapefile.Reader(constraints_file)

        for shape in Constraints.shapes():
            coords = np.array(shape.points).T*scale_factor
            exclusion_zones_shapely.append(Polygon(coords.T))

    lb = [(Boundaries.bbox[0])*scale_factor, (Boundaries.bbox[1])*scale_factor]
    ub = [(Boundaries.bbox[2])*scale_factor, (Boundaries.bbox[3])*scale_factor]
    # print(Boundaries.bbox[0], Boundaries.bbox[1], Boundaries.bbox[2], Boundaries.bbox[3])
    
    return lb, ub, boundary_shapely, exclusion_zones_shapely

# print(Boundaries, boundary_shapely, Constraints, exclusion_zones_shapely)

def plot_spatial_cstr_generation(x, y, obj_function, units, obj_function_value, n_wt, ub, boundary_shapely, exclusion_zones_shapely, max_index="", cg="", ax="", plot_flow_map=False, full_wind_rose=False, save=False, save_name=""):
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
        im = plt.imread('results/WindRose.png')
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

def plot_wake_example(Turbine, D=80):

    def get_flow_map(model=None, grid=XYGrid(x=np.linspace(-200, 500, 200), y=np.linspace(-200, 200, 200), h=70),
                 turbulenceModel=CrespoHernandez()):
        blockage_deficitModel = [None, model][isinstance(model, BlockageDeficitModel)]
        wake_deficitModel = [NoWakeDeficit(), model][isinstance(model, WakeDeficitModel)]
        wfm = All2AllIterative(UniformSite(), Turbine, wake_deficitModel=wake_deficitModel, blockage_deficitModel=blockage_deficitModel,
                            turbulenceModel=turbulenceModel)
        return wfm(x=[0], y=[0], wd=270, ws=10, yaw=0).flow_map(grid)

    def plot_deficit_map(model, cmap='Blues', levels=np.linspace(0, 10, 55)):
        fm = get_flow_map(model)
        data = fm.ws - fm.WS_eff
        fm.plot(data[0,:,:,:,0], clabel='Deficit [m/s]', levels=levels, cmap=cmap, normalize_with=D)
        setup_plot(grid=False, ylabel="Crosswind distance [y/D]", xlabel= "Downwind distance [x/D]",
                xlim=[fm.x.min()/D, fm.x.max()/D], ylim=[fm.y.min()/D, fm.y.max()/D])#, axis='auto')

    def plot_wake_deficit_map(model):
        cmap = np.r_[[[1,1,1,1],[1,1,1,1]],cm.Blues(np.linspace(-0,1,128))] # ensure zero deficit is white
        plot_deficit_map(model,cmap=ListedColormap(cmap))

    plot_wake_deficit_map(BastankhahGaussianDeficit(use_effective_ws=True))
    plt.show()

# plot_wake_example()
    
#########################
# CREATE THE WAKE MODEL #
#########################

def read_csv(WS_BB, WD_BB):
    WS = pd.read_csv(WS_BB, index_col=0)
    WS_column_name = list(WS.columns)
    WD = pd.read_csv(WD_BB, index_col=0)
    WD_column_name = list(WD.columns)
    return WS[WS_column_name[0]], WD[WD_column_name[0]]

# Expected annual production
def aep_func(x, y, fmGROSS, WS_BB, WD_BB):
    cg = fmGROSS(x, y, ws=WS_BB, wd=WD_BB, time=True)
    eap = cg.aep().sum()
    return cg, eap

#################################################################
# UNIFORM SIMULATION WITHIN DOMAIN / EXCLUSION ZONES BOUNDARIES # 
#################################################################

def spiral_function(theta, x_center, y_center, scale_factor):
    x = theta*np.cos(theta)*scale_factor
    y = theta*np.sin(theta)*scale_factor
    return x+x_center, y+y_center

def distance_matrix(x, y):
    nb_wt = len(x)
    dist_matrix = [[0 for _ in range(nb_wt)] for _ in range(nb_wt)]
    for i in range(nb_wt):
        for j in range(i, nb_wt):
            if (j != i):
                dist_matrix[i][j] = np.sqrt((x[i] - x[j])**2 + (y[i] - y[j])**2) 
    return np.array(dist_matrix)

def total_distance_between_wt(x, y):
    dist_matrix = distance_matrix(x, y)
    d = dist_matrix.sum().sum()
    return d

def spacing_constraint_all(Boundaries, boundary_shapely, exclusion_zones_shapely, ax, x, y, D, scale_factor):
    nb_wt = len(x)
    dist_matrix = distance_matrix(x, y)
    failed = False
    safety_margin_x = 1
    safety_margin_y = 1.4
    x_min = Boundaries.bbox[0]*scale_factor
    x_max = Boundaries.bbox[2]*scale_factor
    y_min = Boundaries.bbox[1]*scale_factor
    y_max = Boundaries.bbox[3]*scale_factor
    x_dim = (x_max - x_min)*safety_margin_x
    y_dim = (y_max - y_min)*safety_margin_y
    for i in range(nb_wt):
        list_try_x = []
        list_try_y = []
        for j in range(i, nb_wt): 
            if (j != i):
                nb_iter = 0
                while (dist_matrix[i][j] < D) or (test_constraints_placing(boundary_shapely, exclusion_zones_shapely, [x[i]], [y[i]])[0] == False):
                    x_proj, y_proj = spiral_function(nb_iter*(np.pi//2), x[i], y[i], 1.5)
                    # list_try_x.append(x_proj)
                    # list_try_y.append(y_proj)

                    if (x_proj < x_min - x_dim) or (x_proj > x_max + x_dim) or (y_proj < y_min - y_dim) or (y_proj > y_max + y_dim):
                        failed = True
                        # ax.plot(list_try_x, list_try_y, zorder=8, alpha=0.5)    
                        break
                    else:
                        if test_constraints_placing(boundary_shapely, exclusion_zones_shapely, [x_proj], [y_proj])[0]:
                            d_proj_i_k = []
                            for k in range(nb_wt):
                                if k != i:
                                    d_proj_i_k.append((np.sqrt((x_proj - x[k])**2 + (y_proj - y[k])**2) >= D))

                            if (sum(d_proj_i_k) != nb_wt-1):
                                nb_iter += 1                        
                            else:
                                for k in range(nb_wt):
                                    if k != i:
                                        dist_matrix[i][k] = np.sqrt((x_proj - x[k])**2 + (y_proj - y[k])**2)
                                x[i], y[i] = x_proj, y_proj 
                        else:
                            nb_iter += 1 
            elif (i == nb_wt - 1):
                nb_iter = 0
                while (test_constraints_placing(boundary_shapely, exclusion_zones_shapely, [x[i]], [y[i]])[0] == False):
                    x_proj, y_proj = spiral_function(nb_iter*(np.pi//2), x[i], y[i], 1.5)
                    # list_try_x.append(x_proj)
                    # list_try_y.append(y_proj)

                    if (x_proj < x_min - x_dim) or (x_proj > x_max + x_dim) or (y_proj < y_min - y_dim) or (y_proj > y_max + y_dim):
                        failed = True
                        # ax.plot(list_try_x, list_try_y, zorder=8, alpha=0.5)    
                        break
                    else:
                        if test_constraints_placing(boundary_shapely, exclusion_zones_shapely, [x_proj], [y_proj])[0]:
                            d_proj_i_k = []
                            for k in range(nb_wt):
                                if k != i:
                                    d_proj_i_k.append((np.sqrt((x_proj - x[k])**2 + (y_proj - y[k])**2) >= D))

                            if (sum(d_proj_i_k) != nb_wt-1):
                                nb_iter += 1                        
                            else:
                                for k in range(nb_wt):
                                    if k != i:
                                        dist_matrix[i][k] = np.sqrt((x_proj - x[k])**2 + (y_proj - y[k])**2)
                                x[i], y[i] = x_proj, y_proj 
                        else:
                            nb_iter += 1 
        # ax.plot(list_try_x, list_try_y, zorder=8, alpha=0.5)            
    return x, y, dist_matrix, failed

def compute_number_wrong_wt(boundary_shapely, exclusion_zones_shapely, D, x, y):
    dist_matrix = distance_matrix(x, y)
    nb_wt = len(x)
    unfeasible_wt = []
    for i in range(nb_wt):
        if (test_constraints_placing(boundary_shapely, exclusion_zones_shapely, [x[i]], [y[i]])[0] == False):
            unfeasible_wt.append(i)
        for j in range(i, nb_wt):
            if (j != i):
                if (dist_matrix[i][j] < D):
                    unfeasible_wt.append((i, j))      
    return unfeasible_wt

def projection_point(boundary_shapely, exclusion_zones_shapely, ax, x, y):
    for i in range(len(x)):
        position_bool = test_constraints_placing(boundary_shapely, exclusion_zones_shapely, [x[i]], [y[i]])[0]
        k = 0
        list_try_x = []
        list_try_y = []
        while position_bool == False:
            x_proj, y_proj = spiral_function(k*(np.pi//2), x[i], y[i], 10)
            list_try_x.append(x_proj)
            list_try_y.append(y_proj)
            if test_constraints_placing(boundary_shapely, exclusion_zones_shapely, [x_proj], [y_proj])[0]:
                x[i], y[i] = x_proj, y_proj
                position_bool = True
            else:
                k += 1
        # ax.plot(list_try_x, list_try_y)
    return x, y

def test_constraints_placing(boundary_shapely, exclusion_zones_shapely, x, y):
    # initialize test to True
    test = np.repeat(True, len(x))

    for i in range(len(x)):
        # points must be in windfarm boundaries
        for polygon in boundary_shapely:
            test[i] *= polygon.contains( Point( x[i], y[i] ) ) 
        # points must not be in exclusion zones
        for polygon in exclusion_zones_shapely:
            test[i] *= not polygon.contains( Point( x[i], y[i] ) )
    return test 

def test_constraints_placing_spacing(D, boundary_shapely, exclusion_zones_shapely, x, y, x_acc, y_acc):
    """
    Tests whether 2D points respect spatial constraints

    Parameters
    ----------
    x : (n) array
        x cordinates
    y : (n) array
        y cordinates

    Returns
    -------
    test (n) boolean array
        test[i] is True iff (x[i],y[i]) respect spatial constraints
    """
    test = test_constraints_placing(boundary_shapely, exclusion_zones_shapely, x, y)
    # points must be at a reasonable distance (no overlap) 
    for i in range(len(x)): 
        for j in range(i,len(x)):  
            if (j != i) and test[i]:
                d = np.sqrt((x[i] - x[j])**2 + (y[i] - y[j])**2)
                test[i] *= (d >= D)

    if len(x_acc) != 0:
        for i in range(len(x)):
            for j in range(len(x_acc)):
                if test[i]:
                    d = np.sqrt((x[i] - x_acc[j])**2 + (y[i] - y_acc[j])**2)
                    test[i] *= (d >= D)
    return test

def constrained_random(Boundaries, boundary_shapely, exclusion_zones_shapely, nb_wt_min=1, nb_wt_max=10, D=80, scale_factor=0.1, placing=False, placing_spacing=False, verbose=False):
    """
    Simulate a random number of wind turbines under spatial constraints

    Parameters
    ----------
    verbose : bool, optional
        Whether to print convergence messages. The default is False
    TODO
    Returns
    -------
    x_gen : (nsimu) array
        x coordinates 
    y_gen : (nsimu) array
        y coordinates 
    TODO
    """
    n_wt = np.random.randint(nb_wt_min, nb_wt_max)
    x_gen, y_gen = np.zeros((2,n_wt))
    nacc = 0
    x_acc, y_acc = [], []
    max_limit = 120
    start = time.time()

    while time.time() - start < max_limit:
        if nacc < n_wt:
            x = np.random.uniform( Boundaries.bbox[0], Boundaries.bbox[2], n_wt-nacc )*scale_factor
            y = np.random.uniform( Boundaries.bbox[1], Boundaries.bbox[3], n_wt-nacc )*scale_factor

            if placing_spacing:
                acc = test_constraints_placing_spacing(D, boundary_shapely, exclusion_zones_shapely, x, y, x_acc, y_acc)

            if placing:
                acc = test_constraints_placing(boundary_shapely, exclusion_zones_shapely, x, y)

            nnew = acc.sum()
            for xi in x[acc]:
                x_acc.append(xi)
            for yi in y[acc]:
                y_acc.append(yi)

            x_gen[nacc:nacc+nnew], y_gen[nacc:nacc+nnew] = x[acc], y[acc]
            nacc += nnew
            if verbose:
                print("nb. accepted = %s"%nacc)
        else:
            return x_gen, y_gen, n_wt
    
    for i in range(len(x_gen)):
        if x_gen[i] == 0:
            x_new = np.random.uniform( Boundaries.bbox[0], Boundaries.bbox[2], 1)*scale_factor
            y_new = np.random.uniform( Boundaries.bbox[1], Boundaries.bbox[3], 1)*scale_factor
            # while not test_constraints_placing(boundary_shapely, exclusion_zones_shapely, [x_new], [y_new])[0]:
            #     x_new = np.random.uniform( Boundaries.bbox[0], Boundaries.bbox[2], 1)*scale_factor
            #     y_new = np.random.uniform( Boundaries.bbox[1], Boundaries.bbox[3], 1)*scale_factor 
            x_gen[i] = x_new
            y_gen[i] = y_new
    return x_gen, y_gen, nacc

def constrained_random_2(D, Boundaries, boundary_shapely, exclusion_zones_shapely, scale_factor=0.1, ax="", n_wt=10):
    x = np.random.uniform( Boundaries.bbox[0], Boundaries.bbox[2], n_wt)*scale_factor
    y = np.random.uniform( Boundaries.bbox[1], Boundaries.bbox[3], n_wt)*scale_factor
    x_gen, y_gen = projection_point(boundary_shapely, exclusion_zones_shapely, ax, x, y)
    x_final, y_final, dist_matrix, failed_points = spacing_constraint_all(Boundaries, boundary_shapely, exclusion_zones_shapely, ax, x_gen, y_gen, D, scale_factor)
    # x_gen, y_gen, dist_matrix, failed_points = spacing_constraint_all(Boundaries, boundary_shapely, exclusion_zones_shapely, ax, x, y, D, scale_factor)
    return x_final, y_final, dist_matrix, failed_points
    # return x, y

#################
# ONE GENERATION #
#################

def cost_dependent_plot_random(boundary_file='Test2.shp', constraints_file='Test2_Constraints.shp', TS_path='data', TS="TS_era.csv", TS_column_name='0', DS_path='data', DS="TS_era_direct.csv", DS_column_name='0', obj_function="EAP", nb_wt_min=1, nb_wt_max=10, cost_wt_kw=2500, cost_factor=1, n_years=25, c_MWh=60, scale_factor=0.1, plot_flow_map=False, plot_generation=False):
    """
    Calls the constrained_MC_random function. Plots the result configuration.

    Parameters
    ----------
    TODO

    Returns
    -------
    TODO
    """
    # Generating site and fmGROSS
    site, fmGROSS, TS, DS = site_model(TS_path=TS_path, TS=TS, TS_column_name=TS_column_name, DS_path=DS_path, DS=DS, DS_column_name=DS_column_name)

    # Generating the spatial constraint
    Boundaries, boundary_shapely, exclusion_zones_shapely = spatial_constraints(boundary_file=boundary_file, constraints_file=constraints_file, scale_factor=scale_factor)

    # Generating a configuration
    x, y, n_wt = constrained_random(Boundaries, boundary_shapely, exclusion_zones_shapely, nb_wt_min=nb_wt_min, nb_wt_max=nb_wt_max, scale_factor=scale_factor, placing=True)

    # Expected annual production
    cg, eap = aep_func(x, y, fmGROSS=fmGROSS, TS=TS, DS=DS, TS_column_name=TS_column_name, DS_column_name=DS_column_name) 

    units = ""
    obj_function_value = 0

    if obj_function == "EAP":
        units = "MWh"
        obj_function_value = np.round(1e3*float(eap),3)
        
    elif obj_function == "ROI":
        units = "€"
        cost = ((1 + n_wt)**cost_factor - 1) * cost_wt_kw * Pmax
        obj_function_value = np.round((1e3*float(eap)*n_years*c_MWh) - cost,1)

    elif obj_function == "LCOE":
        units = "€/MWh"
        cost = ((1 + n_wt)**cost_factor - 1) * cost_wt_kw * Pmax
        # Levelized Cost Of Energy
        obj_function_value = np.round(cost/(1e3*float(eap)*n_years*c_MWh),1)

    else:
        return "Error: obj_function must be EAP, ROI or LCOE."
    
    if plot_flow_map:
        cg.flow_map().plot_wake_map() 
    if plot_generation:
        plot_spatial_cstr_generation(x, y, obj_function, units, obj_function_value, n_wt, boundary_shapely, exclusion_zones_shapely, save=True, save_name="Test") 
    
    return obj_function_value, n_wt
    # plt.savefig( "results/random_design_" + obj_function + ".png" )

# obj_function_value, n_wt = cost_dependent_plot_random(boundary_file='Test2.shp', constraints_file='Test2_Constraints.shp', TS="TS_era.csv", TS_column_name='0', DS="TS_era_direct.csv", obj_function="EAP", nb_wt_min=19, nb_wt_max=20, cost_factor=0.7, plot_flow_map=False, plot_generation=True)

############################
# Monte-Carlo Optimization #
############################

def storing_value(fmGROSS, TS, TS_column_name, DS, DS_column_name, Boundaries, boundary_shapely, exclusion_zones_shapely, test_set, cg_set, eap_set, wt_set, nb_wt_min, nb_wt_max, scale_factor):
    """
    Calls the constrained_MC_random function and store the values it returns. Used for each iteration of the monte_carlo_cost_dependent function.

    Parameters
    ----------
    TODO

    Returns
    -------
    TODO
    """
    x, y, n_wt = constrained_random(Boundaries, boundary_shapely, exclusion_zones_shapely, nb_wt_min=nb_wt_min, nb_wt_max=nb_wt_max, scale_factor=scale_factor, placing_spacing=True)
    test_set.append( (x, y) )

    # Expected annual production and storing it
    cg, eap = aep_func(x, y, fmGROSS=fmGROSS, TS=TS, DS=DS, TS_column_name=TS_column_name, DS_column_name=DS_column_name)
    eap_set.append(eap)
    cg_set.append(cg)
    
    # Storing number of wind turbines for each generation
    wt_set.append(n_wt)

    return n_wt, eap, test_set, cg_set, eap_set, wt_set

def monte_carlo_cost_dependent(boundary_file='Test2.shp', constraints_file='Test2_Constraints.shp', TS_path='data', TS="TS_era.csv", TS_column_name='0', DS_path='data', DS="TS_era_direct.csv", DS_column_name='0', nsimu=10, obj_function="EAP", nb_wt_min=1, nb_wt_max=500, cost_wt_kw=2500, cost_factor=1, n_years=25, c_MWh=60, scale_factor=0.1, plot_generation=False, plot_flow_map=False, build_graph=False):
    """
    Generates a number nsimu of configuration and returns the best one.

    Parameters
    ----------
    TODO

    Returns
    -------
    TODO
    """
    test_set = []
    cg_set = []
    eap_set = []
    wt_set = []
    nb_wt = 0

    # Generating site and fmGROSS
    site, fmGROSS, TS, DS, max_index, wd_max = site_model(TS_path=TS_path, TS=TS, TS_column_name=TS_column_name, DS_path=DS_path, DS=DS, DS_column_name=DS_column_name)

    # Generating the spatial constraint
    Boundaries, boundary_shapely, exclusion_zones_shapely = spatial_constraints(boundary_file=boundary_file, constraints_file=constraints_file, scale_factor=scale_factor)

    if obj_function == "EAP":
        eap_max = 0
        
        for i in range(nsimu):
            n_wt, eap, test_set, cg_set, eap_set, wt_set = storing_value(fmGROSS, TS, TS_column_name, DS, DS_column_name, Boundaries, boundary_shapely, exclusion_zones_shapely, test_set, cg_set, eap_set, wt_set, nb_wt_min, nb_wt_max, scale_factor)
            
            # clear_output(wait=True)
            if eap_max <= float(eap):
                eap_max = float(eap)
                nb_wt = n_wt
            else:
                pass

            # display('Simu %s, EAP max = %s Mwh, Nb wind turbines = %s '%(i+1, eap_max, nb_wt))

        iopt = np.argmax(eap_set)
        units = "GWh"
        xopt, yopt = test_set[iopt]
        nopt = wt_set[iopt]
        cg_opt = cg_set[iopt]
        
        if build_graph:
            eap_set = np.array(eap_set)
            plt.figure()
            plt.scatter(wt_set, eap_set)
            plt.xlabel("Number of Wind Turbines")
            plt.ylabel("EAP (GWh)")
            plt.show()
            # plt.savefig( "results/EAP_behaviour.png" )

        if plot_generation:
            if plot_flow_map:
                cg_opt.flow_map().plot_wake_map()
            plot_spatial_cstr_generation(xopt, yopt, obj_function, units, np.round(float(eap_set[iopt]),3), nopt, boundary_shapely, exclusion_zones_shapely, save=True, save_name='MC') 

        return xopt, yopt, nopt, eap_set[iopt], eap_set

    elif obj_function == "ROI":
        roi_set = []
        roi_max = 0

        for i in range(nsimu):
            n_wt, eap, test_set, cg_set, eap_set, wt_set = storing_value(fmGROSS, TS, TS_column_name, DS, DS_column_name, Boundaries, boundary_shapely, exclusion_zones_shapely, test_set, cg_set, eap_set, wt_set, nb_wt_min, nb_wt_max, scale_factor)

            cost = ((1 + n_wt)**cost_factor - 1) * cost_wt_kw * Pmax

            roi = np.round(1e3*float(eap)*n_years*c_MWh - cost,1)
            roi_set.append(roi)

            clear_output(wait=True)
            if roi_max <= roi:
                roi_max = roi
                nb_wt = n_wt
            else:
                pass

            display('Simu %s, ROI max = %s €, Nb wind turbines = %s '%(i+1, roi_max, nb_wt))
        
        iopt = np.argmax(roi_set)
        units = "€"
        xopt, yopt = test_set[iopt]
        nopt = wt_set[iopt]
        cg_opt = cg_set[iopt]

        if build_graph:
            roi_set = np.array(roi_set)
            plt.figure()
            plt.scatter(wt_set, roi_set)
            plt.xlabel("Number of Wind Turbines")
            plt.ylabel("Return on Investment (€)")
            plt.show()
            # plt.savefig( "results/ROI_behaviour.png" )

        if plot_generation:
            if plot_flow_map:
                cg_opt.flow_map().plot_wake_map()
            plot_spatial_cstr_generation(xopt, yopt, obj_function, units, np.round(roi_set[iopt], 2), nopt, boundary_shapely, exclusion_zones_shapely) 
        
        return xopt, yopt, nopt, roi_set[iopt]

    elif obj_function == "LCOE":
        lcoe_set = []
        cost_set = []
        lcoe_min = 1e12

        for i in range(nsimu):
            n_wt, eap, test_set, cg_set, eap_set, wt_set = storing_value(fmGROSS, TS, TS_column_name, DS, DS_column_name, Boundaries, boundary_shapely, exclusion_zones_shapely, test_set, cg_set, eap_set, wt_set, nb_wt_min, nb_wt_max, scale_factor)
            
            cost = ((1 + n_wt)**cost_factor - 1) * cost_wt_kw * Pmax
            cost_set.append(cost)

            # Levelized Cost Of Energy
            lcoe = np.round(cost/(1e3*float(eap)*n_years*c_MWh),3)
            lcoe_set.append(lcoe)

            clear_output(wait=True)
            if lcoe_min >= lcoe:
                lcoe_min = lcoe
                nb_wt = n_wt
            else:
                pass

            display('Simu %s, LCOE min = %s €/MWh, Nb wind turbines = %s '%(i+1, lcoe_min, nb_wt))

        iopt = np.argmin(lcoe_set)
        units = "€/MWh"
        xopt, yopt = test_set[iopt]
        nopt = wt_set[iopt]
        cg_opt = cg_set[iopt]

        if build_graph:
            eap_set = np.array(eap_set)
            plt.figure()
            plt.scatter(cost_set, eap_set)
            for i, txt in enumerate(lcoe_set):
                plt.text(cost_set[i]+.2, eap_set[i]-.2, txt, fontsize=12)
            plt.xlabel("Cost (€)")
            plt.ylabel("EAP (GWh)")
            plt.show()
            # plt.savefig( "results/LCOE_behaviour.png" )
        
        if plot_generation:
            if plot_flow_map:
                cg_opt.flow_map().plot_wake_map()
            plot_spatial_cstr_generation(xopt, yopt, obj_function, units, np.round(lcoe_set[iopt],2), nopt, boundary_shapely, exclusion_zones_shapely) 

        return xopt, yopt, nopt, lcoe_set[iopt]

    else:
        return "Error: obj_function must be EAP, ROI or LCOE."
    
# xopt, yopt, nopt, obj_function_opt, useful_set = monte_carlo_cost_dependent(nsimu=10, obj_function="EAP", nb_wt_min=15, nb_wt_max=16, cost_factor=0.7, plot_generation=True, plot_flow_map=False, build_graph=False)
# print(obj_function_opt)

import fnmatch
def uncertainties_calc(xopt, yopt, nopt, simu_ws_dir='simu_ws', TS_column_name='simulated_0',simu_wd_dir='simu_wd', DS_column_name='0',cost_wt_kw=2500, cost_factor=1, n_years=25, c_MWh=60):
    eap_unc = []
    cost_unc = []
    lcoe_unc = []
    cost = ((1 + nopt)**cost_factor - 1) * cost_wt_kw * Pmax

    # simu_wd_dir = os.path.join('Rework',simu_wd_dir)
    # print(simu_wd_dir)
    # simu_ws_dir = os.path.join('Rework',simu_ws_dir)

    list_name_wd = fnmatch.filter(os.listdir(simu_wd_dir), '*.csv')
    list_name_ws = fnmatch.filter(os.listdir(simu_ws_dir), '*.csv')

    file_nb_wd = len(list_name_wd)
    file_nb_ws = len(list_name_ws)
    
    if file_nb_wd == 1:
        for name_ws in list_name_ws:
            cg, eap = aep_func(xopt, yopt, TS_path=simu_ws_dir, TS=name_ws, TS_column_name=TS_column_name, DS_path=simu_wd_dir, DS=list_name_wd[0], DS_column_name=DS_column_name)
            eap_unc.append(np.round(float(eap),3))
            
            eap_convert = 1e3*float(eap)*n_years*c_MWh
            cost_unc.append(np.round(eap_convert - cost,1))
            lcoe_unc.append(np.round(cost/eap_convert,3))

    elif file_nb_ws == 1:
        for name_wd in list_name_wd:
            cg, eap = aep_func(xopt, yopt, TS_path=simu_ws_dir, TS=list_name_ws[0], TS_column_name=TS_column_name, DS_path=simu_wd_dir, DS=name_wd, DS_column_name=DS_column_name)
            eap_unc.append(np.round(float(eap),3))
            
            eap_convert = 1e3*float(eap)*n_years*c_MWh
            cost_unc.append(np.round(eap_convert - cost,1))
            lcoe_unc.append(np.round(cost/eap_convert,3))

    elif file_nb_ws == file_nb_wd:
        pass
    else:
        return 'Error: You can set a unique wind direction/wind speed csv file or have the same number of random generation for wind speed and wind direction.'

    return eap_unc, cost_unc, lcoe_unc

# eap_unc, cost_unc, lcoe_unc = uncertainties_calc(xopt, yopt, nopt)
# print(eap_unc, cost_unc, lcoe_unc)

def plot_uncertainties(useful_opt, useful_set, set_name):
    useful_set = np.round(np.array(useful_set), 3)
    mean_set = np.round(useful_set.mean(), 3)
    median_set = np.round(np.median(useful_set),3)
    q_10_set = np.round(np.quantile(useful_set, 0.1),3)
    q_90_set = np.round(np.quantile(useful_set, 0.9),3)
    abs = np.array([i for i in range(len(useful_set))])

    
    plt.scatter(abs, useful_set)
    plt.plot(len(useful_set)+1, useful_opt, 'ro') 
    plt.axhline(y = mean_set, color = 'r', linestyle = '--', label = 'Mean') 
    plt.axhline(y = median_set, color = 'b', linestyle = 'dashed', label = 'Median') 
    plt.axhline(y = q_10_set, color = 'g', linestyle = ':', label = 'Quantile 10')
    plt.axhline(y = q_90_set, color = 'purple', linestyle = '-', label = 'Quantile 90')  
    plt.title(str(set_name + " mean: %s, Median " + set_name + ": %s, Quantile 10 " + set_name + ": %s, Quantile 90 " + set_name + ": %s")%(mean_set, median_set, q_10_set, q_90_set))
    plt.show()
    plt.savefig( str("results/" + set_name + "_uncertainties.png") )

# plot_uncertainties(obj_function_opt, eap_unc, "EAP")