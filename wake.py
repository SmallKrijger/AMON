# Install PyWake if needed
import py_wake

# setup site, wind turbines and wind farm model with the corresponding wake models
import numpy as np
from ipywidgets import interact
from ipywidgets import IntSlider
import matplotlib.pyplot as plt
from py_wake.wind_farm_models import All2AllIterative, PropagateDownwind

from py_wake.flow_map import HorizontalGrid
from py_wake.examples.data.iea37._iea37 import IEA37Site, IEA37_WindTurbines
from py_wake.deflection_models.jimenez import JimenezWakeDeflection
from py_wake.deficit_models.no_wake import NoWakeDeficit
from py_wake.deficit_models import BastankhahGaussianDeficit
from py_wake.deficit_models.deficit_model import WakeDeficitModel, BlockageDeficitModel

model = BastankhahGaussianDeficit(use_effective_ws=True) 

wake_deficitModel = [NoWakeDeficit(), model][isinstance(model, WakeDeficitModel)]
blockage_deficitModel = [None, model][isinstance(model, BlockageDeficitModel)]
site = IEA37Site(16)
x, y = [0, 600, 1200], [0, 0, 0]
windTurbines = IEA37_WindTurbines()
wfm = All2AllIterative(site, windTurbines, wake_deficitModel=wake_deficitModel, blockage_deficitModel=blockage_deficitModel, deflectionModel=JimenezWakeDeflection())

# define function that plots the flow field and AEP history of 3 wind turbines
def plot_flow_field_and_aep(WT0, WT1, TILT):

    ax1 = plt.figure(figsize=(20,4)).gca()
    ax2 = plt.figure(figsize=(10,3)).gca()

    sim_res = wfm(x, y, yaw=np.reshape([WT0,WT1,0],(3,1,1)), wd=270, ws=10, tilt=TILT)
    sim_res.flow_map(HorizontalGrid(x = np.linspace(0,1400,200), y=np.linspace(-200,200,50))).plot_wake_map(ax=ax1)
    ax1.set_xlim([-200,1400])
    aep.append(sim_res.aep().values[:,0,0])
    aep_arr = np.array(aep)
    for i in range(3):
        ax2.plot(aep_arr[:,i], '.-', label='WT%d, %.2f'%(i,aep_arr[-1,i]))
    ax2.plot(aep_arr.sum(1), '.-', label='Total, %.2f'%aep_arr[-1].sum())
    ax2.axhline(aep_arr[0].sum(),ls='--',c='r')
    ax2.set_ylabel('AEP [GWh]')
    ax2.set_xlabel('Iteration')
    ax2.legend(loc='upper left')
    plt.show()

# Run the plot_flow_field_and_aep function when moving the sliders
aep=[]
# plot_flow_field_and_aep(50,50,10)

#########################
# PLOT WAKE DEFICIT MAP #
#########################

# Downloaded from: https://topfarm.pages.windenergy.dtu.dk/PyWake/notebooks/WakeDeficitModels.html#NOJDeficit

#here we import all wake deficit models available in PyWake
import numpy as np
import matplotlib.pyplot as plt
import os
import py_wake

from py_wake.examples.data.hornsrev1 import V80, Hornsrev1Site

site = Hornsrev1Site()
windTurbines = V80()

from py_wake import NOJ
from py_wake import Fuga
from py_wake import FugaBlockage
from py_wake import BastankhahGaussian

# Path to Fuga look-up tables
lut_path = os.path.dirname(py_wake.__file__)+'/tests/test_files/fuga/2MW/Z0=0.03000000Zi=00401Zeta0=0.00E+00.nc'

models = {'NOJ': NOJ(site,windTurbines),
          'Fuga': Fuga(lut_path,site,windTurbines),
          'FugaBlockage': FugaBlockage(lut_path,site,windTurbines),
          'BGaus': BastankhahGaussian(site,windTurbines)
         }


from py_wake.superposition_models import LinearSum

models['NOJLinear'] = NOJ(site,windTurbines,superpositionModel=LinearSum())

from py_wake.wind_farm_models import All2AllIterative
from py_wake.deficit_models import NOJDeficit, SelfSimilarityDeficit

models['NOJ_ss'] = All2AllIterative(site,windTurbines,
                                          wake_deficitModel=NOJDeficit(),
                                          superpositionModel=LinearSum(),
                                          blockage_deficitModel=SelfSimilarityDeficit())

from matplotlib import cm
from matplotlib.colors import ListedColormap#, LinearSegmentedColormap
from py_wake.deficit_models.deficit_model import WakeDeficitModel, BlockageDeficitModel
from py_wake.deficit_models.no_wake import NoWakeDeficit
from py_wake.site._site import UniformSite
from py_wake.flow_map import XYGrid
from py_wake.turbulence_models import CrespoHernandez
from py_wake.utils.plotting import setup_plot

#turbine diameter
D = 80

def get_flow_map(model=None, grid=XYGrid(x=np.linspace(-200, 500, 200), y=np.linspace(-200, 200, 200), h=70),
                 turbulenceModel=CrespoHernandez()):
    blockage_deficitModel = [None, model][isinstance(model, BlockageDeficitModel)]
    wake_deficitModel = [NoWakeDeficit(), model][isinstance(model, WakeDeficitModel)]
    wfm = All2AllIterative(UniformSite(), V80(), wake_deficitModel=wake_deficitModel, blockage_deficitModel=blockage_deficitModel,
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

from py_wake.deficit_models import NOJDeficit
plot_wake_deficit_map(NOJDeficit())

plt.savefig( "deficit_map.png" )

