# WindFarm Layout Optimization

<p> This is a blackbox optimization benchmarking framework for windfarm layout. </p>

## Installation

<p> Before running the different scripts, you should create a conda environment with <b>Python version 3.9.19</b> using the following command:</p>

```bash
conda create -n <your_env_name> python=3.9.19
conda activate <your_env_name>
```

<p>The requirements.txt file contains the required python librairies to run the scripts. If you want additional librairies in your environment, make sure to edit this file before running the following instruction.</p>

<p>Install the needed python libraries :</p>

```bash
pip install -r requirements.txt
```

## Validation 

<p> Your installation must be validated. For this, type : </p>

```bash
python .\running_test.py
```

The results of this function will be written in the freshly created `tests_results` directory and the file `output_tests.txt`. All tests should result in a 
lign in the `output_tests.txt` file.
If the validation fails, please send an email to josephine.gobert@hotmail.com with the full output.

## Execution

<p> To run a simulation with your own solver, you can either call directly the blackbox functions with the constraints in a Python file using the following code: </p>

```bash
import windfarm_opt
eap, spacing_constraint, placing_constraint = windfarm_opt.aep('.\instance.txt', '.\X.txt')
```


<p> Or you can call it in a terminal using the following command: </p>

```bash
python .\windfarm_opt.py .\instance.txt .\X.txt
```

Your `.\instance.txt` file is the instance parameter file and should follow the format such as the one in the `tests` folder:

```bash
DIMENSION               Number of wind turbines in your problem
WT_DIAMETER             Diameter of your wind turbine
WT_HUB_HEIGHT           Height of your wind turbine hub
SCALE_FACTOR            Scale factor to adjust the size of your terrain (usually 1)
POWER_CURVE             Path to the power curve data of your wind turbine
BOUNDARY_FILE           Path to the boundary shapefile
EXCLUSION_ZONE_FILE     Path to the exclusion zones shapefile (na if you have no exclusion zones shapefile)
WIND_SPEED              Path to the wind speed csv data
WIND_DIRECTION          Path to the wind direction csv data
```
`.\X.txt` is the path to the initial set of points for the blackbox to begin with. It should follow the same format as the one in the `tests` folder:

```bash
[x_0, y_0, x_1, y_1, ..., x_n, y_n]
```

<p> It is to be noted that the first iteration will be slower, but the following one will be quicker due to storage in a cache of the site and the terrain. </p>