# AMON

<p> AMON is a blackbox optimization benchmarking framework for windfarm layout. </p>

## Installation

<p> Before running the different scripts, you should create a conda environment with <b>Python version 3.9.19</b> using the following command:</p>

```bash
conda create -n your_env_name python=3.9.19
conda activate your_env_name
```
where `your_env_name` is the name you want to give to your environment.

**<span style="color:red;"> The `requirements.txt` file contains the required python librairies to run the scripts. </span>**
If you want additional librairies in your environment, make sure to edit this file before running the following instruction. 

**<span style="color:red;"> Please install the needed python libraries in your environment with the following command: </span>**
```bash
pip install -r requirements.txt
```

## Validation 

<p> Your installation must be validated. For this, type : </p>

```bash
python ./windfarm_nomad.py
```

The results of this function will be written in the freshly created `instances_results` directory and the file `output_instances.txt`. All instances should result in a lign in the `output_instances.txt` file.
If the validation fails, please send an email to josephine.gobert@hotmail.com with the full output.

## Execution

<p> To run a simulation with your own solver, you can either call directly the blackbox function with the constraints in a Python file using the following code: </p>

```bash
import windfarm_eval
eap, spacing_constraint, placing_constraint = windfarm_eval.windfarm_eval('path/to/your/param.txt', 'path/to/your/X.txt')
```

<p> Or you can call it in a terminal using the following command: </p>

```bash
python ./windfarm_eval.py 'path/to/your/param.txt' 'path/to/your/X.txt'
```

The `path/to/your/param.txt` file is the instance parameter file and should follow the format such as the one in the `instances` folder:

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
The `path/to/your/X.txt` is the path to the initial set of points for the blackbox to begin with. It should follow the same format as the one in the `instances` folder:

```bash
[x_0, y_0, x_1, y_1, ..., x_n, y_n]
```

<p> It is to be noted that the first iteration will be slower, but the following one will be quicker due to storage in a cache of the site and the terrain. </p>

## Execution example

You can call the blackbox with a parameter and x file from one of the instances with the following command in a terminal:
```bash
python ./windfarm_eval.py 'instances/1/param.txt' 'instances/1/x0.txt'
```

or in a Python file:
```bash
import windfarm_eval
eap, spacing_constraint, placing_constraint = windfarm_eval.windfarm_eval('instances/1/param.txt', 'instances/1/x0.txt')
```

You will find in the `bb_example` folder, a simple Python code for a Monte Carlo type of solver that calls the blackbox a define number of time. To execute it, launch it in a terminal using the following command:

```bash
python bb_example/mc_bb_example.py
```
