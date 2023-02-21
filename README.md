# Monte-Carlo-simulation
This repo contains the code used for Monte Carlo Simulation in DBT molecules using trained machine learning model

## Installation
First clone the repo
```
git clone https://github.com/ChunHou20c/Monte-Carlo-simulation.git
```

then create an virtual environment for the project 
```
python3 -m venv ./env
```
activate the virtual environment
```
source activate ./env/bin/activate
#depends on your shell the environment file might be different
```
for windows:
```
.\env\scripts\activate
```

install all the requirement packages
```
pip install -r requirements.txt
```

## Usage
Run the script run_simulation.py
```
python run_simulation.py
```

## Modifications of cache path
The simulation will create cache file which stores all the coulombic matrix automatically when the simulation object is first constructed.
This will reduce the simulation time when the simulation need to be run again.

By default the cache files are located in cache directory.

Modify CACHE_PATH in .env file to change the path for the cache

## Running simulation for different molecule
The simulation can be run for 3 molecules: DBT1, DBT2 and DBT3 as defined in simulation package
**The complete script is the same as in run_simulation.py, the parts may be modified are shown below

### import the required packages
```python
from simulation import simulation
import pandas as pd
from simulation import DBT1 # you may change to DBT2 or DBT3 
import os
from decouple import config
```
If you want to run for DBT2 or DBT3, make sure the molecule type and prediction model are also for DBT2 and DBT3

### Define the simulation class
```python
Simulation = simulation.Define_simulation_model('cnn_dbt1') # change cnn_dbt1 to use random model and for DBT2 and DBT3
```
options: 
- cnn_dbt1 - trained cnn model for dbt1 molecule using tensorflow
- random_model_dbt1 - random model for dbt1
- random_model_dbt2 - random model for dbt2
- random_model_dbt3 - random model for dbt3

** There are only random model available for DBT2 and DBT3 molecule

### Construct simulation object
```python
test_sim = Simulation('Frames/DBT1-static.gro', DBT1.DBT1, memory_saving=True, cache_path=config('CACHE_PATH'))
```
This is the part that will affected by .env file

### Create a path to store the results
```python
path = "result_dbt1"
# Check whether the specified path exists or not
isExist = os.path.exists(path)
if not isExist:

   # Create a new directory because it does not exist
   os.makedirs(path)
```
### Running the simulation
```python
for i in range(10): #run for 10 times
    distance_list, time_list = test_sim.run()
    data_to_save = pd.DataFrame(list(zip(distance_list, time_list)), columns=['distance', 'time'])
    data_to_save.to_csv(f'{path}/sim_result_{i}.csv')
```

### The result file
The result files are csv file that records all the displacement in nanometer of electron at time in second
