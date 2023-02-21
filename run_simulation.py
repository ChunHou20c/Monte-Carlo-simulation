"""This module run the simulation"""


from simulation import simulation
import pandas as pd
from simulation import DBT1
import os

Simulation = simulation.Define_simulation_model('cnn_dbt1')
test_sim = Simulation('Frames/DBT1-static.gro', DBT1.DBT1, memory_saving=True)

path = "result_dbt1"
# Check whether the specified path exists or not
isExist = os.path.exists(path)
if not isExist:

   # Create a new directory because it does not exist
   os.makedirs(path)

for i in range(10):
    distance_list, time_list = test_sim.run()
    data_to_save = pd.DataFrame(list(zip(distance_list, time_list)), columns=['distance', 'time'])
    data_to_save.to_csv(f'{path}/sim_result_{i}.csv')
