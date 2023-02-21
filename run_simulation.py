"""This module run the simulation"""


from simulation import simulation
import pandas as pd
from simulation import DBT3
test_sim = simulation.Simulation('Frames/DBT3-static.gro', DBT3.DBT3, memory_saving=True)

for i in range(1000):
    distance_list, time_list = test_sim.run()
    data_to_save = pd.DataFrame(list(zip(distance_list, time_list)), columns=['distance', 'time'])
    data_to_save.to_csv(f'result_dbt3/sim_result_{i}.csv')

