"""This script make the cache for the dynamic simulation"""

from decouple import config
from simulation import dynamic_simulation
from simulation import simulation
import  multiprocessing as mp
import os

cache_path:str = config('CACHE_PATH')
dynamic_frame_path:str = config('DYNAMIC_FRAME_PATH')

num_workers = mp.cpu_count()
pool = mp.Pool(num_workers)

for file in os.listdir(dynamic_frame_path):
    
    print(file)
    #simulation.extract_metadata(f'{dynamic_frame_path}{file}', cache_path)
    pool.apply_async(simulation.extract_metadata, args=(f'{dynamic_frame_path}{file}', cache_path))

pool.close()
pool.join()
    

