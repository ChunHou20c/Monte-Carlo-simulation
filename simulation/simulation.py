"""
This module is the simulation of electron transfer
this simulation object take a gro file and a config file as input
The config file should contain the parameters that will be used in the simulation,
including the stopping condition, model to use and cut off distance for molecular pair
"""
from math import log
import pickle
import os
from typing import Callable
from electron_coupling import marcus_equation

from simulation import constructors, algorithm
from simulation import DBT1
from simulation.molecule_graph import Molecule_graph as Graph
from simulation.molecule_graph import Molecule_vertex

#for model
from simulation.model import prediction_model
import tensorflow as tf

from itertools import combinations

import random
import numpy as np
from simulation.molecule_relation import Relation

import gc


class Simulation:
    """This class define how the simulation should be run from a single frame"""

    prediction_model = prediction_model('/home/chunhou/Documents/FYP/simulation/model/cnn2d-deployment/')
    #prediction_model = tf.keras.models.load_model('model/ANN1')
    
    def __init__(self, gro_file:str, cache_path:str= './cache', memory_saving:bool=False) -> None:
        """
        gro_file - the file path of the gromac file
        cache_path - directory to cache the data (the file size is large)
        memory_saving - set to true if a lot of memory is required
        """
        
        self.__file__ = gro_file
        self.graph, self.box_width, self.timestamp = extract_metadata(gro_file, cache_path)
        self.total_jump = 100000
        self.initial_box = np.array([0,0,0])
        self.current_box = np.array([0,0,0])
        self.time = 0
        self.electron_coupling_list, self.electron_coupling_key = make_cache_prediction(self.graph, self.prediction_model)
        
        if memory_saving:

            self.delete_coulomb_matrix()
            print('memory saving enabled, coulomb matrix are deleted!')

    def predicted_electron_coupling(self, key1, key2):

        index = self.electron_coupling_key[(key1, key2)]
        return self.electron_coupling_list[index]

    def single_jump(self, key):
        """
        this method do a single jump from a selected key
        """

        current_key = key
        new_key, jumping_time = jump(self.graph, current_key, self.predicted_electron_coupling)
        translation = self.graph.get_vertex(current_key).get_weight(new_key).translation

        initial_molecule = self.graph.get_vertex(current_key).molecule
        final_molecule = self.graph.get_vertex(new_key).molecule

        x0, y0, z0 = initial_molecule.center_coordinate(self.box_width, (0, 0, 0))
        x1, y1, z1 = final_molecule.center_coordinate(self.box_width, translation)
        
        Vector = (x1-x0, y1-y0, z1-z0)
        
        print(f'{new_key=}, {jumping_time=}, {Vector=}')

        return new_key, jumping_time, Vector
    
    def run(self):
        """this method run the simulation for electron jumps in the periodic space"""
        
        #start by choosing a random molecule - vertex

        initial_key = random.choice([ i for i in self.graph.get_vertices()])
        initial_box = np.array([0,0,0])
        
        print('simulation starting from molecule {}'.format(self.graph.get_vertex(initial_key).id))

        current_key = initial_key

        current_time = 0
        current_box = np.array([0,0,0])
        
        displacement_list = []
        time_list = []

        for i in range(self.total_jump):
                
            new_key, jumping_time = jump(self.graph, current_key, self.predicted_electron_coupling)
            
            current_time += jumping_time

            translation = self.graph.get_vertex(current_key).get_weight(new_key).translation

            current_box += translation
            #if np.sum(abs(np.array(translation)))>=1:

            #    print(f'large translation found! {translation}')

            current_key = new_key

#            if np.sum(abs(current_box)>=10):
#                
#                print(abs(current_box))
#                print(abs(current_box)>=10)
#                print(np.sum(abs(current_box)>=10))
#
#                break
#            if current_time >= 2e-10:
#                
#                print('simulation time limit reached')
#                break

            initial_molecule = self.graph.get_vertex(initial_key).molecule
            current_molecule = self.graph.get_vertex(current_key).molecule

            Coord1 = initial_molecule.center_coordinate(self.box_width, tuple(initial_box))
            Coord2 = current_molecule.center_coordinate(self.box_width, tuple(current_box))

            distance = DBT1.cartesian_distance(Coord1[0], Coord2[0], Coord1[1], Coord2[1], Coord1[2], Coord2[2])

            displacement_list.append(distance)
            time_list.append(current_time)

            if current_time >= 1e-9:
            
                print(f'simulation time reached! number of jumps = {i}')
                break

        print(f'final box = {current_box}')
        
        initial_molecule = self.graph.get_vertex(initial_key).molecule
        final_molecule = self.graph.get_vertex(current_key).molecule

        Coord1 = initial_molecule.center_coordinate(self.box_width, tuple(initial_box))
        Coord2 = final_molecule.center_coordinate(self.box_width, tuple(current_box))

        distance = DBT1.cartesian_distance(Coord1[0], Coord2[0], Coord1[1], Coord2[1], Coord1[2], Coord2[2])

        print(f'distance travelled = {distance}')
        print(f'time taken = {current_time}')

        return displacement_list, time_list
    
    def delete_coulomb_matrix(self)->None:
        """This method is used to delete coulomb matrix from the graph so that memory can be freed"""

        for vertex in self.graph:

            for relation in vertex.adjacent.values():

                del relation.coulomb_matrix
        
        gc.collect

def make_cache_prediction(graph:Graph, prediction_model):
    """This function make the cache for the model prediction"""
    
    list_of_coulomb_matrix = []
    keys = {}
    index = 0 #index use to access the predicted value
    for vertex in graph.vert_dict.values():
        
        for key, relation in vertex.adjacent.items():
            
            if (vertex.id, key) not in keys.keys():
                list_of_coulomb_matrix.append(relation.coulomb_matrix)
                keys[(vertex.id, key)] = index
                keys[(key, vertex.id)] = index
                index += 1

    array_of_coulomb_matrix = np.array(list_of_coulomb_matrix)
    predictions = prediction_model.predict(array_of_coulomb_matrix)
    
    return predictions.flatten(), keys

def jump(graph:Graph, key, func:Callable):
    """This function perform a jump and return the next vertex"""

    #find all possible neighbour
    _vertex = graph.get_vertex(key)
    neigbours = _vertex.get_connections()
    #print(neigbours)
    #calculate total rate
    #jump_time = []
    options = []
    rates = []

    for neighbour_key in neigbours:
        
        electron_coupling = func(_vertex.id, neighbour_key)
        reorganization_energy = 0.180
        temperature = 300
        Eij = 0
        
        options.append(neighbour_key)
    
        #random_number = -log(random.uniform(0,1))
        rates.append(marcus_equation.transfer_rate(electron_coupling, reorganization_energy, temperature, Eij))

    new_key, jumping_rate = random_weight_selector(options, rates)
    
    #min_time_index = np.argmin(jump_time)
    #new_key = options[min_time_index]
    random_jumping_rate = np.random.exponential(jumping_rate)
    jumping_time = 1/random_jumping_rate
    #print('{:e}'.format(jumping_rate))

    #print(f'{random_number=}')
    #jumping_time = random_number/jumping_rate
    #this part is reserved for the calculation of the jumping time

    #for now we only care about the box tracking
    #print(new_key)

    return new_key, jumping_time

def random_weight_selector(keys:list[str], _weights:list[float]):
    """A wrapper to random choice of new key"""

    list_of_tuple = [i for i in zip(keys, _weights)]
    choice = random.choices(list_of_tuple, weights=tuple(_weights), k = 1)
    return choice[0]

def extract_metadata(file_path:str, cache_directory:str='./cache')->tuple[Graph,float, float]:
    """This method build the graph for the simulation, (for internal use only, not to be called in runtime)"""
    
    file_name = file_path.split('/')[-1]
    cache_path = f'{cache_directory}/{file_name}'
    if os.path.isfile(cache_path):

        with open(cache_path, 'rb') as f:

            graph, boundary_data, timestamp = pickle.load(f)
    else:
        graph = Graph()
        file_data = constructors._Gro_file_parser(file_path, graph.molecule_length)
        molecules = file_data['list_of_molecules']
        boundary_data = file_data['boundary']
        timestamp = file_data['timestamp']

        for index,molecule in enumerate(molecules):
            molecules[index] = algorithm.complete_molecule(molecule, boundary_data)

        for m1, m2 in combinations(molecules, 2):

            relation = algorithm.molecular_pair_relation(m1, m2, 1.2, boundary_data)
            if (relation is not None):

                graph.add_edge(m1,m2,relation)
        
        with open(cache_path, 'wb') as f:
            
            object_to_save = (graph, boundary_data, timestamp)
            pickle.dump(object_to_save, f)

    return graph, boundary_data, timestamp


