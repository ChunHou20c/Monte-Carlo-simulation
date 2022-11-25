"""
This module is the simulation of electron transfer
this simulation object take a gro file and a config file as input
The config file should contain the parameters that will be used in the simulation,
including the stopping condition, model to use and cut off distance for molecular pair
"""

from simulation import constructors, algorithm
from simulation import DBT1
from simulation.molecule_graph import Molecule_graph as Graph
from simulation.molecule_graph import Molecule_vertex
from itertools import combinations

import random
import numpy as np
from simulation.molecule_relation import Relation

class Simulation:
    """This class define how the simulation should be run from a single frame"""

    def __init__(self, gro_file:str) -> None:
        """gro_file - the file path of the gromac file"""
        
        self.__file__ = gro_file
        self.graph, self.box_width = build_graph(gro_file, 56)
        self.total_jump = 10000
        self.initial_box = np.array([0,0,0])
        self.current_box = np.array([0,0,0])

        print(self.box_width)

    def run(self):
        """this method run the simulation for electron jumps in the periodic space"""
        
        #start by choosing a random molecule - vertex

        initial_key = random.choice([ i for i in self.graph.get_vertices()])
        
        initial_vertex = self.graph.get_vertex(initial_key)
        print(initial_vertex.id)

        current_vertex = initial_vertex
        
        def jump(vertex:Molecule_vertex):
            """This function perform a jump and return the next vertex"""

            #choose a neigbour
            neigbours = vertex.get_connections()
            key = random.choice([i for i in neigbours])

            relation:Relation = vertex.get_weight(key)

            #this part is reserved for the calculation of the jumping time

            #for now we only care about the box tracking

            return key, np.array(relation.translation)
        
        for i in range(self.total_jump):

            new_key, translation = jump(current_vertex)
            
            current_vertex = self.graph.get_vertex(new_key)

            self.current_box += translation

        print(f'final box = {self.current_box}')

        print(type(initial_vertex))
        print(type(current_vertex))
        
        initial_molecule = initial_vertex.molecule
        final_molecule = current_vertex.molecule

        Coord1 = initial_molecule.N_S1_S2_coordinates(self.box_width, tuple(self.initial_box))
        Coord2 = final_molecule.N_S1_S2_coordinates(self.box_width, tuple(self.current_box))

        distance = DBT1.DBT1_distance(*Coord1, *Coord2)

        print(f'distance travelled = {distance}')



def build_graph(file_path:str, chunk_size:int)->Graph:
    """This method build the graph for the simulation, (for internal use only, not to be called in runtime)"""

    file_data = constructors._Gro_file_parser(file_path, chunk_size)
    molecules = file_data['list_of_molecules']
    boundary_data = file_data['boundary']
    graph = Graph()
    print(type(graph))


    for index,molecule in enumerate(molecules):
        molecules[index] = algorithm.complete_molecule(molecule, boundary_data)

    for m1, m2 in combinations(molecules, 2):

        relation = algorithm.molecular_pair_relation(m1, m2, 1.2, boundary_data)
        if (relation is not None):

            graph.add_edge(m1,m2,relation)

    return graph, boundary_data




