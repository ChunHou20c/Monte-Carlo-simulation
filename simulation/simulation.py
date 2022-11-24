"""
This module is the simulation of electron transfer
this simulation object take a gro file and a config file as input
The config file should contain the parameters that will be used in the simulation,
including the stopping condition, model to use and cut off distance for molecular pair
"""

from simulation import constructors, algorithm
from simulation.molecule_graph import Molecule_graph as Graph
from itertools import combinations


class Simulation:
    """This class define how the simulation should be run from a single frame"""

    def __init__(self, gro_file:str) -> None:
        """gro_file - the file path of the gromac file"""
        
        self.__file__ = gro_file
        self.graph = build_graph(gro_file, 56)



def build_graph(file_path:str, chunk_size:int)->Graph:
    """This method build the graph for the simulation, (for internal use only, not to be called in runtime)"""

    file_data = constructors._Gro_file_parser(file_path, chunk_size)
    molecules = file_data['list_of_molecules']
    boundary_data = file_data['boundary']
    graph = Graph()


    for index,molecule in enumerate(molecules):
        molecules[index] = algorithm.complete_molecule(molecule, boundary_data)

    for m1, m2 in combinations(molecules, 2):

        relation = algorithm.molecular_pair_relation(m1, m2, 1.2, boundary_data)
        if (relation is not None):

            graph.add_edge(m1,m2,relation)

    return graph




