"""This is the graph that will be used to track relationship between molecules in the simulation"""
from data_structure.graph import Graph, Vertex
from data_structure.cube import Cube
from typing import Any
from simulation import DBT1, algorithm, molecule_relation

class Molecule_vertex(Vertex):
    """This class stores the id of the molecule(name), the molecule itself, and the neighbour of the molecule in a dict
    the class is used for this project only"""

    def __init__(self, node:DBT1.DBT1) -> None:
        
        self.id = node.__name__
        self.molecule = node
        self.adjacent = {}

    def __str__(self):
        
        return f'Vertex : molecule {self.id}'

    def add_neighbor(self, neighbor: str, weight:molecule_relation.Relation) -> None:
        """
        parameters:
            
        neigbour - id of the neighbor
        weight - relation generated from the pair
        """
        
        self.adjacent[neighbor] = weight

class Molecule_graph(Graph):
    """This class stores the molecules as node and the relation as weight"""

    def __init__(self):
        super().__init__()

    def add_edge(self, frm:DBT1.DBT1, to:DBT1.DBT1, relation:molecule_relation.Relation):
        """
        This method add the relation vertex to the graph, as this graph use DBT1 molecule as input,
        the name of the molecules and the value of the keys need to be handled.

        ** this class object does not store the molecule object in its dictionary, the molecular data is store in the vertex

        parameters:

        frm - DBT1 molecule
        to - DBT1 molecule
        relation - the relation between the two molecule (the relation should be handle outside before adding the edge)
        """

        if frm.__name__ not in self.vert_dict.keys():
            self.add_vertex(frm.__name__)

        if to.__name__ not in self.vert_dict.keys():
            self.add_vertex(to.__name__)

        self.vert_dict[frm.__name__].add_neighbor(to.__name__, relation)
        self.vert_dict[to.__name__].add_neighbor(frm.__name__, relation.conjugate())

def build_graph():
    """This is the algorithm that is used to build the molecular graph"""
