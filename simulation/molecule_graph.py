"""This is the graph that will be used to track relationship between molecules in the simulation"""
from typing import Any, Type
from data_structure.graph import Graph, Vertex
from simulation import DBT, algorithm, molecule_relation

class Molecule_vertex(Vertex):
    """This class stores the id of the molecule(name), the molecule itself, and the neighbour of the molecule in a dict
    the class is used for this project only"""

    def __init__(self, node:DBT.DBT) -> None:
        
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

    def __init__(self, molecule:Type[DBT.DBT]):
        super().__init__()
        self.molecule_length = molecule.length

    def add_vertex(self, node:DBT.DBT):
        """method to add a new vertex to the graph"""

        self.num_vertex +=1
        
        new_vertex = Molecule_vertex(node)
        self.vert_dict[new_vertex.get_id()] = new_vertex

    def add_edge(self, frm:Any, to:Any, relation:molecule_relation.Relation):
        """
        This method add the relation vertex to the graph, as this graph use molecule as input,
        the name of the molecules and the value of the keys need to be handled.

        ** this class object does not store the molecule object in its dictionary, the molecular data is store in the vertex

        parameters:

        frm - molecule
        to - molecule
        relation - the relation between the two molecule (the relation should be handle outside before adding the edge)
        """

        if frm.__name__ not in self.vert_dict.keys():
            self.add_vertex(frm)

        if to.__name__ not in self.vert_dict.keys():
            self.add_vertex(to)

        self.vert_dict[frm.__name__].add_neighbor(to.__name__, relation)
        self.vert_dict[to.__name__].add_neighbor(frm.__name__, relation.conjugate())

