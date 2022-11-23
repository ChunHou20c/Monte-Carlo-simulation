"""This is the graph that will be used to track relationship between molecules in the simulation"""
from data_structure.graph import Graph
from data_structure.cube import Cube
from typing import Any
from simulation import DBT1, algorithm, molecule_relation

class Molecule_graph(Graph):
    """This class stores the molecules as node and the relation as weight"""

    def __init__(self):
        super().__init__()

    def add_edge(self, frm:Any, to: Any, relation:molecule_relation.Relation):
        """This is the algorithm that will use to check the relation between molecules and make the relation"""


        if frm not in self.vert_dict.keys():
            self.add_vertex(frm)

        if to not in self.vert_dict.keys():
            self.add_vertex(to)

        self.vert_dict[frm].add_neighbor(to, relation)
        self.vert_dict[to].add_neighbor(frm, relation.conjugate())

def build_graph():
    """This is the algorithm that is used to build the molecular graph"""
