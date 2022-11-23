"""this module consist of the general undirectional graph data structure"""

from typing import Any, Optional

class Vertex:
    """Vertex in the graph data structure. This is a general vertex for the structure"""
    
    hint = 'vertex'

    def __init__(self, node:Any) -> None:
        """the node is a molecule, the molecule name will be used as the key to access the vertex
        adjacent is a dictionary that stores the keys of the adjacent molecule"""

        self.id = node
        self.adjacent = {}
    
    def __str__(self):

        return 'Vertex'

    def add_neighbor(self, neighbor:str, weight:Any) -> None:
        """Weight can be anything"""

        self.adjacent[neighbor] = weight # this weight might change later or will do calculation of weight first

    def get_connections(self):
        """getter method to get the neighbor keys"""

        return self.adjacent.keys()

    def get_id(self) -> str:
        """getter method to get the molecule's key"""

        return self.id

    def get_weight(self, neighbor:str) -> float:
        """getter method to get the weight connect to the selected neighbor
        currently doesn't handle key doesn't exist error"""

        return self.adjacent[neighbor]

class Graph:
    """The data structure of the undirectional weighted graph, this is a general graph"""

    def __init__(self):

        self.vert_dict = {}
        self.num_vertex = 0
        
        #at the beginning there is no node and the vertex dictionary is empty

    def __iter__(self):

        return iter(self.vert_dict.values())

    def add_vertex(self, node:Any):
        """method to add a new vertex to the graph"""

        self.num_vertex +=1
        
        new_vertex = Vertex(node)
        self.vert_dict[new_vertex.get_id()] = new_vertex

        return new_vertex

    def get_vertex(self, key)-> Optional[Vertex]:
        """getter method to get the vertex base on the key"""
        
        if key in self.vert_dict:
            
            return self.vert_dict[key]

        else:

            return None

    def add_edge(self, frm:Any, to:Any, weight:Any):
        """method to add edge (connection) to two vertex, overide this when use to modify the behaviour"""

        if frm not in self.vert_dict.keys():
            self.add_vertex(frm)

        if to not in self.vert_dict.keys():
            self.add_vertex(to)

        self.vert_dict[frm].add_neighbor(to, weight)
        self.vert_dict[to].add_neighbor(frm, weight)

    def get_vertices(self):
        """getter method to get all the keys in the graph"""

        return self.vert_dict.keys()
