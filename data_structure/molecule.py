"""this module contain the general data structure molecule and atom"""

class atom:
    """atom have x, y, z coordinate and charge"""
    
    def __init__(self, x: float, y: float, z: float, charge: float) -> None:

        self.x = x
        self.y = y
        self.z = z
        self.charge = charge
        self.__name__ = 'atom'

    def __str__(self) -> str:
        
        return self.__name__

class molecule:
    """a molecule consist of a list of atoms"""

    def __init__(self, atoms: list[atom]) -> None:
        
        self.atoms = atoms
        self.__name__ = 'molecule'

    def __str__(self) -> str:
        
        return self.__name__

def atomic_distance(a1:atom, a2:atom)->float:
    """This function return the distance between 2 atom"""

    x1, y1, z1 = a1.x, a1.y, a1.z
    x2, y2, z2 = a2.x, a2.y, a2.z

    distance = ((x1-x2)**2 + (y1-y2)**2 + (z1 - z2)**2)**0.5

    return distance
