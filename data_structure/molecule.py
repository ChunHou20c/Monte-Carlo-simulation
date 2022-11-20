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
