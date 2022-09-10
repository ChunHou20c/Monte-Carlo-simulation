class atom():
    """Molecule is consist of atoms, this class is mainly used for analysing,
            every line of in the molecule object is an atom"""

    def __init__(self,line:str):
        """Constructor of the atom class, define the atomic mass and xyz position of the atom from the line"""

        self.data = line

    def print_data(self):

        print(self.data)


