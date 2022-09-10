class atom():
    """Molecule is consist of atoms, this class is mainly used for analysing,
            every line of in the molecule object is an atom"""

    def __init__(self,line:str):
        """Constructor of the atom class, define the atomic mass and xyz position of the atom from the line"""

        data = line.split()
        self.name = data[0]
        self.element = data[1]
        
        #as the index of the atom will get mix up with the element at high number, it is ignored here
        self.x = data[-3]
        self.y = data[-2]
        self.z = data[-1]

    def print_data(self):

        print(f'name = {self.name}')
        print(f'atom = {self.element}')
        print(f'xyz coordinate = {self.x},{self.y},{self.z}')


