class atom():
    """Molecule is consist of atoms, this class is mainly used for analysing,
            every line of in the molecule object is an atom"""

    def __init__(self,line:str):
        """Constructor of the atom class, define the atomic mass and xyz position of the atom from the line"""

        data = line.split()
        #this line extract the element from the string, currently it is fine to take the first character only as this will be faster
        #and there is no element with 2 character currently
        self.element = data[1][0]
        
        #as the index of the atom will get mix up with the element at high number, it is ignored here
        self.x = data[-3]
        self.y = data[-2]
        self.z = data[-1]

    def print_data(self):

        print(f'atom = {self.element}')
        print(f'xyz coordinate = {self.x},{self.y},{self.z}')


