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
        self.x = float(data[-3])
        self.y = float(data[-2])
        self.z = float(data[-1])
    
    def get_coordinates(self)->tuple(float):
        """getter function that return the xyz coordinates of the atom"""

        return (self.x, self.y, self.z)

    def print_data(self):

        print(f'atom = {self.element}')
        print(f'xyz coordinate = {self.x},{self.y},{self.z}')

def distance(a1:atom,a2:atom)->float:
    """atomic distance base on xyz coordinate of atom 1 and atom 2

    return: float"""

    x1, y1, z1 = a1.get_coordinates()
    x2, y2, z2 = a2.get_coordinates()

    distance = ((x1-x2)**2 + (y1-y2)**2 + (z1 - z2)**2)**0.5

    return distance
