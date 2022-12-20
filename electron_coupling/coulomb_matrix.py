import numpy as np

class CoulombMatrix:
    """This class is used to construct the full coulomb matrix object"""

    def __init__(self, x:list[float], y:list[float], z:list[float], numerator_matrix:np.ndarray) -> None:
        """Constructor of the coulomb matrix object, take x, y, and z coordinate matrix to generate the coulomb matrix
        Please make sure the size of the numerator matrix is the same as the requested coulomb matrix size
        the expected x, y, z inputs are a list of x, y, z coordinates"""

        self.numerator_matrix = numerator_matrix
        
        matrix_length = len(x)
        #generate matrix of xyz coordinate
        x_matrix = np.tile(x,(matrix_length,1))
        y_matrix = np.tile(y,(matrix_length,1))
        z_matrix = np.tile(z,(matrix_length,1))

        self.x_list = x
        self.y_list = y
        self.z_list = z

        self.x = x_matrix
        self.y = y_matrix
        self.z = z_matrix

    def gen_full_coulomb_matrix(self):
        """This function generate the numpy array with the calculated coulomb matrix"""
        
        distance_matrix = distance_vect(self.x, self.x.T, self.y, self.y.T, self.z, self.z.T)
        #distance matrix is also the denominator matrix
        
        return self.numerator_matrix/distance_matrix

    def gen_denominator_matrix(self):

        distance_matrix = distance_vect(self.x, self.x.T, self.y, self.y.T, self.z, self.z.T)
        return distance_matrix
    
def dist(x1:float,x2:float,y1:float,y2:float,z1:float,z2:float)->float:
    """This method calculate the distance element in the coulomb matrix"""
    
    if(x1==x2 and y1==y2 and z1==z2):
        return 1

    return ((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**0.5

distance_vect = np.vectorize(dist, otypes=[float], cache = True)
        
def charge_product(arr1:str,arr2:str,is_diagonal:bool=False)->int:
    """function to convert char to charge value"""
    
    charge = {'H':1,'C':6,'N':7,'O':8,'S':16}
    if(not is_diagonal):
        return charge[arr1]*charge[arr2]
    
    return 0.5*(charge[arr1]**2.4)

charge_vect = np.vectorize(charge_product,otypes=[float],cache=True)
"""vectorized charge function"""

def gen_numerator_matrix(path_to_molecule_data:str)->np.ndarray:
    """This function generate numerator matrix from molecule data"""

    with open(path_to_molecule_data, 'r') as f:

        content = f.readlines()

        matrix_size = len(content)

        element_list = [i.replace('\n', '') for i in content]

        element_matrix = np.tile(element_list,(matrix_size,1))
        
        numerator_matrix = charge_vect(element_matrix, element_matrix.T, False)

        diagonal_element = charge_vect(element_list, element_list, True)

        np.fill_diagonal(numerator_matrix, diagonal_element)

        return numerator_matrix
