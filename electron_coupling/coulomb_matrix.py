import numpy as np

class CoulombMatrix:
    """This class is used to construct the full coulomb matrix object"""

    def __init__(self, x1:list[float], x2:list[float],
            y1:list[float], y2:list[float],
            z1:list[float], z2:list[float],
            numerator_matrix:np.ndarray,
            molecular_size:int,
            Use_full_CM = True) -> None:
        """Constructor of the coulomb matrix object, take x, y, and z coordinate matrix to generate the coulomb matrix
        Please make sure the size of the numerator matrix is the same as the requested coulomb matrix size
        the expected x, y, z inputs are a list of x, y, z coordinates"""

        self.numerator_matrix = numerator_matrix

        if(Use_full_CM):
            x_list = x1 + x2 
            y_list = y1 + y2
            z_list = z1 + z2

            x = np.tile(x_list,(molecular_size*2,1))
            y = np.tile(y_list,(molecular_size*2,1))
            z = np.tile(z_list,(molecular_size*2,1))

            self.x = x
            self.xT = x.T
            self.y = y
            self.yT = y.T
            self.z = z
            self.zT = z.T

        else:
            
            x1_array = np.tile(x1,(molecular_size,1))
            x2_array = np.tile(x2,(molecular_size,1))
            y1_array = np.tile(y1,(molecular_size,1))
            y2_array = np.tile(y2,(molecular_size,1))
            z1_array = np.tile(z1,(molecular_size,1))
            z2_array = np.tile(z2,(molecular_size,1))

            self.x = x1_array
            self.xT = x2_array.T
            self.y = y1_array
            self.yT = y2_array.T
            self.z = z1_array
            self.zT = z2_array.T

    def gen_CM(self):
        """This function generate the numpy array with the calculated coulomb matrix"""
        
        distance_matrix = distance_vect(self.x, self.xT, self.y, self.yT, self.z, self.zT)
        #distance matrix is also the denominator matrix
        
        return self.numerator_matrix/distance_matrix
    
def dist(x1:float,x2:float,y1:float,y2:float,z1:float,z2:float)->float:
    """This method calculate the distance element in the coulomb matrix"""
    
    if(x1==x2 and y1==y2 and z1==z2):
        return 1

    return ((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**0.5

distance_vect = np.vectorize(dist, otypes=[float], cache = True)
        
def charge(arr1:str,arr2:str,is_diagonal:bool=False)->int:
    """function to convert char to charge value"""
    
    charge = {'H':1,'C':6,'N':7,'O':8,'S':16}
    if(not is_diagonal):
        return charge[arr1]*charge[arr2]
    
    return 0.5*(charge[arr1]**2.4)

charge_vect = np.vectorize(charge,otypes=[float],cache=True)
"""vectorized charge function"""

def gen_cm(distance, charge1, charge2)->float:
    """function to do calculation"""
    if(distance != -1):
        return (charge1*charge2)/distance
    
    return 0.5*charge1**2.4

gen_cm_vect = np.vectorize(gen_cm,otypes=[float],cache=True)
