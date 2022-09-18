import numpy as np

class CoulombMatrix:
    """This class is used to construct the full coulomb matrix object"""

    def __init__(self, x:np.ndarray, y:np.ndarray, z:np.ndarray, numerator_matrix:np.ndarray, Use_full_CM = True) -> None:
        """Constructor of the coulomb matrix object, take x, y, and z coordinate matrix to generate the coulomb matrix

        the expected x, y, z inputs are a 2 dimensional full array"""
        
        self.x = x
        self.xT = x.T
        self.y = y
        self.yT = y.T
        self.z = z
        self.zT = z.T
        self.numerator_matrix = numerator_matrix
        
        #currently the use_full_CM is not handled
        
    def gen_CM(self):
        """This function generate the numpy array with the calculated coulomb matrix"""
        
        distance_matrix = distance_vect(self.x, self.xT, self.y, self.yT, self.z, self.z)
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
