import sys
import pandas as pd
import numpy as np

MATRIX_SIZE = 112

print('generating',len(sys.argv)-1,'matrices')
# command line argument, pass the file to generate coulomb matrix

# vectorized function to get atomic charge
def char_to_charge(array):
    charge = {'H':1,'C':6,'N':7,'O':8,'S':16}
    return charge[array[0]]

def gen_CM(position_matrix_row,position_matrix_column,coordx1,coordx2,coordy1,coordy2,coordz1,coordz2,charge1,charge2):
    if(position_matrix_row==position_matrix_column):
        return 0.5*charge1**2.4
    return (charge1*charge2)/((coordx1-coordx2)**2+(coordy1-coordy2)**2+(coordz1-coordz2)**2)**0.5

#vectorize the functions
char_to_charge = np.vectorize(char_to_charge,otypes=[object],cache=False)
gen_CM = np.vectorize(gen_CM,otypes=[float],cache=False)

char_to_charge = np.vectorize(char_to_charge,otypes=[object],cache=False)
pos_matrix_x=np.indices((MATRIX_SIZE,MATRIX_SIZE))[0] 


CM_list = []

for files in sys.argv[1:]:
    data = pd.read_csv(files,delimiter='\\s+',on_bad_lines='skip',names=[i for i in range(9)],skipfooter=1,engine='python')
    data['charge'] = char_to_charge(data[1])

    #create the matrices for vectorized function gen_matrix
    x = np.tile(list(data.iloc[:][3]),(MATRIX_SIZE,1))
    y = np.tile(list(data.iloc[:][4]),(MATRIX_SIZE,1))
    z = np.tile(list(data.iloc[:][5]),(MATRIX_SIZE,1))
    Charge = np.tile(list(data['charge']),(MATRIX_SIZE,1))
    
    CM_list.append(gen_CM(pos_matrix_x, pos_matrix_x.T, x, x.T, y, y.T, z, z.T, Charge, Charge.T))


np.save('input_feature.npy',CM_list)
