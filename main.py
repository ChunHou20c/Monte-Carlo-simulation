from frame_data_processing import frame
import  multiprocessing as mp
import numpy as np
from os import listdir

def save_cm(frame:frame.frame, numerator_matrix:np.ndarray)->None:
    """Function that is used to save the npy file at the end of coulomb matrix generation, use mainly for multiprocessing"""

    list_of_cm = frame.gen_cm_list(numerator_matrix)
    file_name = f'{frame.dir_name}'
    np.save(file_name,list_of_cm)

def main():

    #project setup
    working_dir = '/lustre/user/chunhou/FYP_data'
    data_dir = '/home/user/chunhou/time'
    molecule_size=56
    files = listdir(data_dir)
    #setting up multiprocessing
    num_workers = mp.cpu_count()
    pool = mp.Pool(num_workers)
    
    #this part is to prevent recalculation of numerator matrix
    frame0 = frame.frame(f'{data_dir}/{files[0]}', working_dir=working_dir,
            molecule_size=molecule_size, cut_off_distance= 1.2,
            Use_full_CM = False)

    numerator_matrix = frame.gen_numerator_matrix(frame0, False)
    save_cm(frame0, numerator_matrix)
    for file in files:
        print(file)       
        Frame = frame.frame(f'{data_dir}/{file}', working_dir=working_dir,
                molecule_size=molecule_size,
                cut_off_distance=1.2,
                Use_full_CM=False)
        
        pool.apply_async(save_cm, args=(Frame, numerator_matrix))

    pool.close()
    pool.join()

if(__name__=='__main__'):

    main()
