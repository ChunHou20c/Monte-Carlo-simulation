from frame_data_processing import frame
import  multiprocessing as mp
import numpy as np
from os import listdir

def save_cm(frame:frame.frame, numerator_matrix:np.ndarray, filter:str="close to border")->None:
    """Function that is used to save the npy file at the end of coulomb matrix generation, use mainly for multiprocessing
    filter: str = "close to border" or "cut by boundary" """

    if (filter == "close to border"):

        func = frame.is_close_to_boundary

    elif (filter == "cut by boundary"):

        func = frame.is_cut_by_boundary

    elif (filter == "no filter"):
        
        def _func(*arg):
            
            return False

        func = _func

    else:

        print("invalid function!")
        return

    list_of_cm = frame.gen_cm_list(numerator_matrix, func)
    file_name = f'{frame.dir_name}'
    np.save(file_name,list_of_cm)

def check_molecule(frame:frame.frame, filter:str)->None:
    """Function to check the molecular data"""

    if (filter == "close to border"):

        func = frame.is_close_to_boundary

    elif (filter == "cut by boundary"):

        func = frame.is_cut_by_boundary

    elif (filter == "no filter"):
        
        def _func(*arg):

            return False
            
        func = _func

    else:

        print("invalid function!")
        return
    
    count = 0
    string_to_save = ''
    for molecule in frame.molecules:
        
        if (func(molecule)):

            string_to_save += molecule.print_data()
            count+=1

    print(count)
    print(string_to_save)
    
    with open("tmp.txt", mode = 'w') as f:

        f.write(string_to_save)


def main():

    #project setup
    working_dir = '/home/chunhou/Documents/FYP'
    molecule_size=56
    files = listdir('Frames')
    
    #setting up multiprocessing
    num_workers = mp.cpu_count()
    pool = mp.Pool(num_workers)
    
    #this part is to prevent recalculation of numerator matrix
    frame0 = frame.frame(f'Frames/{files[0]}', working_dir=working_dir,
            molecule_size=molecule_size, cut_off_distance= 1.2,
            Use_full_CM = False)

    numerator_matrix = frame.gen_numerator_matrix(frame0, False)

    for file in files:
        print(file)       
        Frame = frame.frame(f'Frames/{file}', working_dir=working_dir,
                molecule_size=molecule_size,
                cut_off_distance=1.2,
                Use_full_CM=False)
        
        #pool.apply_async(save_cm, args=(Frame, numerator_matrix, "cut by boundary"))
        pool.apply_async(check_molecule, args=(Frame, "cut by boundary"))

    pool.close()
    pool.join()

if(__name__=='__main__'):

    main()
