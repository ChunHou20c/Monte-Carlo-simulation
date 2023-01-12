"""
This script delete unwanted element from the coulomb matrix
"""
import numpy as np

def delete_elements(matrix:np.ndarray)->np.ndarray:
    """
    this function convert the original 112x112 array of DBT1 into array of 66x66 by removing unwanted element
    """
    H_index = [0,3,8,10,13,17,19,23,25,28,30,33,34,35,38,40,45,47,49,51,55]
    C_index = [32]
    O_index = [31]

    index_removal = np.array(H_index + C_index + O_index)
    pair_index_removal = np.concatenate((index_removal, index_removal+56))

    reduced_matrix_list = []
    for i in matrix:

        reduced_matrix = np.delete(i, pair_index_removal, axis=0)
        reduced_matrix = np.delete(reduced_matrix, pair_index_removal, axis=1)
        
        reduced_matrix_list.append(reduced_matrix)

    return np.array(reduced_matrix_list)

def convert_to_inter_matrix(matrix:np.ndarray):
    """This function convert the 66x66 array into 33x33 array by removing the intramolecular element in the matrix"""

    inter_feature = []

    for i in matrix:

        reduced_matrix = np.delete(i, [j for j in range(33)], axis=0)
        reduced_matrix = np.delete(reduced_matrix,[k for k in range(33,66)], axis=1)

        inter_feature.append(reduced_matrix)

    feature = np.array(inter_feature)

    return feature
