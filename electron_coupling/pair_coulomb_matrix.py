"""
This module provide method to generate coulomb matrix for a pair of molecule
"""

from __future__ import annotations
from electron_coupling import coulomb_matrix
import numpy as np

class pair_coulomb_matrix(coulomb_matrix.CoulombMatrix):
    """
    Using this class means that the coulomb matrix is generated from a pair of molecule
    """
    def __init__(self, x: list[float], y: list[float], z: list[float], numerator_matrix: np.ndarray) -> None:
        super().__init__(x, y, z, numerator_matrix)

    def gen_inter_matrix(self)->np.ndarray:
        """
        this method return only the inter molecular element of the coulomb matrix (upper right quadrant)
        to get lower left quadrant, make a tranpost of the returned matrix

        return: np.ndarray
        """
            
        matrix_length = len(self.x)

        full_matrix = self.gen_full_coulomb_matrix()
        
        upper_half = full_matrix[:int(matrix_length/2)]

        upper_right = np.array([i[int(matrix_length/2):] for i in upper_half])
        
        return upper_right

    def conjugate(self)->pair_coulomb_matrix:

        matrix_length = len(self.x_list)
        cutting_index = int(matrix_length/2)

        x = self.x_list[cutting_index:] + self.x_list[:cutting_index]
        y = self.y_list[cutting_index:] + self.y_list[:cutting_index]
        z = self.z_list[cutting_index:] + self.z_list[:cutting_index]
        
        return pair_coulomb_matrix(x, y, z, self.numerator_matrix)
