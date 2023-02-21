"""This module define the prediction model, any customization of model should be done here"""

from simulation.model import delete_element
import numpy as np

class random_model:
    """This is the model that returns random electronic coupling from training dataset"""

    def __init__(self, data_file:str,*args) -> None:
        
        self.data = np.load(data_file)

    def predict(self, matrix):
        """
        This function randomly return electronic coupling
        """
        
        Size = len(matrix)
        
        return np.random.choice(self.data, size=Size)

class prediction_model:
    """This is a wrapper class for the tensorflow ml model"""
    
    def __init__(self, model_path:str) -> None:
        
        import tensorflow as tf

        self.model = tf.keras.models.load_model(model_path, compile = False)

    def summary(self)->None:
        """a wrapper to the summary of model summary"""

        self.model.summary()
    
    def predict(self, matrix):
        """This method is a wrapper of the original model predict, with some preprocessing and postprocessing"""

        matrix_to_predict = process_matrix(matrix)

        result = self.model.predict(matrix_to_predict)

        return abs(result)

def process_matrix(matrix):
    """This is processing of matrix before feeding into the model"""

    processed_matrix = nomad(matrix)\
            .apply(delete_element.delete_elements)\
            .apply(delete_element.convert_to_inter_matrix)\
            .value

    return processed_matrix

class nomad:
    """
    This class is use to capture the error when converting the matrix
    """

    def __init__(self, value)-> None:
        
        self.value = value

    def apply(self, func):
        
        value = func(self.value)

        return nomad(value)
