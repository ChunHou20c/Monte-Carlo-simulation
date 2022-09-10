from frame_data_processing.atom import atom

class molecule():
    """This is a inner class to store the molecular data"""

    def __init__(self,molecular_data:list[str])->None:
        """Constructor of the molecule object.
            Take file object as argument and separates the molecule into atoms for processing"""

        self.raw_data = molecular_data
        self.atoms = [atom(line) for line in self.raw_data]

    def print_data(self):
        """Function to check the functionality of the molecular class"""

        print(len(self.raw_data))


