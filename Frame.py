"""This module is used to process the frame object into data that can be loaded"""

class frame():
    """The frame class contains the frame number, all the molecules number and also the sides of the box of the frame"""
    def __init__(self,path:str) -> None:
        """The constructor to read all the data from the text file, take path of the file as argument"""

        with open(path) as frame_data:
            
            for lines in frame_data:
                print(lines)
                break

def main():

    frame0 = frame('Frames/DBT1-00')

if(__name__=='__main__'):

    main()
