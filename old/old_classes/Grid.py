import numpy as np

# the Grid class is used to initialize a 2D grid and to lay the protein over the grid
class Grid(object):
    """ Initializes a 2D grid to visualize a protein folded in a grid """

    # initialize grid based on size
    def __init__(self, size):
        """ Initializes grid based on input size """

        self.size = size
        self.grid = np.zeros((self.size, self.size), dtype=np.int32)

    # method to place protein over grid
    def placeProtein(self, protein):
        """ Places a protein over the grid based on the aminoacid coordinates """

        # loop over aminoacids in protein and put its code at aminoacids coordinates
        for aminoacid in protein:
            self.grid[aminoacid.x][aminoacid.y] = aminoacid.code

    # method to print grid
    def printGrid(self):
        """ Prints the grid """

        print(self.grid)