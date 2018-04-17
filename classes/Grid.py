import numpy as np


# the Grid class is used to initialize a 2D grid and to lay the protein over the grid
class Grid(object):

    # initialize grid based on size
    def __init__(self, size):

        self.size = size
        self.grid = np.zeros((self.size, self.size), dtype=np.int32)

    # method to place protein over grid
    def placeProtein(self, protein):

        # loop over aminoacids in protein and put its code at aminoacids coordinates
        for aminoacid in protein:
            self.grid[aminoacid.x][aminoacid.y] = aminoacid.code

    # method to print grid
    def printGrid(self):
        print(self.grid)