import numpy as np
# test
# the Protein class has AminoAcid as inner class and reflects a list of aminoacids
class Protein(object):

    # input should be a string reflecting the protein (e.g. "HPHHHPHH")
    def __init__(self, protein_string):

        self.protein_string = protein_string

        # create empty protein list to be filled with aminoacids
        self.protein_list = []

        # determine length of protein and set starting coordinates at 0
        length_protein = len(self.protein_string)
        x, y = 0, 0

        # loop over chars/aminoacids in protein string
        for aa_index in range(length_protein):

            # create aminoacid of class AminoAcid and append to protein list
            aminoacid = self.AminoAcid(self.protein_string[aa_index], x, y)
            self.protein_list.append(aminoacid)

            # increment y coordinate with 1 for next aminoacid
            y += 1

    # the AminoAcid class sets up the properties of each aminoacid
    class AminoAcid(object):

        # uses a char (H/P) and coordinates as input
        def __init__(self, type, x, y):

            # set type and type specific code
            self.type = type
            if self.type == 'H':
                self.code = 1
            elif self.type == 'P':
                self.code = 2

            # set coordinates
            self.x = x
            self.y = y

            # link to previous and next aminoacid
            if self.y != 0:
                self.prev = y - 1
            else:
                self.prev = None

            self.nex = y + 1

    # print method
    def __str__(self):
        return self.protein_string

# the Grid class is used to initialize a 2D grid and to lay the protein over the grid
class Grid(object):

    # initialize grid based on size
    def __init__(self, size):

        self.size = size
        self.grid = np.zeros((self.size, self.size), dtype=int)

    # method to place protein over grid
    def placeProtein(self, protein):

        # loop over aminoacids in protein and put its code at aminoacids coordinates
        for aminoacid in protein:
            self.grid[aminoacid.x][aminoacid.y] = aminoacid.code

    # method to print grid
    def printGrid(self):
        print(self.grid)

def main():

    # define protein as string
    protein_string = "HPHHHPHH"

    # create protein of class Protein
    protein = Protein(protein_string)

    # create grid
    grid = Grid(len(protein_string))

    # place protein in grind
    grid.placeProtein(protein.protein_list)

    # print protein and grid
    print()
    print(protein)
    print()
    grid.printGrid()
    print()

if __name__ == "__main__":
    main()
