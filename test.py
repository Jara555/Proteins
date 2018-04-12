import numpy as np

class AminoAcid(object):

    def __init__(self, type, x, y):

        self.type = type
        if self.type == 'H':
            self.code = 1
        elif self.type == 'P':
            self.code = 2

        self.x = x
        self.y = y

        self.prev = y - 1
        self.nex = y + 1

    def __str__(self):
        return self.type


class Grid(object):

    def __init__(self, size):
        self.size = size

    def initGrid(self):
        self.matrix = np.zeros((self.size, self.size), dtype=int)

    def placeAmino(self, proteins):

        for aminoacid in proteins:
            self.matrix[aminoacid.x][aminoacid.y] = aminoacid.code

    def printGrid(self):
        print(self.matrix)

def main():

    aa1 = AminoAcid('H', 0, 0)
    aa2 = AminoAcid('H', 0, 1)
    aa3 = AminoAcid('P', 0, 2)
    aa4 = AminoAcid('H', 0, 3)
    aa5 = AminoAcid('H', 0, 4)
    aa6 = AminoAcid('H', 0, 5)
    aa7 = AminoAcid('P', 0, 6)
    aa8 = AminoAcid('H', 0, 7)

    protein = [aa1, aa2, aa3, aa4, aa5, aa6, aa7, aa8]

    print()
    for aminoAcid in protein:
        print("{} - ".format(aminoAcid), end="")
    print("\n")

    grid = Grid(len(protein))

    grid.initGrid()

    grid.placeAmino(protein)

    grid.printGrid()

    print()


if __name__ == "__main__":
    main()
