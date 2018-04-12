
class AminoAcid(object):

    def __init__(self, type, x, y):
        self.type = type
        self.x = x
        self.y = y

        self.prev = y - 1
        self.nex = y + 1

    def __str__(self):
        return (self.type + " -")

class Grid(object):

    def __init__(self, size):
        self.size = size

    def printGrid(self):
        for row in range(self.size):
            for col in range(self.size):
                print(' 0 ', end="")
            print()

def main():

    h1 = AminoAcid('H', 0, 0)
    h4 = AminoAcid('H', 0, 3)
    h2 = AminoAcid('H', 0, 1)
    p3 = AminoAcid('P', 0, 2)

    protein = [h1, h2, p3, h4]

    print()
    for aminoAcid in protein:
        print(aminoAcid, end="")
    print()

    grid = Grid(4)
    grid.printGrid()


if __name__ == "__main__":
    main()