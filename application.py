from classes.Grid import Grid
from classes.Protein import Protein


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
