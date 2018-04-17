from classes.Grid import Grid
from classes.Protein import Protein
from functions.Ufold import Ufold


def main():

    # open protein text file
    with open('data/protein1.txt', 'r') as file:
        protein_string = file.read()

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

    Ufold(protein)

    print()
    print(protein)
    print()
    grid.printGrid()
    print()

if __name__ == "__main__":
    main()
