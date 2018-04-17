from classes.Grid import Grid
from classes.Protein import Protein
from functions.Ufold import Ufold
from functions.foldprotein import foldprotein

def main():
    """ Implements all algorithms in order to most efficiently fold a protein """

    # open protein text file
    with open('data/protein1.txt', 'r') as file:
        protein_string = file.read()

    # create protein of class Protein
    protein = Protein(protein_string)

    # create grid of class Grid
    grid = Grid(len(protein_string))

    # place protein in grid
    grid.placeProtein(protein.protein_list)

    # print initial protein and grid
    print()
    print(protein)
    print()
    grid.printGrid()
    print()

    # get folding pattern for folding type Ufold
    folding_pattern = Ufold(protein_string)

    # fold protein based on folding pattern
    folded_protein = foldprotein(protein.protein_list, folding_pattern)

    # initiate new grid and place folded protein in it
    grid2 = Grid(len(protein_string))
    grid2.placeProtein(folded_protein)

    print()
    grid2.printGrid()
    print()


if __name__ == "__main__":
    main()
