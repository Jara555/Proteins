from classes.Grid import Grid
from classes.Protein import Protein
from functions.Ufold import Ufold
from functions.foldprotein import foldprotein
from functions.multipleUfold import multipleUfold
from functions.visualize import visualize
from functions.stability import stability

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
    folding_pattern = Ufold(protein_string,0)

    # fold protein based on folding pattern
    folded_protein = foldprotein(protein.protein_list, folding_pattern)

    # initiate new grid and place folded protein in it
    grid2 = Grid(len(protein_string))
    grid2.placeProtein(folded_protein)

    print()
    grid2.printGrid()
    print()

    multipleUfold(protein_string)

    # make new folding pattern based on multipleUfold
    patterns = multipleUfold(protein_string)
    folding_pattern2=patterns[0]
    folding_pattern3=patterns[1]
    folded_protein2 = foldprotein(protein.protein_list, folding_pattern2)
    folded_protein3 = foldprotein(protein.protein_list, folding_pattern3)

    print("stability protein2 is: " + str(stability(folded_protein2)))


    # initiate new grid and place folded protein in it
    grid3 = Grid(len(protein_string))
    grid3.placeProtein(folded_protein2)
    grid4 = Grid(len(protein_string))
    grid4.placeProtein(folded_protein3)

    print()
    grid3.printGrid()
    print()

    print()
    grid4.printGrid()
    print()

    visualize(folded_protein)
    visualize(folded_protein2)
    visualize(folded_protein3)



if __name__ == "__main__":
    main()
