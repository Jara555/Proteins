from classes.Grid import Grid
from classes.Protein import Protein
from functions.Ufold import Ufold
from functions.foldprotein import foldprotein
from functions.multipleUfold import multipleUfold
from functions.visualize import visualize
from functions.stability import stability
from functions.randomizer import randomizer
from functions.visualize import visualize
from functions.stability import stability

import matplotlib.pyplot as plt


def main():
    """ Implements all algorithms in order to most efficiently fold a protein """

    # open protein text file
    with open('data/protein2.txt', 'r') as file:
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

    # Ufold: get folding pattern for folding type Ufold and lay over protein
    Ufolding_pattern = Ufold(protein_string, 0)
    Ufolded_protein = foldprotein(protein.protein_list, Ufolding_pattern)

    # visualize Ufold and print stability score
    visualize(Ufolded_protein)
    print(stability(Ufolded_protein))

    # test for multipleUfold
    multipleUfold(protein_string)
    # folding_patterns = multipleUfold(protein_string)
    # folded_protein1 = foldprotein(protein.protein_list, folding_patterns[0])
    # visualize(folded_protein1)
    # folded_protein2 = foldprotein(protein.protein_list, folding_patterns[1])
    #
    # visualize(folded_protein2)

    # Random fold: randomize and visualize protein folding patterns
    for i in range(1):
        # get random folding pattern
        Rfolding_pattern = randomizer(protein_string)
        print(Rfolding_pattern)
        Rfolded_protein = foldprotein(protein.protein_list, Rfolding_pattern)

        # visualize protein
        visualize(Rfolded_protein)

        # print stability
        print(stability(Rfolded_protein))

        # TODO: update figures instead of open new figures
        plt.show(block=False)
        plt.pause(0.5)
        plt.close('all')


if __name__ == "__main__":
    main()
