from classes.Grid import Grid
from classes.Protein import Protein
from functions.Ufold import Ufold
from functions.foldprotein import foldprotein
from functions.randomizer import randomizer
from functions.visualize import visualize
from functions.stability import stability

import matplotlib.pyplot as plt


def main():
    """ Implements all algorithms in order to most efficiently fold a protein """

    # open protein text file
    with open('data/protein1.txt', 'r') as file:
        protein_string = file.read()

    # print protein name
    print()
    print(protein_string)
    print()

    # create protein of class Protein
    protein = Protein(protein_string)

    # Ufold: get folding pattern for folding type Ufold and lay over protein
    Ufolding_pattern = Ufold(protein_string, 0)
    Ufolded_protein = foldprotein(protein.protein_list, Ufolding_pattern)

    # visualize Ufold and print stability score
    visualize(Ufolded_protein)
    print(stability(Ufolded_protein))

    # Random fold: randomize and visualize protein folding patterns
    for i in range(1):
        # get random folding pattern
        Rfolding_pattern = randomizer(protein_string)
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
