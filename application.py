from new_classes.Protein import Protein
from new_classes.AminoAcid import AminoAcid

from classes.Grid import Grid
#from classes.Protein import Protein
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
    with open('data/protein1.txt', 'r') as file:
        protein_string = file.read()

    # create protein of class Protein
    protein = Protein(protein_string)

    # make all combinations
    combinations = Protein.makeCombinations(protein)
    print("combinations are: " + str(combinations))

    1000*randomizer

    # # Ufold: get folding pattern for folding type Ufold and lay over protein
    # Ufolding_pattern = Ufold(protein_string, 0)
    # print(Ufolding_pattern)
    # Ufolded_protein = foldprotein(protein.protein_list, Ufolding_pattern)
    #
    # # visualize Ufold and print stability score
    # visualize(Ufolded_protein, 'hi')
    # print(stability(Ufolded_protein))
    # print()

    # test for multipleUfold
    # multipleUfold(protein_string)
    # folding_patterns = multipleUfold(protein_string)
    # folded_protein0 = foldprotein(protein.protein_list, folding_patterns[0])
    # visualize(folded_protein0, '0')
    # stability(folded_protein0)
    # folded_protein1 = foldprotein(protein.protein_list, folding_patterns[1])
    # visualize(folded_protein1, '1')
    # stability(folded_protein1)
    # folded_protein2 = foldprotein(protein.protein_list, folding_patterns[2])
    # visualize(folded_protein2, '2')
    # stability(folded_protein2)
    # folded_protein3 = foldprotein(protein.protein_list, folding_patterns[3])
    # visualize(folded_protein3, '3')
    # stability(folded_protein3)
    #
    # # Random fold: randomize and visualize protein folding patterns
    # for i in range(1):
    #     # get random folding pattern
    #     Rfolding_pattern = randomizer(protein_string)
    #     print(Rfolding_pattern)
    #     Rfolded_protein = foldprotein(protein.protein_list, Rfolding_pattern)
    #
    #     # visualize protein
    #     visualize(Rfolded_protein, "folding {i}", i)
    #
    #     # print stability
    #     print(stability(Rfolded_protein))
    #
    #     # TODO: update figures instead of open new figures
    #     plt.show(block=False)
    #     plt.pause(0.5)
    #     plt.close('all')




if __name__ == "__main__":
    main()
