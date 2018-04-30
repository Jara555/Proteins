from classes.Protein import Protein
from classes.algorithm.Randomizer import Randomizer


def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    # create protein of class Protein
    protein = Protein(1)

    # run random algorithm
    randomAlgorithm = Randomizer(protein, 1000, 0)
    randomAlgorithm.runRandomizer()
    randomAlgorithm.printBest()

if __name__ == "__main__":
    main()