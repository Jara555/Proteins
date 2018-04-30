from classes.Protein import Protein
from classes.algorithm.Randomizer import Randomizer
from classes.algorithm.DepthFirst import DepthFirst


def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    # create protein of class Protein
    protein = Protein(1)

    # run random algorithm
    randomAlgorithm = Randomizer(protein, 1000, 0)
    randomAlgorithm.runRandomizer()
    randomAlgorithm.printBestRandom()

    # run depth first algorithm
    depthFirstAlgorithm = DepthFirst(protein)
    depthFirstAlgorithm.runDepthFirst()
    depthFirstAlgorithm.printBestDepth()


if __name__ == "__main__":
    main()