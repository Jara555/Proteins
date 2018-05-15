from classes.Protein import Protein
from classes.algorithm.Randomizer import Randomizer
from classes.algorithm.DepthFirst import DepthFirst


def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    # TODO: Change these numbers per protein / run !!
    number = 1
    randomIterations = 1000
    dimensions = 3
    writeCsv = "ON"

    # run random algorithm
    protein = Protein(number, dimensions)
    randomAlgorithm = Randomizer(protein, randomIterations, writeCsv)
    randomAlgorithm.runRandomizer()

    # run depth first algorithm
    protein = Protein(number)
    depthFirstAlgorithm = DepthFirst(protein, dimensions)
    depthFirstAlgorithm.runDepthFirst()

    # print solutions
    randomAlgorithm.printBest()
    depthFirstAlgorithm.printBest()

if __name__ == "__main__":
    main()