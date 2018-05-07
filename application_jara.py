from classes.Protein import Protein
from classes.algorithm.Randomizer import Randomizer
from classes.algorithm.DepthFirst import DepthFirst


def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    # TODO: Change these numbers per protein / run !!
    number = 1

    iterations = 10000

    # run random algorithm
    protein = Protein(number)
    randomAlgorithm = Randomizer(protein, iterations, 1)
    randomAlgorithm.runFastRandomizer()

    # run depth first algorithm
    # protein = Protein(number)
    # depthFirstAlgorithm = DepthFirst(protein)
    # depthFirstAlgorithm.runFastDepthFirst()

    # print solutions
    randomAlgorithm.printBest()

    # depthFirstAlgorithm.printBest()

if __name__ == "__main__":
    main()



