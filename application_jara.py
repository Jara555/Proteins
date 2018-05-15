from classes.Protein import Protein
from classes.algorithm.BranchNBound import BranchNBound
from classes.algorithm.Randomizer import Randomizer
from classes.algorithm.DepthFirst import DepthFirst


def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    # TODO: Change these numbers per protein / run !!
    number = 7
    iterations = 10000
    dimensions = 3
    writeCsv = "OFF"

    # run random algorithm
    protein = Protein(number, dimensions)
    randomAlgorithm = Randomizer(protein, writeCsv, iterations)
    randomAlgorithm.runAlgorithm()

    # # run depth first algorithm
    # protein = Protein(number, dimensions)
    # depthFirstAlgorithm = DepthFirst(protein, writeCsv)
    # depthFirstAlgorithm.runAlgorithm()
    #
    # # run branch n bound algorithm
    # protein = Protein(number, dimensions)
    # branchNBoundAlgorithm = BranchNBound(protein, writeCsv)
    # branchNBoundAlgorithm.runAlgorithm()

    # print solutions
    randomAlgorithm.printBest()
    # depthFirstAlgorithm.printBest()
    # branchNBoundAlgorithm.printBest()

if __name__ == "__main__":
    main()