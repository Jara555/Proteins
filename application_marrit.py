from classes.Protein import Protein
from classes.algorithm.BranchNBound import BranchNBound
from classes.algorithm.Randomizer import Randomizer
from classes.algorithm.DepthFirst import DepthFirst

def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    # TODO: Change these numbers per protein / run !!
    number = 2
    iterations = 1000
    dimensions = 2
    writeCsv = "OFF"

    # TODO: Optional parameter for branch n bound or depth first to limit the amount of iterations
    maxIterations = 1000000

    protein = Protein(number, dimensions)
    # protein.findbonds()
    # protein.prune(7, -2)

    # run random algorithm
    randomAlgorithm = Randomizer(protein, writeCsv, iterations)
    randomAlgorithm.runAlgorithm()

    # run depth first algorithm
    depthFirstAlgorithm = DepthFirst(protein, writeCsv)
    depthFirstAlgorithm.runAlgorithm()

    # run branch n bound algorithm
    branchNBoundAlgorithm = BranchNBound(protein, writeCsv)
    branchNBoundAlgorithm.runAlgorithm()

    # print solutions
    randomAlgorithm.printBest()
    depthFirstAlgorithm.printBest()
    branchNBoundAlgorithm.printBest()


if __name__ == "__main__":
    main()