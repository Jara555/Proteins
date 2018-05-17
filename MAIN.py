from classes.Protein import Protein
from classes.algorithm.BranchNBound import BranchNBound
from classes.algorithm.Randomizer import Randomizer
from classes.algorithm.DepthFirst import DepthFirst


def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    # TODO: Change these numbers per protein / run !!
    proteinNumber = 1
    dimensions = 2
    writeCsv = "OFF"
    maxIterations = 1000

    # run random algorithm
    protein = Protein(proteinNumber, dimensions)
    randomAlgorithm = Randomizer(protein, writeCsv, maxIterations)
    randomAlgorithm.runAlgorithm()

    # run depth first algorithm
    protein = Protein(proteinNumber, dimensions)
    depthFirstAlgorithm = DepthFirst(protein, writeCsv, maxIterations=None)
    depthFirstAlgorithm.runAlgorithm()

    # run branch n bound algorithm
    protein = Protein(proteinNumber, dimensions)
    branchNBoundAlgorithm = BranchNBound(protein, writeCsv, maxIterations=None)
    branchNBoundAlgorithm.runAlgorithm()

    #print solutions
    randomAlgorithm.printBest()
    depthFirstAlgorithm.printBest()
    branchNBoundAlgorithm.printBest()

if __name__ == "__main__":
    main()