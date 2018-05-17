from classes.Protein import Protein
from classes.algorithm.BranchNBound import BranchNBound
from classes.algorithm.DepthFirst import DepthFirst
from classes.algorithm.HillClimber import HillClimber
from classes.algorithm.Randomizer import Randomizer


def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    # TODO: Change these numbers per protein run !!
    proteinNumber = 4
    dimensions = 3
    writeCsv = "OFF"
    maxIterations = 10000

    # run random algorithm
    protein = Protein(proteinNumber, dimensions)
    randomAlgorithm = Randomizer(protein, writeCsv, maxIterations=100)
    randomAlgorithm.runAlgorithm()

    # run depth first algorithm
    # protein = Protein(proteinNumber, dimensions)
    # depthFirstAlgorithm = DepthFirst(protein, writeCsv, maxIterations=None)
    # depthFirstAlgorithm.runAlgorithm()

    # run branch n bound algorithm
    # protein = Protein(proteinNumber, dimensions)
    # branchNBoundAlgorithm = BranchNBound(protein, writeCsv, maxIterations=None)
    # branchNBoundAlgorithm.runAlgorithm()

    # TODO: Use a (random?) pattern as start of hillclimber
    startPattern = randomAlgorithm.bestPattern

    # # run hill climber algorithm
    protein = Protein(proteinNumber, dimensions)
    hillClimberAlgorithm = HillClimber(protein, writeCsv, maxIterations, startPattern)
    hillClimberAlgorithm.runAlgorithm()

    # print solutions
    randomAlgorithm.printBest()
    # depthFirstAlgorithm.printBest()
    # branchNBoundAlgorithm.printBest()
    hillClimberAlgorithm.printBest()

if __name__ == "__main__":
    main()