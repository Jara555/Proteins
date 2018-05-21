from classes.Protein import Protein
from algorithms import HillClimber
from algorithms.Randomizer import Randomizer


def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    # TODO: Change these numbers per protein run !!
    proteinNumber = 7
    dimensions = 3
    writeCsv = "OFF"
    maxIterations = 10000

    # run random algorithms
    protein = Protein(proteinNumber, dimensions)
    randomAlgorithm = Randomizer(protein, writeCsv, maxIterations=100)
    randomAlgorithm.runAlgorithm()

    # run depth first algorithms
    # protein = Protein(proteinNumber, dimensions)
    # depthFirstAlgorithm = DepthFirst(protein, writeCsv, maxIterations=None)
    # depthFirstAlgorithm.runAlgorithm()

    # run branch n bound algorithms
    # protein = Protein(proteinNumber, dimensions)
    # branchNBoundAlgorithm = BranchNBound(protein, writeCsv, maxIterations=None)
    # branchNBoundAlgorithm.runAlgorithm()

    # TODO: Use a (random?) pattern as start of hillclimber
    startPattern = randomAlgorithm.bestPattern

    # # run hill climber algorithms
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