from classes.Protein import Protein
from algorithms.BranchNBound import BranchNBound
from algorithms.DepthFirst import DepthFirst
from algorithms.HillClimber import HillClimber
from algorithms.Randomizer import Randomizer
from Experiment import Experiment

def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    # TODO: Change these numbers per protein run !!
    # proteinNumber = 3
    # dimensions = 3
    # writeCsv = "OFF"
    # maxIterations = 10000
    #
    # # run random algorithms
    # protein = Protein(proteinNumber, dimensions)
    # randomAlgorithm = Randomizer(protein, writeCsv, maxIterations)
    # randomAlgorithm.runAlgorithm()

    # # run depth first algorithms
    # protein = Protein(proteinNumber, dimensions)
    # depthFirstAlgorithm = DepthFirst(protein, writeCsv, maxIterations=None)
    # depthFirstAlgorithm.runAlgorithm()

    # # run branch n bound algorithms
    # protein = Protein(proteinNumber, dimensions)
    # branchNBoundAlgorithm = BranchNBound(protein, writeCsv, maxIterations=None)
    # branchNBoundAlgorithm.runAlgorithm()

    # TODO: Use a (random?) pattern as start of hillclimber
    # startPattern = randomAlgorithm.bestPattern
    #
    # # run hill climber algorithms
    # protein = Protein(proteinNumber, dimensions)
    # hillClimberAlgorithm = HillClimber(protein, writeCsv, maxIterations, startPattern)
    # hillClimberAlgorithm.runAlgorithm()

    # print solutions
    # randomAlgorithm.printBest()
    # depthFirstAlgorithm.printBest()
    # branchNBoundAlgorithm.printBest()
    # hillClimberAlgorithm.printBest()

    protein = 3
    dimensions = 3
    algorithm = "R"

    # Experiment
    experiment = Experiment(protein, dimensions, algorithm)
    experiment.setFile()
    experiment.readFile()
    experiment.visualize()

if __name__ == "__main__":
    main()