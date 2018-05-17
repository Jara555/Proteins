from classes.Protein import Protein
from classes.algorithm.BranchNBound import BranchNBound
from classes.algorithm.Randomizer import Randomizer
from classes.algorithm.DepthFirst import DepthFirst
#from classes.algorithm.HillClimber import HillClimber
from Experiment import Experiment


def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    # TODO: Change these numbers per protein / run !!
    number = 6
    iterations = 10000
    dimensions = 2
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

    # run hill climber algorith
    # protein = Protein(number, dimensions)
    # hillclimberAlgirhtm = HillClimber(protein, writeCsv)
    # hillclimberAlgirhtm.runHillClimber()
    #
    # # print solutions
    randomAlgorithm.printBest()
    # depthFirstAlgorithm.printBest()
    # branchNBoundAlgorithm.printBest()
    # hillclimberAlgirhtm.printBest()

    #analyzer = Experiment("random1.2", "random1.3", "run2", "run1")
    #analyzer.stability()

if __name__ == "__main__":
    main()