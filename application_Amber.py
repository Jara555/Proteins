from classes.Protein import Protein
from classes.algorithm.Randomizer import Randomizer


def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    # TODO: Change these numbers per protein / run !!
    number = 3
    iterations = 100000
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
    #
    # # print solutions
    randomAlgorithm.printBest()
    # depthFirstAlgorithm.printBest()
    # branchNBoundAlgorithm.printBest()

    #analyzer = Experiment("random1.2", "random1.3", "run2", "run1")
    #analyzer.stability()

if __name__ == "__main__":
    main()