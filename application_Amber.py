from classes.Protein import Protein
from classes.algorithm.Randomizer import Randomizer
from classes.algorithm.DepthFirst import DepthFirst


def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    # TODO: Change these numbers per protein / run !!
    number = 1
    iterations = 2
    dimensions = 2
    write = "OFF"

    # run random algorithm
    protein = Protein(number, dimensions)
    randomAlgorithm = Randomizer(protein, iterations, write)
    randomAlgorithm.runRandomizer()

    # run depth first algorithm
    # protein = Protein(number)
    # depthFirstAlgorithm = DepthFirst(protein, dimensions)
    # depthFirstAlgorithm.runDepthFirst()

    # print solutions
    randomAlgorithm.printBest()
    #depthFirstAlgorithm.printBest()

if __name__ == "__main__":
    main()

# from classes.Protein import Protein
# from classes.algorithm.HillClimber import HillClimber
# from classes.algorithm.Randomizer import Randomizer
# from StabilityAnalyzer import StabilityAnalyzer
#
#
# def main():
#     """ Implements random algorithms in order to most efficiently fold a protein """
#
#     # protein number
#     number = 1
#     #
#     # # initialize values for randomizer
#     # iterationsRand = 10000
#     # writeOptionsRand = 1
#     # dimensions = 3
#     #
#     # # initialize values for hill climber
#     # # iterationsHill = 1000000
#     # # writeOptionsHill = 1
#     #
#     # create protein
#     #
#     # # run randomizer and return best protein
#     # randomAlgorithm = Randomizer(protein, iterationsRand, writeOptionsRand, dimensions)
#     # randomAlgorithm.runRandomizer()
#     # bestPattern = randomAlgorithm.bestPattern
#     # print(bestPattern)
#     # protein.visualize3D("Test", dimensions)
#
#     # run random algorithm
#     # bestHillClimber = HillClimber(protein, bestPattern, iterationsHill, writeOptionsHill)
#     # bestHillClimber.runHillClimber()
#     # bestHillClimber.printBestHill()
#
#     # analyze stability scores
#     # filename = "hillclimber3"
#     # StabilityAnalyzer(filename)
#
#     protein = Protein(number)
#     foldingPattern = ['0', '+Y', '-X', '-Y', '+Z', '-Y', '-Z', '+X', '-Y']
#
#     protein.fold(foldingPattern, 3)
#
# if __name__ == "__main__":
#     main()