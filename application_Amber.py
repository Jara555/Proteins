from classes.Protein import Protein
from classes.algorithm.HillClimber import HillClimber
from classes.algorithm.Randomizer import Randomizer
from StabilityAnalyzer import StabilityAnalyzer


def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    # protein number
    number = 3

    # initialize values for randomizer
    iterationsRand = 100
    writeOptionsRand = 1

    # initialize values for hill climber
    iterationsHill = 1000000
    writeOptionsHill = 1

    # create protein
    protein = Protein(number)

    # run randomizer and return best protein
    randomAlgorithm = Randomizer(protein, iterationsRand, writeOptionsRand)
    randomAlgorithm.runFastRandomizer()
    bestPattern = randomAlgorithm.bestPattern

    # run random algorithm
    bestHillClimber = HillClimber(protein, bestPattern, iterationsHill, writeOptionsHill)
    bestHillClimber.runHillClimber()
    bestHillClimber.printBestHill()

    # analyze stability scores
    # filename = "hillclimber3"
    # StabilityAnalyzer(filename)

if __name__ == "__main__":
    main()