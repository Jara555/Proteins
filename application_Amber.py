from classes.Protein import Protein
from classes.algorithm.HillClimber import HillClimber
from classes.algorithm.Randomizer import Randomizer
from StabilityAnalyzer import StabilityAnalyzer


def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    # protein number
    number = 4

    # initialize values for randomizer
    iterationsRand = 1000
    writeOptionsRand = 1

    # initialize values for hill climber
    iterationsHill = 100000
    writeOptionsHill = 0

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

    # analyze stabilityscores
    filename = "hillclimber4"
    StabilityAnalyzer(filename)




if __name__ == "__main__":
    main()