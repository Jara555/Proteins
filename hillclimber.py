from classes.Protein import Protein
from classes.algorithm.Randomizer import Randomizer
from classes.algorithm.DepthFirst import DepthFirst
from classes.algorithm.HillClimber import HillClimber
import time


def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    # create protein of class Protein
    protein = Protein(2)

    # run random algorithm
    randomAlgorithm = Randomizer(protein, 1000, 0)

    startRand = time.time()
    randomAlgorithm.runRandomizer()
    endRand = time.time()
    elapsedRand = endRand - startRand

    randomAlgorithm.printBestRandom()

    print('Random Time: ' + str(elapsedRand))

    # run depth first algorithm

    depthFirstAlgorithm = DepthFirst(protein)

    startDepth = time.time()
    depthFirstAlgorithm.runDepthFirst()
    endDepth = time.time()
    elapsedDepth = endDepth - startDepth

    depthFirstAlgorithm.printBestDepth()

    print('Depth First Time: ' + str(elapsedDepth))

    bestHillClimber = HillClimber(protein)
    bestHillClimber.printBestHill()


if __name__ == "__main__":
    main()