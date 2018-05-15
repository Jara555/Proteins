from classes.Protein import Protein
from classes.algorithm.HillClimber import HillClimber
from classes.algorithm.BranchNBound import BranchNBound
from classes.algorithm.DepthFirst import DepthFirst
from classes.algorithm.Randomizer import Randomizer
from classes.algorithm.HillClimber import HillClimber
import time

def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    # create protein of class Protein
    protein = Protein(4)
    iterationsRandom = 1000
    iterationsHill = 100000
    dimensions = 3
    writeoptions = 0

    randomizer = Randomizer(protein, iterationsRandom, 1, dimensions)
    randomizer.runRandomizer()

    hillclimber = HillClimber(protein, randomizer.bestPattern, iterationsHill, writeoptions, dimensions)
    hillclimber.runHillClimber()
    hillclimber.printBestHill()

    # foldpattern = ['0', '+Y', '+Z', '-Y', '-Y', '-Y', '-X', '+Y']
    #
    # protein.fold(foldpattern, dimensions)
    # protein.stability(protein.length, dimensions)
    # print(protein.stabilityScore)
    # protein.visualize("protein 3, bestScore = " + str(protein.stabilityScore), dimensions)



    # depthFirst = DepthFirst(protein, dimensions)
    # depthFirst.runDepthFirst()
    # depthFirst.printBest()

    # depthFirst = DepthFirst(protein)
    # depthFirst.runDepthFirst()
    # depthFirst.printBest()

    # randomizer = Randomizer(protein, iterations, 1, dimensions)
    # randomizer.runRandomizer()
    # randomizer.printBest()

    # maxStability = randomizer.maxStability
    #
    # branchNBound = BranchNBound(protein, dimensions)
    # branchNBound.runBranchNBound()
    # branchNBound.printBest()


    #
    # protein.fold(['0', '+Y', '-X', '-X', '-X', '-X', '-X', '-Y', '-Y', '-X', '-Y', '-Y', '-X', '+Y', '+X']


    # protein.stability(protein.length)
    # protein.fold(['0', '+Y', '-X', '-X', '-Y', '+X', '+X', '+X', '+X', '+Y', '+X', '-Y', '+X', '+Y'])
    # protein.visualize("5th aminoacid random '+X'")
    #
    # # startHill = time.time()
    # # bestHillClimber = HillClimber(protein)
    # # endHill = time.time()
    # # elapsedHill = endHill - startHill
    # # print("elapsed time Hill = " + elapsedHill)
    # # bestHillClimber.printBestHill()

if __name__ == "__main__":
    main()