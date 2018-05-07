from classes.Protein import Protein
from classes.algorithm.HillClimber import HillClimber
from classes.algorithm.BranchNBound import BranchNBound
from classes.algorithm.DepthFirst import DepthFirst
from classes.algorithm.Randomizer import Randomizer
import time

def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    # create protein of class Protein
    protein = Protein(2)

    protein.findHbonds()

    branchNbound = BranchNBound(protein)
    branchNbound.runBranchNBound()
    branchNbound.printBest()

    # depthFirst = DepthFirst(protein)
    # depthFirst.runDepthFirst()
    # depthFirst.printBest()

    # randomizer = Randomizer(protein, 1, 1)
    # randomizer.runRandomizer()
    # randomizer.printBest()
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