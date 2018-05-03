from classes.Algorithms import Algorithms
from algorithm.Randomizer import Randomizer
from random import randint


class HillClimber(Algorithms):
    """ Implements hill climber algorithms in order to most efficiently fold a protein """

    def __init__(self, protein):
        Algorithms.__init__(self, protein)
        self.foldingPattern = ['0', '+Y', '+Y', '-X', '-X', '+Y', '-X', '-Y', '-Y', '+X', '+X', '-Y', '-Y', '+X']
        self.lengthProtein = len(protein)

    def runHillClimber(self, protein):
        """ Runs algorithm and finds best pattern with highest stability """

        stabilityScore = self.protein.stability();

        for i in range(1000):
            randomDirection = randint(0, 3)
            randomAmino = randint(0, 7)
            direction = ['+X', '-X', '+y', '-Y']

            tempProtein = protein
            foldingTemp = tempProtein[randomDirection] = direction[randomAmino]

            if self.tempProtein.checkOverlap(self.protein.length):
                if self.protein.stability() < stabilityScore:
                    protein = tempProtein

        return protein


    def printBestHill(self):

        print()
        print('Hill maximal stability: ' + str(self.maxStability))
        print(self.bestPattern)
        print()

        # plot protein
        self.protein.fold(self.bestPattern)
        self.protein.visualize(('Best hill solution ' + str(self.maxStability)))



























