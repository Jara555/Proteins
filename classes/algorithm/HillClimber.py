from random import randint
from classes.Algorithms import Algorithms
import copy
import time
import csv


class HillClimber(Algorithms):
    """ Implements hill climber algorithms in order to most efficiently fold a protein """

    def __init__(self, protein, bestPattern, iterations, writeOptions):
        """ Declare all variables """
        Algorithms.__init__(self, protein)
        self.iterations = iterations
        self.writeOptions = writeOptions
        self.bestPattern = bestPattern
        self.maxStability = 0
        self.bestRun = 0
        self.overlapCount = 0
        self.run = 0
        self.elapsedTime = 0
        self.initialPattern = bestPattern
        self.foldPattern = []
        self.copyPattern = []
        self.initialStability = 0

    def runHillClimber(self):
        """ Runs algorithm and finds best pattern with highest stability """

        print()
        print("------------  Hill Climber started ----------------")
        print()

        start = time.time()

        # create csv file to write output to
        write_file = ('results/hillclimber' + str(self.protein.number) + '.csv')
        with open(write_file, 'w') as csvfile:
            fieldnames = ['run', 'stability', 'foldingPattern']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            # initialize folding pattern
            self.proteinFold(self.bestPattern)

            # initialize stability values
            self.initialValues()

            # visualize protein
            self.protein.visualize('protein')

            # initializes directions
            direction = ['+X', '-X', '+Y', '-Y']

            # runs hill climber i times
            for i in range(self.iterations):

                for j in range(len(direction)):

                    # prints progress per 1000 iterations
                    self.run += 1
                    if self.run % 1000 == 0:
                        print("Hill climber iteration: " + str(self.run) + "(stability " + str(self.maxStability) + ")")

                    # initialize random amino acid and direction
                    #randomDirection = randint(0, 3)
                    randomAmino = randint(0, self.protein.length - 1)

                    # saves copy of current folding pattern
                    self.copyPattern = copy.copy(self.foldPattern)

                    # changes direction of a random amino acid
                    self.foldPattern[randomAmino] = direction[j]
                    self.protein.fold(self.foldPattern)

                    # skip if overlap detected
                    if self.protein.checkOverlap(self.protein.length):
                        self.overlapCount += 1
                        continue

                    # check stability
                    self.protein.stability(self.protein.length)
                    self.stabilityCheck(i)

                    # if write all is on, write every solution to csv
                    if self.writeOptions == 0:
                        writer.writerow(
                            {'run': i, 'stability': self.protein.stabilityScore, 'foldingPattern': self.foldPattern})

        end = time.time()
        self.elapsedTime = end - start

        print("------------  Hill Climber finished ----------------")

    def printBestHill(self):

        # print info
        print()
        print(' HILL CLIMBER')
        print(' Maximal stability: ' + str(self.maxStability))
        print(' Total runs: ' + str(self.iterations))
        print(' First found in run: ' + str(self.bestRun))
        print(' Total overlap: ' + str(self.overlapCount))
        print(' Elapsed time: ' + "{0:.4f}".format(self.elapsedTime))
        print()

        # print fold pattern
        print('Initial fold pattern (stability ' + str(self.initialStability) + '): ')
        print(self.initialPattern)
        print()
        print('Best fold pattern: (stability ' + str(self.maxStability) + '): ')
        print(self.bestPattern)
        print()

        # plot protein
        self.protein.fold(self.bestPattern)
        self.protein.visualize(('Best random solution ' + str(self.maxStability)))

    def proteinFold(self, bestPattern):
        """ Folds protein. """

        self.foldPattern = []
        for k in range(self.protein.length):
            self.foldPattern.append(bestPattern[k])
        self.protein.fold(self.foldPattern)

    def stabilityCheck(self, i):
        """ Checks whether stability score is equal or better than max stability save new. """

        if self.protein.stabilityScore <= self.maxStability:
            self.maxStability = self.protein.stabilityScore
            self.bestPattern = copy.copy(self.foldPattern)
            self.bestRun = i
        else:
            for k in range(self.protein.length):
                self.foldPattern[k] = self.copyPattern[k]
            self.protein.fold(self.foldPattern)

    def initialValues(self):
        """ Stores initial values. """

        self.protein.stability(self.protein.length)

        # stores initial stablility
        self.initialStability = self.protein.stabilityScore

        # stores initial stability as maximum stability
        self.maxStability = self.protein.stabilityScore
