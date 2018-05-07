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
        self.initialPattern = []
        self.foldPattern = []
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

            # folds protein according to best fold from randomizer
            self.initialPattern = self.bestPattern
            print('Initial pattern: ' + str(self.initialPattern))
            self.generator()
            for k in range(self.protein.length):
                self.foldPattern[k] = self.bestPattern[k]
            self.protein.fold(self.foldPattern)

            print('Foldpattern: ' + str(self.foldPattern))

            # determines stability
            self.protein.stability()
            self.maxStability = self.protein.stabilityScore

            self.protein.stability()
            self.initialStability = self.protein.stabilityScore

            # initializes directions
            direction = ['+X', '-X', '+Y', '-Y']

            # runs hill climber i times
            for i in range(self.iterations):

                # prints progress per 1000 iterations
                self.run += 1
                if self.run % 1000 == 0:
                    print()
                    print("Hill climber iteration: " + str(self.run) + "(stability " + str(self.maxStability) + ")")
                    print()

                # initialize random amino acid and direction
                randomDirection = randint(0, 3)
                randomAmino = randint(1, 7)

                # saves copy of current folding pattern
                copyFold = self.foldPattern

                # changes direction of a random amino acid
                self.foldPattern[randomAmino] = direction[randomDirection]
                self.protein.fold(self.foldPattern)

                # skip if overlap detected
                if self.protein.checkOverlap(self.protein.length):
                    self.overlapCount += 1
                    continue

                # check stability
                self.protein.stability()

                # if write all is on, write every solution to csv
                if self.writeOptions == 0:
                    writer.writerow(
                        {'run': i, 'stability': self.protein.stabilityScore, 'foldingPattern': self.foldPattern})

                # if stability score is equal or better than max stability save new
                if self.protein.stabilityScore <= self.maxStability:
                    self.maxStability = self.protein.stabilityScore
                    self.bestPattern = copy.copy(self.foldPattern)
                    self.bestRun = i
                else:
                    for k in range(self.protein.length):
                        self.foldPattern[k] = copyFold[k]
                    self.protein.fold(self.foldPattern)

        end = time.time()
        self.elapsedTime = end - start

        print()
        print()
        print("------------  Hill Climber finished ----------------")
        print()

    def printBestHill(self):

        # print info
        print()
        print('HILL CLIMBER')
        print('Maximal stability: ' + str(self.maxStability))
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

    def generator(self):
        """ Creates random folding pattern """

        # get size protein and create empty folding pattern list
        size = self.protein.length
        self.foldPattern = []

        # index counter for number of items in folding pattern list
        i = 0

        # iterate over required folding pattern length
        while i < size:

            # get random orientation index
            orientation = randint(1, 4)

            # first element should always start at 0, 0
            if i == 0:
                self.foldPattern.append('0')
                self.foldPattern.append('+Y')
                i += 2
            # previous element cannot be inverse orientation
            elif orientation == 1 and self.foldPattern[i - 1] != '-X':
                self.foldPattern.append('+X')
                i += 1
            elif orientation == 2 and self.foldPattern[i - 1] != '+X':
                self.foldPattern.append('-X')
                i += 1
            elif orientation == 3 and self.foldPattern[i - 1] != '-Y':
                self.foldPattern.append('+Y')
                i += 1
            elif orientation == 4 and self.foldPattern[i - 1] != '+Y':
                self.foldPattern.append('-Y')
                i += 1
