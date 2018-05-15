from random import randint
from classes.Algorithms import Algorithms
import copy
import time
import csv


class HillClimber(Algorithms):
    """ Implements hill climber algorithms in order to efficiently fold a protein."""

    def __init__(self, protein, bestPattern, iterations, writeOptions, dimensions):
        """ Set and initiate all properties.

        :param protein: protein being folded
        :param bestPattern: the pattern to start with
        :param iterations: how many times should the hillclimber run
        :param writeOptions: 0 for write all solutions to .CSV-file, 1 for write only best solutions to .CSV-file
        :param dimensions: 2 for 2D or 3 for 3D
        """
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
        #self.foldPattern = bestPattern
        self.copyPattern = []
        self.initialStability = 0
        self.dimensions = dimensions

    def runHillClimber(self):
        """ Runs algorithm and find pattern with higher stability, can be local minimum.

        :return: .csv file with the best folding patterns and associated stability's
        """

        print()
        print("------------  Hill Climber started ----------------")
        print()

        start = time.time()

        # create csv file to write output to
        write_file = ('results/hillclimber' + str(self.protein.number) + '.' + str(self.dimensions) + 'D.csv')
        with open(write_file, 'w') as csvfile:
            fieldnames = ['run', 'stability', 'foldingPattern']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            # initialize folding pattern
            self.proteinFold(self.bestPattern)

            # initialize stability values
            self.initialValues()

            # initializes directions
            if self.dimensions == 2:
                direction = ['+X', '-X', '+Y', '-Y']
            if self.dimensions == 3:
                direction = ['+X', '-X', '+Y', '-Y', '+Z', '-Z']

            # runs hill climber i times
            for i in range(self.iterations):

                for j in range(len(direction)):

                    # prints progress per 1000 iterations
                    self.run += 1
                    if self.run % 50000 == 0:
                        print("Hill climber iteration: " + str(self.run) + "(stability " + str(self.maxStability) + ")")

                    # initialize random amino acid and direction
                    #randomDirection = randint(0, 3)
                    randomAmino = randint(0, self.protein.length - 1)

                    # saves copy of current folding pattern
                    self.copyPattern = copy.copy(self.foldPattern)

                    # changes direction of a random amino acid
                    self.foldPattern[randomAmino] = direction[j]
                    self.protein.fold(self.foldPattern, self.dimensions)

                    # skip if overlap detected
                    if self.protein.checkOverlap(self.protein.length):
                        self.overlapCount += 1
                        continue

                    # check stability
                    self.protein.stability(self.protein.length, self.dimensions)
                    self.stabilityCheck(i)

                    # if write all is on, write every solution to csv
                    if self.writeOptions == 0:
                        writer.writerow(
                            {'run': i, 'stability': self.protein.stabilityScore, 'foldingPattern': self.foldPattern})

        end = time.time()
        self.elapsedTime = end - start

        print("------------  Hill Climber finished ----------------")

    def printBestHill(self):
        """ Print and visualise the initial and best foldingPattern found.

        :return: print foldingPatterns and associated stability's to terminal, visualise in plot
        """

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

        # plot proteins
        self.protein.fold(self.initialPattern, self.dimensions)
        self.protein.visualize(
            str(self.protein.number) + '. Initial folding (stability ' + str(self.initialStability) + ')',
            self.dimensions)
        self.protein.fold(self.bestPattern, self.dimensions)
        self.protein.visualize((str(self.protein.number) + '. Best HillClimber solution ' + str(self.maxStability)),
            self.dimensions)

    def proteinFold(self, bestPattern):
        """ Folds protein.

        :param bestPattern: the best folding pattern so far
        """

        for k in range(self.protein.length):
            self.foldPattern.append(bestPattern[k])
        self.protein.fold(self.foldPattern, self.dimensions)

    def stabilityCheck(self, iteration):
        """ Checks whether stability score is equal or better than max stability save new.

        !:param iteration: the iteration of this check
        """

        if self.protein.stabilityScore <= self.maxStability:
            self.maxStability = self.protein.stabilityScore
            self.bestPattern = copy.copy(self.foldPattern)
            self.bestRun = iteration
        else:
            for k in range(self.protein.length):
                self.foldPattern[k] = self.copyPattern[k]
            self.protein.fold(self.foldPattern, self.dimensions)

    def initialValues(self):
        """ Stores initial values. """

        self.protein.fold(self.foldPattern, self.dimensions)
        self.protein.stability(self.protein.length, self.dimensions)

        # stores initial stablility
        self.initialStability = self.protein.stabilityScore

        # stores initial stability as maximum stability
        self.maxStability = self.protein.stabilityScore
