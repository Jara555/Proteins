import copy
import csv
import time
from random import randint
from classes.Algorithms import Algorithms


class Randomizer(Algorithms):
    """ Randomize algorithm: finding best protein solution based on random patterns """

    def __init__(self, protein, iterations, writeOptions):
        Algorithms.__init__(self, protein)
        self.iterations = iterations
        self.writeOptions = writeOptions
        self.bestPattern = []
        self.maxStability = 0
        self.firstHit = 0
        self.overlapCount = 0

    def runRandomizer(self):
        """ Runs the randomizer and finds best pattern with highest stability """

        start = time.time()

        # create csv file to write output to
        write_file = ('results/random' + str(self.protein.number) + '.csv')
        with open(write_file, 'w') as csvfile:
            fieldnames = ['run', 'stability', 'foldingPattern']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            i = 0

            # iterate till the repetition amount is reached
            while i <= self.iterations:
                # get random folding pattern and fold protein according to this pattern
                self.generator()
                self.protein.fold(self.foldPattern)

                # skip if overlap detected
                if self.protein.checkOverlap(self.protein.length):
                    self.overlapCount += 1
                    continue

                # get stability score of input protein
                self.protein.stability()

                # if write all is on, write every solution to csv
                if self.writeOptions == 0:
                    writer.writerow(
                        {'run': i, 'stability': self.protein.stabilityScore, 'foldingPattern': self.foldPattern})

                # if stability score is equal or better than max stability save new
                if self.protein.stabilityScore <= self.maxStability:
                    self.maxStability = self.protein.stabilityScore
                    self.bestPattern = copy.copy(self.foldPattern)
                    self.bestProtein = copy.copy(self.protein)

                    # if write all is off, write only best solutions to csv
                    if self.writeOptions == 1:
                        writer.writerow(
                            {'run': i, 'stability': self.protein.stabilityScore, 'foldingPattern': self.foldPattern})

                # next iteration
                i += 1

            end = time.time()
            self.elapsed = end - start

    def runFastRandomizer(self):
        """ Runs the randomizer and finds best pattern with highest stability
            Without writing to a csv! """

        start = time.time()

        i = 0

        # iterate till the repetition amount is reached
        while i <= self.iterations:
            # get random folding pattern and fold protein according to this pattern
            self.fastGenerator()

            # skip if overlap detected
            if self.protein.checkOverlap(self.protein.length):
                self.overlapCount += 1
                continue

            # get stability score of input protein
            self.protein.stability()

            # if stability score is equal or better than max stability save new
            if self.protein.stabilityScore < self.maxStability:
                self.maxStability = self.protein.stabilityScore
                self.bestPattern = copy.copy(self.foldPattern)
                self.bestProtein = copy.copy(self.protein)
                self.bestRun = i

            # next iteration
            i += 1

        end = time.time()
        self.elapsed = end - start

    def generator(self):
        """ Creates random folding pattern """

        # get size protein and create empty folding pattern list
        size = self.protein.length
        self.foldPattern = []

        # index counter for number of items in folding pattern list
        i = 0

        # iterate over required folding pattern length
        while i <= size:

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

    def fastGenerator(self):
        """ Creates random folding pattern """

        # create folding pattern list
        self.foldPattern = [0] * self.protein.length
        self.foldPattern[1] = 3

        # index counter for number of items in folding pattern list
        x, y = 0, 1
        i = 2

        self.protein.list[0].setCoordinates(0, 0)
        self.protein.list[1].setCoordinates(x, y)

        # iterate over required folding pattern length
        while i < self.protein.length:

            orientation = randint(1, 4)

            if orientation == 1 and self.foldPattern[i - 1] != 2:
                self.foldPattern[i] = orientation
                x = x + 1
                self.protein.list[i].setCoordinates(x, y)
                i += 1
            elif orientation == 2 and self.foldPattern[i - 1] != 1:
                self.foldPattern[i] = orientation
                x = x - 1
                self.protein.list[i].setCoordinates(x, y)
                i += 1
            elif orientation == 3 and self.foldPattern[i - 1] != 4:
                self.foldPattern[i] = orientation
                y = y + 1
                self.protein.list[i].setCoordinates(x, y)
                i += 1
            elif orientation == 4 and self.foldPattern[i - 1] != 3:
                self.foldPattern[i] = orientation
                y = y - 1
                self.protein.list[i].setCoordinates(x, y)
                i += 1

    def printBest(self):
        """ Prints the best found solution """

        # print info
        print()
        print('RANDOMIZER')
        print(' Maximal stability: ' + str(self.maxStability))
        print(' First found in run: ' + str(self.bestRun))
        print(' Total overlap: ' + str(self.overlapCount))
        print(' Elapsed time: ' + "{0:.4f}".format(self.elapsed))
        print()

        # print fold pattern
        print('Fold pattern: ')
        print(self.bestPattern)
        print()

        # plot protein
        self.bestProtein.visualize(('Best random solution ' + str(self.maxStability)))



