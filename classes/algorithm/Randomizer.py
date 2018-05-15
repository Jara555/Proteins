import copy
import csv
import time
from random import random

from classes.Algorithms import Algorithms


class Randomizer(Algorithms):
    """ Randomize algorithm: finding best protein solution based on random patterns """

    def __init__(self, protein, iterations, writeCsv):
        """ Set and initiate all properties.

        :param protein: protein being folded
        :param iterations: how many random folding patterns should be generated
        :param writeOptions: 0 for write all solutions to .CSV-file, 1 for write only best solutions to .CSV-file
        :param dimensions: 2 for 2D or 3 for 3D
        """

        # input properties
        Algorithms.__init__(self, protein)
        self.iterations = iterations
        self.writeCsv = writeCsv

        # declare variables
        self.bestPattern = []
        self.bestStability = 0
        self.overlapCount = 0
        self.elapsed = 0
        self.bestRun = 0
        self.writer = None

        # get possible orientations based on dimension
        if self.protein.dimensions == 2:
            self.orientations = ['+X', '-X', '+Y', '-Y']
        elif self.protein.dimensions == 3:
            self.orientations = ['+X', '-X', '+Y', '-Y', '+Z', '-Z']

        #TODO: Make property of protein
        # set coordinates of first 2 aminoacids
        self.protein.list[0].setCoordinates(0, 0, 0)
        self.protein.list[1].setCoordinates(0, 1, 0)

        # starting fold pattern
        self.foldPattern = ['0'] + (['+Y'] * (self.protein.length - 1))

    def runRandomizer(self):
        """ Runs the randomizer to find a pattern with highest stability """

        print()
        print("------------  Randomizer started ----------------")
        print()

        # start timer
        start = time.time()

        # write to csv file
        if self.writeCsv == "ON":
            write_file = ('results/random' + str(self.protein.number) + '-' + str(self.protein.dimensions) + 'D.csv')
            with open(write_file, 'w') as csvfile:
                fieldnames = ['run', 'stability', 'foldingPattern']
                self.writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                self.writer.writeheader()

                # start iteration
                self.iterateRandomizer()

        # do not write to csv file
        elif self.writeCsv == "OFF":
            # start iteration
            self.iterateRandomizer()

        # end timer
        end = time.time()
        self.elapsed = end - start

        print()
        print("------------  Randomizer finished ----------------")
        print()

    def iterateRandomizer(self):
        """ Runs the iteration of the randomizer """

        for i in range(0, self.iterations):
            # generate random fold in protein
            self.generator()

            # check overlap
            if self.protein.checkOverlap(self.protein.length):
                self.overlapCount += 1
                continue

            # get stability score
            self.protein.stability(self.protein.length)

            # if ON write every pattern to csv
            if self.writeCsv == "ON":
                self.writer.writerow(
                    {'run': i, 'stability': self.protein.stabilityScore, 'foldingPattern': self.foldPattern})

            # if stability score is better write to best
            if self.protein.stabilityScore < self.bestStability:
                self.bestStability = self.protein.stabilityScore
                self.bestPattern = copy.copy(self.foldPattern)
                self.bestRun = i

            if i % (self.iterations * 0.05) == 0:
                print('Random iteration: ' + str(i) + '     (stability ' + str(
                    self.bestStability) + ')' + ' (foldpattern ' + str(self.bestPattern) + ')')

            # next iteration
            i += 1

    def generator(self):
        """ Creates random fold in protein based on pattern """

        # fold starts at aminoacid at index 2 (first are always ['0', '+Y' ... ] )
        i = 2

        # starting coordinates (of 2nd aminoacid)
        x, y, z = 0, 1, 0

        # iterate over folding pattern
        while i <= self.protein.length:
            # pick random orientation
            orientation = random.choice(self.orientations)

            # set coordinates to orientation + quick overlap check with previous aminoacid
            if orientation == '+X' and self.foldPattern[i - 1] != '-X':
                x += 1
            elif orientation == '-X' and self.foldPattern[i - 1] != '+X':
                x -= 1
            elif orientation == '+Y' and self.foldPattern[i - 1] != '-Y':
                y += 1
            elif orientation == '-Y' and self.foldPattern[i - 1] != '+Y':
                y -= 1
            elif orientation == '+Z' and self.foldPattern[i - 1] != '-Z':
                z += 1
            elif orientation == '-Z' and self.foldPattern[i - 1] != '+Z':
                z -= 1
            else:
                continue

            # save orientation in fold pattern and in coordinates of aminoacid
            self.foldPattern[i] = orientation
            self.protein.list[i].setCoordinates(x, y, z)
            i += 1

    def printBest(self):
        """ Prints the best found solution """

        # print info
        print()
        print('RANDOMIZER')
        print(' Maximal stability: ' + str(self.bestStability))
        print(' Total runs: ' + str(self.iterations))
        print(' First found in run: ' + str(self.bestRun))
        print(' Total overlap: ' + str(self.overlapCount))
        print(' Elapsed time: ' + "{0:.4f}".format(self.elapsed))
        print()

        # print fold pattern
        print('Fold pattern: ')
        print(self.bestPattern)
        print()

        # plot protein
        self.protein.fold(self.bestPattern)
        self.protein.visualize(('Best random solution ' + str(self.bestStability)))



