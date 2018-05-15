import copy
import random

from classes.Algorithm import Algorithm


class Randomizer(Algorithm):
    """ Implements Randomizer algorithm in order to efficiently fold a protein """

    def __init__(self, protein, writeCsv, iterations):
        """
        Set and initiate all properties.

        :param protein: protein to be fold
        :param writeCsv: writes solutions to .csv file when ON
        :param iterations: number of random folding patterns to be generated

        """

        # initialize input variables
        self.iterations = iterations

        # set class properties
        Algorithm.__init__(self, protein, writeCsv)
        self.name = "Randomizer"

    def runRandomizer(self):
        """ Runs the iterative core of the randomizer """

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
        while i < self.protein.length:
            # pick random orientation
            orientation = self.orientations[random.randrange(len(self.orientations))]

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




