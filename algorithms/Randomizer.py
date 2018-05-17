import random

from classes.Algorithm import Algorithm


class Randomizer(Algorithm):
    """ Subclass of Algorithm:
    Implements Randomizer algorithms in order to efficiently fold a protein """

    def __init__(self, protein, writeCsv, maxIterations):
        """
        Set and initiate all properties.

        :param protein: protein to be fold
        :param writeCsv: writes solutions to .csv file when ON
        :param maxIterations: number of random folding patterns to be generated

        """

        # initialize input variables
        self.maxIterations = maxIterations

        # set class properties
        Algorithm.__init__(self, protein, writeCsv)
        self.name = "Randomizer"

    def run(self, param=None):
        """ Runs a maximal amount of iterations
        in which random folding patterns are created """

        # loop over max iterations
        for self.iterations in range(self.maxIterations):

            # generate random fold in protein
            self.generator()

            # print progress in terminal
            self.printProgress()

            # skip if overlap detected
            if self.skipOverlap():
                continue

            # check for lower stability
            self.checkBest()

            # if ON write every pattern to csv
            self.writeCsvRow()

    def generator(self):
        """ Creates random fold in protein and saves pattern """

        # fold starts at aminoacid 2 (always ['0', '+Y' ... ] )
        i = 2

        # starting coordinates of 2nd aminoacid)
        x, y, z = 0, 1, 0

        # iterate over folding pattern
        while i < self.protein.length:
            # pick random orientation out of list
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

            # iterate 1
            i += 1




