from random import randint
from classes.Algorithm import Algorithm
import copy


class HillClimber(Algorithm):
    """ Subclass of Algorithm:
    Implements hill climber algorithms in order to efficiently fold a protein."""

    def __init__(self, protein, writeCsv, maxIterations, startPattern):
        """ Set and initiate all properties.
        :param protein: protein being folded
        :param bestPattern: the pattern to start with
        :param iterations: how many times should the hillclimber run
        :param writeOptions: 0 for write all solutions to .CSV-file, 1 for write only best solutions to .CSV-file
        :param dimensions: 2 for 2D or 3 for 3D
        """

        # initialize input variables
        self.maxIterations = maxIterations
        self.startPattern = startPattern

        # set class properties
        Algorithm.__init__(self, protein, writeCsv)
        self.name = "HillClimber"

        # declare other variables
        self.foldPattern = []
        self.copyPattern = []
        self.initialStability = 0

    def run(self, param=None):
        """ Runs algorithms and find pattern with higher stability, can be local minimum.
        :return: .csv file with the best folding patterns and associated stability's
        """

        # initialize patterns
        self.bestPattern = self.startPattern
        self.proteinFold()

        # initialize stability values
        self.initialValues()

        iterationRange = int(self.maxIterations / (self.protein.dimensions * 2))

        # loop over max iterations
        for i in range(iterationRange):

            # initialize random amino acid
            randomAmino = randint(0, self.protein.length - 1)

            # loop over orientations
            for j in range(len(self.orientations)):

                # keep track of iterations and print progress
                self.iterations += 1
                self.printProgress()

                # save copy of current folding pattern
                self.copyPattern = copy.copy(self.foldPattern)

                # change orientation of a random amino acid
                self.foldPattern[randomAmino] = self.orientations[j]
                self.protein.fold(self.foldPattern)

                # skip if overlap detected
                if self.skipOverlap():
                    continue

                # check for lower stability otherwise change back
                self.stabilityCheck()

                # if ON write every pattern to csv
                self.writeCsvRow()

    def proteinFold(self):
        """ Folds protein.
        :param bestPattern: the best folding pattern so far
        """

        # loop over protein length
        for i in range(self.protein.length):
            # append best pattern to fold pattern
            self.foldPattern.append(self.bestPattern[i])

        # fold protein
        self.protein.fold(self.foldPattern)

    def stabilityCheck(self):
        """ Returns true when stability score is lower or equal to max stability.
        !:param iteration: the iteration of this check
        :return boolean: true when stability score is lower
        """

        # if stability is lower return true
        if self.protein.stabilityScore <= self.bestStability:
            self.checkBest()
        else:
            # change back
            for k in range(self.protein.length):
                self.foldPattern[k] = self.copyPattern[k]
            self.protein.fold(self.foldPattern)

    def initialValues(self):
        """ Stores initial values. """

        # fold protein and check stability
        self.protein.fold(self.foldPattern)
        self.protein.stability(self.protein.length)

        # stores initial stablility
        self.initialStability = self.protein.stabilityScore

        # stores initial stability as maximum stability
        self.bestStability = self.protein.stabilityScore