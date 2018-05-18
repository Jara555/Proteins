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

        self.nonOverlapPattern = startPattern

    def run(self, k=0):
        """ Runs algorithms and find pattern with higher stability, can be local minimum.

        :return: .csv file with the best folding patterns and associated stability's
        """

        # initialize patterns
        self.foldPattern = self.startPattern

        # initialize stability values
        self.initialValues()

        iterationRange = int(self.maxIterations / (self.protein.dimensions * 2))

        track = 0

        # loop over max iterations
        for i in range(iterationRange):

            # initialize random amino acid
            randomAmino = self.pickAminoAcid()

            # loop over orientations
            for orientation in self.orientations:

                # keep track of iterations and print progress
                self.iterations += 1
                self.printProgress()

                # save copy of current folding pattern
                self.copyPattern = copy.copy(self.foldPattern)

                # change orientation of a random amino acid
                self.foldPattern[randomAmino] = orientation
                self.protein.fold(self.foldPattern)

                # if self.skipOverlap():
                #     continue

                # allow overlap 1 time in a row
                if self.skipOverlap():
                    track += 1
                    # if more than 1 time in a row overlap go back to non overlap pattern
                    if track > 1:
                        track = 0
                        self.foldPattern = copy.copy(self.nonOverlapPattern)
                        continue
                # if no overlap detected, save overlap pattern
                else:
                    self.nonOverlapPattern = copy.copy(self.foldPattern)

                # check for lower of equal stability
                self.protein.stability(self.protein.length)
                self.stabilityCheck()

                # if ON write every pattern to csv
                self.writeCsvRow()

        self.checkEndState(k)

    def stabilityCheck(self):
        """ Returns true when stability score is lower or equal to max stability.

        !:param iteration: the iteration of this check
        :return boolean: true when stability score is lower
        """

        # if stability is lower or equal save new
        if self.protein.stabilityScore <= self.bestStability:
            self.bestStability = self.protein.stabilityScore
            self.bestPattern = copy.copy(self.foldPattern)
            self.bestRun = self.iterations
        else:
            for k in range(self.protein.length):
                self.foldPattern[k] = self.copyPattern[k]
            self.protein.fold(self.foldPattern)

    def initialValues(self):
        """ Stores initial values. """

        # fold protein and check stability
        self.protein.fold(self.startPattern)
        self.protein.stability(self.protein.length)

        # stores initial stablility
        self.initialStability = self.protein.stabilityScore

        # stores initial stability as maximum stability
        self.bestStability = self.protein.stabilityScore

    def pickAminoAcid(self):

        # pick random
        randomAmino = randint(0, self.protein.length - 1)

        # if a range with 4 in a row, pick that one!
        for j in range(2, self.protein.length - 2):
            if self.foldPattern[j - 1] == self.foldPattern[j] and \
                    self.foldPattern[j + 1] == self.foldPattern[j] and \
                    (self.foldPattern[j + 2] == self.foldPattern[j] or self.foldPattern[j - 2]):
                randomAmino = j

        return randomAmino

    def checkEndState(self, k):

        # prevent end state with overlap
        self.protein.fold(self.bestPattern)
        if self.protein.checkOverlap(self.protein.length):
            # if overlap in end state: repeat run for maximal 10 times
            if k < 10:
                k = k + 1

                # calculate stability of non overlap pattern
                self.protein.fold(self.nonOverlapPattern)
                self.protein.stability(self.protein.length)
                self.nonOverlapScore = self.protein.stabilityScore

                # if non overlap had lower or equal score than best stability
                # set to start pattern else keep original start pattern
                if self.nonOverlapScore <= self.initialStability:
                    self.startPattern = copy.copy(self.nonOverlapPattern)

                # run again
                self.run(k)

            # if end state is reached with overlap, take last non overlap pattern
            else:
                # calculate stability of non overlap pattern
                self.protein.fold(self.nonOverlapPattern)
                self.protein.stability(self.protein.length)
                self.nonOverlapScore = self.protein.stabilityScore

                # if non overlap had lower or equal score than best stability
                # set to start pattern else keep original start pattern
                if self.nonOverlapScore <= self.initialStability:
                    self.bestPattern = copy.copy(self.nonOverlapPattern)
                else:
                    self.bestPattern = self.startPattern

                self.protein.fold(self.bestPattern)
                self.bestStability = self.protein.stabilityScore
                self.bestRun = self.iterations
