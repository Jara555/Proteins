import random
from classes.Algorithm import Algorithm
import copy


class HillClimber(Algorithm):
    """ Subclass of Algorithm:
    Implements hill climber algorithms in order to efficiently fold a protein.
    The algorithm picks aminoacids of a protein and changes their orientations
    while checking for any improvements in stability score. """

    def __init__(self, protein, writeCsv, maxIterations, startPattern):
        """ Set and initiate all properties.

        :param protein: protein to be fold
        :param writeCsv: writes solutions to .csv file when ON
        :param maxIterations: stop after maxIterations
        :param startPattern: start pattern of algorithm
        """

        # initialize input variables
        self.maxIterations = maxIterations

        # set class properties
        Algorithm.__init__(self, protein, writeCsv)
        self.name = "HillClimber"

        # declare lists
        self.foldPattern = startPattern
        self.tempPattern = []

    def run(self, param=None):
        """ Runs algorithms and find pattern with higher stability, can be local minimum.
        :return: .csv file with the best folding patterns and associated stability's (if ON)
        """

        # initialize start values
        self.startValues()

        print()
        print("Started with a stability of " + str(self.bestStability))
        print()

        # calculate iteration range based on orientations to loop over
        iterationRange = self.getIterationRange()

        # loop over max iterations
        for i in range(iterationRange):

            # pick amino acid (starting from 2nd)
            amino = self.getAminoAcid()

            # shuffle orientations
            orientationsShuffled = self.orientations
            random.shuffle(orientationsShuffled)

            # loop over shuffled orientations
            for orientation in orientationsShuffled:

                # keep track of iterations and print progress
                self.printProgress()
                self.iterations += 1

                # cool down temperature (only relevant in SA)
                self.coolDown()

                # save copy of current folding pattern
                self.tempPattern = copy.copy(self.foldPattern)

                # change orientation of a random amino acid
                self.foldPattern[amino] = orientation
                self.protein.fold(self.foldPattern)

                # check overlap and deal with it (differs for HC and SA)
                if self.skipOverlap():
                    self.handleOverlap()
                    continue

                # deal with stability score
                self.stabilityCheck()

                # if ON write every pattern to csv
                self.writeCsvRow()

    def startValues(self):
        """ Sets initial values and as best """

        # fold protein to pattern and check stability
        self.protein.fold(self.foldPattern)
        self.protein.stability()

        # store values as best values
        self.bestPattern = copy.copy(self.foldPattern)
        self.bestStability = self.protein.stabilityScore
        self.startStability = self.protein.stabilityScore

    def getIterationRange(self):
        """ Calculate amount of iterations based on max iterations input
        and orientations to loop over. """
        return int(self.maxIterations / (self.protein.dimensions * 2))

    def getAminoAcid(self):
        """ Picks an amino acid between 2nd and last at random ,
        or when 3 aminoacids in a row have the same orientation"""

        # pick random
        amino = random.randint(2, self.protein.length - 1)

        # if a range with 3 in a row, pick that one
        for j in range(1, self.protein.length - 1):
            if self.foldPattern[j - 1] == self.foldPattern[j] and \
                    self.foldPattern[j + 1] == self.foldPattern[j]:
                amino = j

        return amino

    def stabilityCheck(self):
        """ Checks if stability is lower or equal to best stability, otherwise
         handles the degradation in stability"""

        # get stability scores
        self.protein.stability()

        # if stability is lower or equal make new best
        if self.protein.stabilityScore <= self.bestStability:
            self.bestStability = self.protein.stabilityScore
            self.bestPattern = copy.copy(self.foldPattern)
            if self.protein.stabilityScore < self.bestStability:
                self.bestRun = self.iterations
        else:
            # deal with degradation (differs for HC and SA)
            self.handleDegradation()

    def handleOverlap(self):
        """ Stores temp pattern back to fold pattern to undo overlap
        !Overridden in simulated annealing! """

        self.foldPattern = copy.copy(self.tempPattern)
        self.protein.fold(self.foldPattern)

    def handleDegradation(self):
        """ Stores temp pattern back to fold pattern to undo degradation in stability
        !Overridden in simulated annealing! """

        self.foldPattern = copy.copy(self.tempPattern)
        self.protein.fold(self.foldPattern)

    def coolDown(self):
        """ Method is only relevant in SA algorithm """
        pass
