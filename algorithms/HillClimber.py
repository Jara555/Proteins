from random import randint
from classes.Algorithm import Algorithm
import copy
import random


class HillClimber(Algorithm):
    """ Subclass of Algorithm:
    Implements hill climber algorithms in order to efficiently fold a protein."""

    def __init__(self, protein, writeCsv, maxIterations, startPattern):
        """ Set and initiate all properties.

        :param protein: protein to be fold
        :param writeCsv: writes solutions to .csv file when ON
        :param maxIterations: stop after maxIterations
        :param startPattern: start pattern of algorithm
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

        # save starting stability
        self.protein.fold(self.startPattern)
        self.protein.stability()
        self.startStability = self.protein.stabilityScore

    def run(self, T):
        """ Runs algorithms and find pattern with higher stability, can be local minimum.
        :param T: only relevant for simulated annealing algorithm
                Temperature of the algorithm (cooling down every iteration)
        :return: .csv file with the best folding patterns and associated stability's
        """

        # initialize fold pattern and values
        self.foldPattern = copy.copy(self.startPattern)
        self.startValues()

        print()
        print("Started with a stability of " + str(self.bestStability))
        print()

        # calculate iteration range based on orientations to loop over
        iterationRange = self.getIterationRange()

        # loop over max iterations
        for i in range(iterationRange):

            # pick random amino acid (starting from 2nd)
            amino = self.getAminoAcid()

            # loop over orientations in random order
            orientationsShuffled = self.orientations
            random.shuffle(orientationsShuffled)
            for orientation in orientationsShuffled:

                # keep track of iterations and print progress
                self.printProgress()
                self.iterations += 1

                # save copy of current folding pattern
                self.tempPattern = copy.copy(self.foldPattern)

                # change orientation of a random amino acid
                self.foldPattern[amino] = orientation
                self.protein.fold(self.foldPattern)

                # standard overlap method (!overridden for SA!)
                if self.skipOverlap():
                    self.handleOverlap()
                    continue

                # deal with stability score
                self.stabilityCheck()

                # if ON write every pattern to csv
                self.writeCsvRow()

        # set (and check) end state
        self.setEndState(T)

    def getIterationRange(self):
        """ Calculate amount of iterations based on max iterations input and
        orientations to loop over. """
        return int(self.maxIterations / (self.protein.dimensions * 2))

    def stabilityCheck(self):
        """ Checks if stability is lower or equal to best stability """

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
        """ Stores temp pattern back to fold pattern
        to undo overlap
        Overridden in simulated annealing! """

        self.foldPattern = copy.copy(self.tempPattern)
        self.protein.fold(self.foldPattern)

    def handleDegradation(self):
        """ Stores temp pattern back to fold pattern
        to undo a degradation in stability score
        Overridden in simulated annealing! """

        self.foldPattern = copy.copy(self.tempPattern)
        self.protein.fold(self.foldPattern)

    def startValues(self):
        """ Stores initial values and sets to best """

        # fold protein to pattern and check stability
        self.protein.fold(self.foldPattern)
        self.protein.stability()

        # stores start values as best values
        self.bestPattern = copy.copy(self.foldPattern)
        self.bestStability = self.protein.stabilityScore

    def getAminoAcid(self):
        """ Picks an amino acid number at random ,
        or when 3 aminoacids in a row have the same orientation"""

        # pick random
        amino = randint(2, self.protein.length - 1)

        # if a range with 3 in a row, pick that one!
        for j in range(1, self.protein.length - 1):
            if self.foldPattern[j - 1] == self.foldPattern[j] and \
                    self.foldPattern[j + 1] == self.foldPattern[j]:
                amino = j

        return amino

    def setEndState(self, T):
        """ Folds protein to best found pattern in order to
         end in the best found state
         Overridden in simulated annealing! """

        # make sure correct end values are saved
        self.protein.fold(self.bestPattern)
        self.protein.stability()
        self.bestStability = self.protein.stabilityScore

