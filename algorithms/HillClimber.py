from random import randint
from classes.Algorithm import Algorithm
import copy


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
        self.startStability = 0

        self.foldPattern = copy.copy(self.startPattern)
        self.startValues()
        self.randStartStability = self.startStability

    def run(self, k=0):
        """ Runs algorithms and find pattern with higher stability, can be local minimum.
        :param k: only relevant for simulated annealing algorithm
                  keeps track of repetition number of run method
        :return: .csv file with the best folding patterns and associated stability's
        """
        # initialize start pattern
        self.foldPattern = copy.copy(self.startPattern)
        self.startValues()

        print("Started with a stability of " + str(self.bestStability))
        print()

        # calculate iteration range based on orientations to loop over
        iterationRange = int(self.maxIterations / (self.protein.dimensions * 2))

        # loop over max iterations
        for i in range(iterationRange):

            # pick random amino acid (starting from 2nd)
            randomAmino = self.pickAminoAcid()

            # loop over orientations
            for orientation in self.orientations:

                # keep track of iterations and print progress
                self.iterations += 1
                self.printProgress()

                # save copy of current folding pattern
                self.tempPattern = copy.copy(self.foldPattern)

                # change orientation of a random amino acid
                self.foldPattern[randomAmino] = orientation
                self.protein.fold(self.foldPattern)

                # if HC: skip if overlap detected
                if self.name == "HillClimber":
                    if self.skipOverlap():
                        continue
                # if SA: allow overlap for x times in a row
                elif self.name == "SimulatedAnnealing":
                    self.allowOverlap()

                # check for lower stability otherwise change back (HC) or allow x times (SA)
                self.stabilityCheck()

                # if ON write every pattern to csv
                self.writeCsvRow()

        # when SA finished: check end state for overlap or degradations
        if self.name == "SimulatedAnnealing":
            self.checkEndState(k)

        # make sure correct end values are saved
        self.protein.fold(self.bestPattern)
        self.protein.stability()
        self.bestStability = self.protein.stabilityScore
        self.bestRun = self.iterations

    def stabilityCheck(self):
        """ Checks if stability is lower or equal to best stability """

        # get stability scores
        self.protein.stability()

        # if stability is lower or equal make new best
        if self.protein.stabilityScore <= self.bestStability:
            self.bestStability = self.protein.stabilityScore
            self.bestPattern = copy.copy(self.foldPattern)
            self.bestRun = self.iterations
        else:
            # if HC: change back to previous
            if self.name == "HillClimber":
                self.foldPattern = copy.copy(self.tempPattern)
                self.protein.fold(self.foldPattern)
            # if SA: allow degradation in stability for x times
            elif self.name == "SimulatedAnnealing":
                self.allowDegrade()

    def startValues(self):
        """ Stores initial values. """

        # fold protein to start pattern and check stability
        self.protein.fold(self.foldPattern)
        self.protein.stability()

        # store start stablility
        self.startStability = self.protein.stabilityScore

        # stores start values as best values
        self.bestStability = self.protein.stabilityScore
        self.bestPattern = copy.copy(self.foldPattern)
        self.bestRun = self.iterations

    def pickAminoAcid(self):
        """ Picks a amino acid number at random ,
        or when 3 aminoacids in a row have the same orientation"""

        # pick random
        randomAmino = randint(2, self.protein.length - 1)

        # if a range with 3 in a row, pick that one!
        for j in range(1, self.protein.length - 1):
            if self.foldPattern[j - 1] == self.foldPattern[j] and \
                    self.foldPattern[j + 1] == self.foldPattern[j]:
                randomAmino = j

        return randomAmino