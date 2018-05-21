import copy

from algorithms.HillClimber import HillClimber


class SimulatedAnnealing(HillClimber):
    """ Subclass of HillClimber algorithm:
    Implements SimulatedAnnealing algorithm in order to efficiently fold a protein
    This algorithm allows overlap and a degradation of stability for several times
    in order to escape local maxima. The algorithm tries to end with a state without
    overlap and with a better solution than the starting state (if possible)"""

    def __init__(self, protein, writeCsv, maxIterations, startPattern):
        """ Set and initiate all properties.

        :param protein: protein to be fold
        :param writeCsv: writes solutions to .csv file when ON
        :param maxIterations: stop after maxIterations
        :param startPattern: start pattern of algorithm
        """

        HillClimber.__init__(self, protein, writeCsv, maxIterations, startPattern)
        self.name = "SimulatedAnnealing"

        # set overlap and stability degrading trackers
        self.trackOverlap = 0
        self.trackDegrade = 0

        # set non overlap pattern to start pattern
        self.nonOverlapPattern = copy.copy(self.startPattern)

        # TODO: Fine tune these parameters
        self.maxIterations = 1000
        self.maxDegrade = 1
        self.maxOverlap = 10
        self.maxTries = 10

    def allowOverlap(self):
        """ Allows overlap for a maxOverlap times in a row """

        # checks and tracks if overlap
        if self.skipOverlap():
            self.trackOverlap += 1
            # if more than .. times in a row overlap reset to non overlap pattern
            if self.trackOverlap > self.maxOverlap:
                self.trackOverlap = 0
                self.foldPattern = copy.copy(self.nonOverlapPattern)
                self.resetProtein()
        # if no overlap detected, save overlap pattern
        else:
            # only if no degrade in stability
            if self.trackDegrade == 0:
                self.nonOverlapPattern = copy.copy(self.foldPattern)

    def allowDegrade(self):
        """ Allows degrade in stability for a maxDegrade times in a row """

        # track stability degradation
        self.trackDegrade += 1

        # if max amount of degrading is reached reset to best pattern
        if self.trackDegrade > self.maxDegrade:
            self.trackDegrade = 0
            self.foldPattern = copy.copy(self.bestPattern)
            self.resetProtein()

    def resetProtein(self):
        """ Resets protein and stability values to current folding pattern """

        self.protein.fold(self.foldPattern)
        self.protein.stability()
        self.bestStability = self.protein.stabilityScore
        self.bestPattern = copy.copy(self.foldPattern)

    def checkEndState(self, k):
        """ Checks end state of protein and tries to prevent overlap in end state
         :param k: amount of times the run method is repeated """

        # fold protein to best found pattern
        self.protein.fold(self.bestPattern)

        # if overlap in best pattern try to eliminate this
        if self.protein.checkOverlap():

            # calculate stability of non overlap pattern
            self.protein.fold(self.nonOverlapPattern)
            self.protein.stability()
            self.nonOverlapScore = self.protein.stabilityScore

            # repeat run() several times in order to eliminate overlap
            if k < self.maxTries:
                print("\n>> OVERLAP IN END STATE: Get rid of overlap <<\n")

                # if non overlap had lower or equal score than best stability
                # set to start pattern else keep original start pattern
                if self.nonOverlapScore <= self.startStability:
                    self.startPattern = copy.copy(self.nonOverlapPattern)

                # run again
                k = k + 1
                self.run(k)

            # if end state is reached with overlap, take last non overlap pattern
            else:
                print("\n>> COULD NOT GET RID OF OVERLAP: Finished with last known best non-overlap state <<\n")

                # if non overlap had lower or equal score than best stability
                # set to start pattern else keep original start pattern
                if self.nonOverlapScore <= self.startStability:
                    self.bestPattern = copy.copy(self.nonOverlapPattern)
                else:
                    self.bestPattern = copy.copy(self.startPattern)

        # if no overlap, but no improvement: try again
        elif self.startStability <= self.bestStability:
            # try again
            if k < self.maxTries:
                print("\n>> COULD NOT FIND BETTER SOLUTION THAN START SOLUTION: Go again <<\n")

                # run again
                k = k + 1
                self.run(k)

            # end with start state
            else:
                print("\n>> COULD NOT FIND BETTER SOLUTION THAN START SOLUTION: Finished with start solution <<\n")
                self.bestPattern = copy.copy(self.startPattern)

