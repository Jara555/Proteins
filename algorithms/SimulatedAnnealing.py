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
        self.maxDegrade = 1
        self.maxOverlap = 20
        self.maxTries = 10

        # set starting temperature
        self.temperature = 10

    def setParameter(self):
        """ Returns the starting parameter for variable T temperature,
        reflecting the starting temperature of the simulated annealing algorithm
        :return T: starting temperature """
        return self.temperature

    def skipOverlap(self, length=None):
        """ Overrides the standard method of the algorithm super class.
        Instead: Allows overlap for maxOverlap times in a row, in order to
         avoid local minima/maxima """

        # if overlap track and allow for maxOverlap times
        if self.protein.checkOverlap():
            self.overlapCount += 1
            self.trackOverlap += 1

            # if more than max times in a row: reset back to previous non-overlap pattern
            if self.trackOverlap > self.maxOverlap:
                # reset tracker
                self.trackOverlap = 0

                # reset protein
                self.foldPattern = copy.copy(self.nonOverlapPattern)
                self.resetProtein()

        # if no overlap and no degrade in stability: save non-overlap pattern
        else:
            if self.trackDegrade == 0:
                self.nonOverlapPattern = copy.copy(self.foldPattern)

    def handleDegradation(self):
        """ Overrides the standard method of the Hill climber class.
        Instead: Does not set pattern back to temp pattern when a degradation in
        stability is observed, but allows a degradation in stability for
        maxDegrade times in a row """

        # track stability degradation
        self.trackDegrade += 1

        # if max amount of degrading is reached reset to best pattern
        if self.trackDegrade > self.maxDegrade:
            # reset tracker
            self.trackDegrade = 0

            # reset protein
            self.foldPattern = copy.copy(self.bestPattern)
            self.resetProtein()

    def resetProtein(self):
        """ Resets best protein and stability values to folding pattern """

        self.protein.fold(self.foldPattern)
        self.protein.stability()
        self.bestStability = self.protein.stabilityScore
        self.bestPattern = copy.copy(self.foldPattern)

    def setEndState(self, T):
        """ Overrides the standard method of the Hill Climber class.
        Before folding a protein to best found pattern, makes sure there is no overlap
        and there is a better stability than the start stability.
        :param k: amount of times the run method is repeated """

        # check if endstate contains no overlap or worse stability score
        self.checkEndState(T)

        # make sure correct end values are saved
        self.protein.fold(self.bestPattern)
        self.protein.stability()
        self.bestStability = self.protein.stabilityScore
        self.bestRun = self.iterations

    def checkEndState(self, T):
        """ Checks end state of protein and tries to prevent overlap
        and degradations in end state
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
            if T > 0:
                print("\n>> OVERLAP IN END STATE: Get rid of overlap <<\n")

                # take pattern with best stability
                if self.nonOverlapScore <= self.startStability:
                    self.startPattern = copy.copy(self.nonOverlapPattern)

                # run again
                T = T - 1
                self.run(T)

            # if end state is reached with overlap, take last non overlap pattern
            else:
                print("\n>> COULD NOT GET RID OF OVERLAP: Finished with last known best non-overlap state <<\n")

                # if take pattern with best stability
                if self.nonOverlapScore <= self.startStability:
                    self.bestPattern = copy.copy(self.nonOverlapPattern)
                else:
                    self.bestPattern = copy.copy(self.startPattern)

        # if no overlap, but no improvement: try again
        elif self.startStability <= self.bestStability:
            if T > 0:
                print("\n>> COULD NOT FIND BETTER SOLUTION THAN START SOLUTION: Go again <<\n")

                # run again
                T = T - 1
                self.run(T)

            # end with start state
            else:
                print("\n>> COULD NOT FIND BETTER SOLUTION THAN START SOLUTION: Finished with start solution <<\n")
                self.bestPattern = copy.copy(self.startPattern)

        else:

            T = T - 1
            self.run(T)

