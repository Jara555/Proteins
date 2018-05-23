import copy
import math

from algorithms.HillClimber import HillClimber


class SimulatedAnnealing(HillClimber):
    """ Subclass of HillClimber algorithm:
    Implements SimulatedAnnealing algorithm in order to efficiently fold a protein
    This algorithm allows overlap and a degradation of stability for a maximal amount of
    times in order to escape local minima/maxima. Temperature is cooling down every run,
    as is the allowed amount of overlap and degradation."""

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

        # set starting temperature
        self.temperature = 100

        # set start of max degradation and overlap scores based on temperature
        self.maxOverlap = self.temperature / 2
        self.maxDegrade = self.temperature / 10

    def getParameter(self):
        """ Returns the starting parameter for variable T temperature,
        reflecting the starting temperature of the simulated annealing algorithm
        :return T: starting temperature """
        return self.temperature

    def getIterationRange(self):
        """ Calculate amount of iterations based on max iterations input and
        orientations to loop over. Divide amount by 10, because at 1/10th
        the temperature of the simulated annealing has to drop."""
        return int((self.maxIterations / (self.protein.dimensions * 2)) / 10)

    def handleOverlap(self):
        """ Overrides the standard method of the Hill climber class.
        Instead: Does not set pattern back to temp pattern when overlap
        is observed, but allows this for maxOverlap times in a row """

        self.trackOverlap += 1

        # if max amount of overlap is reached reset to best pattern
        if self.trackOverlap >= self.maxOverlap:
            # reset tracker and protein
            self.trackOverlap = 0
            self.foldPattern = copy.copy(self.bestPattern)

    def handleDegradation(self):
        """ Overrides the standard method of the Hill climber class.
        Instead: Does not set pattern back to temp pattern when a degradation in
        stability is observed, but allows a degradation in stability for
        maxDegrade times in a row """

        self.trackDegrade += 1

        # if max amount of degrading is reached reset to best pattern
        if self.trackDegrade >= self.maxDegrade:
            # reset tracker and protein
            self.trackDegrade = 0
            self.foldPattern = copy.copy(self.bestPattern)

    def setEndState(self, T):
        """ Overrides the standard method of the Hill Climber class.
        Adjust temperature and calculates new max overlap and degradation levels.
        :param T: Temperature of the algorithm (cooling down every iteration) """

        # cooling down temperature according to temperature function and go again
        T = self.getTemperature(T)

        # every time T is cooling down overlap and degradation borders move with it
        self.maxDegrade = self.getMaxDegrade(T)
        self.maxOverlap = self.getMaxOverlap(T)

        # as long as temperature is not lower than 0 keep on going
        if T > 0:
            # start with best pattern so far
            self.startPattern = copy.copy(self.bestPattern)

            # reset trackers
            self.trackOverlap = 0
            self.trackDegrade = 0

            # go again
            self.run(T)

        # make sure correct end values are saved
        self.protein.fold(self.bestPattern)
        self.protein.stability()
        self.bestStability = self.protein.stabilityScore

    def getTemperature(self, T):
        """ Calculates the cooling down in temperature
        :param T: Temperature of the algorithm (cooling down every iteration)
        :return T: New temperature"""

        return T - (self.temperature / 10)

    def getMaxOverlap(self, T):
        """ Calculates the exponential decay of maxOverlap based on temperature
        :param T: Temperature of the algorithm (cooling down every iteration)
        :return maxOverlap: maximal allowed overlap level"""

        if T > 0:
            return int(math.log(T) * 10)
        else:
            return int(self.maxOverlap / 2)

    def getMaxDegrade(self, T):
        """ Calculates the decay of maxDegrade based on temperature
        :param T: Temperature of the algorithm (cooling down every iteration)
        :return maxDegrade: maximal allowed degradation level"""

        return T / (self.temperature / 10)



