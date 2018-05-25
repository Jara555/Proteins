import copy
import random

from algorithms.HillClimber import HillClimber


class SimulatedAnnealing(HillClimber):
    """ Subclass of HillClimber algorithm:
    Implements SimulatedAnnealing algorithm in order to efficiently fold a protein.
    This algorithm tries to escape local maxima/minima by calculating a probability score
    for the acceptance of degradations. Degradations can be in the form of overlap or
    stability. Together with a temperature which is cooling down every iteration,
    the probability scores for both degradations can be calculated."""

    def __init__(self, protein, writeCsv, maxIterations, startPattern):
        """ Set and initiate all properties.

        :param protein: protein to be fold
        :param writeCsv: writes solutions to .csv file when ON
        :param maxIterations: stop after maxIterations
        :param startPattern: start pattern of algorithm
        """

        HillClimber.__init__(self, protein, writeCsv, maxIterations, startPattern)
        self.name = "SimulatedAnnealing"

        # set overlap tracker and maximal allowed overlap
        self.trackOverlap = 0
        self.maxOverlap = 50

        # set starting temperature
        self.temperature = 100

    def coolDown(self):
        """ Calculates the temperature for every iteration. The cooling down of
        the temperature follows a linear pattern."""

        self.temperature = self.temperature - (self.temperature - self.maxIterations)

    def handleOverlap(self):
        """ Overrides the standard method of the Hill climber class.
        Instead: Does not set pattern back to temp pattern when overlap
        is observed, accepts overlap based on a probability score """

        # track overlap
        self.trackOverlap += 1

        # calculate probability score and random probability
        probabilityScore = self.getOverlapProbability()
        randomProbability = random.random()

        # only accepts overlap if probability is higher than random
        if probabilityScore < randomProbability:
            # reset tracker
            self.trackOverlap = 0

            # go back to best pattern (without overlap)
            self.foldPattern = copy.copy(self.bestPattern)
            self.protein.fold(self.foldPattern)

    def handleDegradation(self):
        """ Overrides the standard method of the Hill climber class.
        Instead: Does not set pattern back to temp pattern when a degradation in
        stability is observed, but allows a degradation in stability for
        maxDegrade times in a row """

        # calculates probability score and random probability
        probabilityScore = self.getStabilityProbability()
        randomProbability = random.random()

        # only accepts overlap if probability is higher than random
        if probabilityScore < randomProbability:
            # reset tracker
            self.trackOverlap = 0

            # go back to best pattern (without overlap)
            self.foldPattern = copy.copy(self.bestPattern)
            self.protein.fold(self.foldPattern)

    def getOverlapProbability(self):
        """ Calculates the probability of acceptance, based on an overlap score and
        the current temperature. The overlap score depends on the amount of overlap
        tracked in a row and the maximal allowed overlap.
        :return: probability of acceptance """

        # score of overlap tracker and max allowed overlap
        score = (self.maxOverlap - self.trackOverlap) / self.maxOverlap

        # probability is based on overlap score and temp
        return (score * self.temperature) / 100

    def getStabilityProbability(self):
        """ Calculates the probability of acceptance, based on an stability score and
        the current temperature. The stability score depends on the difference between
        the new stability and the best found stability.
        :return: probability of acceptance """

        # score of difference in best stability and current stability score
        score = 1 - ((self.bestStability - self.protein.stabilityScore)
                     / self.bestStability)

        # probability is based on stability score and temp
        return (score * self.temperature) / 100
