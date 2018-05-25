from algorithms.DepthFirst import DepthFirst
from algorithms.Randomizer import Randomizer


class BranchNBound(DepthFirst):
    """ Subclass of DepthFirst algorithms:
    Implements Branch 'N Bound algorithms in order to efficiently fold a protein """

    def __init__(self, protein, writeCsv, maxIterations):
        """ Set and initiate all properties.

        :param protein: protein to be fold
        :param writeCsv: writes solutions to .csv file when ON
        :param maxIterations: stop after maxIterations
        """

        # run randomizer to have a initialStability
        randomAlgorithm = Randomizer(protein, "OFF", 10000)
        randomAlgorithm.run()
        bestStability = randomAlgorithm.bestStability
        bestPattern = randomAlgorithm.bestPattern

        # set class properties
        DepthFirst.__init__(self, protein, writeCsv, maxIterations=None)
        self.name = "BranchNBound"

        #initialize input varialbes
        self.maxIterations = maxIterations

        # set initial stability
        self.bestStability = bestStability
        self.bestPattern = bestPattern

        self.protein.findBonds()

    def pruneStability(self, k):
        """ Prunes on basis of the expected stability of the protein
                Calculated in the prune method of the protein class

        :param k: the aminoacid currently being placed
        :return: True if branch is pruned, False if not pruned
        """

        # if stability target can not be reached anymore: prune!
        if self.protein.prune(k, self.bestStability):
            self.pruneCount += 1
            if self.pruneCount % 10000 == 0:
                print(">> " + str(self.pruneCount) + " times pruned    ----     Stability: "
                      + str(self.bestStability) + " <<")
            return True
        else:
            return False

    def pruneStraight(self, k):
        """ Prunes when 3 aminoacids in a row have the same orientation

        :param k: the aminoacid currently being placed
        """

        # if a range with 3 in a row, prune!
        for j in range(1, k - 1):
            if self.foldPattern[j - 1] == self.foldPattern[j] and \
                    self.foldPattern[j + 1] == self.foldPattern[j]:
                return True

        return False

    def checkOptimum(self):
        """ Checks if there is a known optimum for the protein.
        And returns true if that optimum is reached. """

        if self.protein.optimum:
            if self.bestStability == self.protein.optimum:
                return True
            else:
                return False





