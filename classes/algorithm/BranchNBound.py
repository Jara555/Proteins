from classes.algorithm.DepthFirst import DepthFirst


class BranchNBound(DepthFirst):
    """ Subclass of DepthFirst algorithm:
    Implements Branch 'N Bound algorithm in order to efficiently fold a protein """

    def __init__(self, protein, writeCsv, maxIterations=None):
        """
        Set and initiate all properties.

        :param protein: protein to be fold
        :param writeCsv: writes solutions to .csv file when ON
        :param maxIterations: stop after maxIterations
        """

        #initialize input varialbes
        self.maxIterations = maxIterations

        # set class properties
        DepthFirst.__init__(self, protein, writeCsv)
        self.name = "BranchNBound"

        self.protein.findbonds()

    def prune(self, k):

        # go to next orientation
        if self.protein.prune(k, self.bestStability):
            self.pruneCount += 1
            return True
        else:
            return False





