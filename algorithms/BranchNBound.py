from algorithms.DepthFirst import DepthFirst


class BranchNBound(DepthFirst):
    """ Subclass of DepthFirst algorithms:
    Implements Branch 'N Bound algorithms in order to efficiently fold a protein """

    def __init__(self, protein, writeCsv, maxIterations=None):
        """ Set and initiate all properties.

        :param protein: protein to be fold
        :param writeCsv: writes solutions to .csv file when ON
        :param maxIterations: stop after maxIterations
        """

        #initialize input variables
        self.maxIterations = maxIterations

        # set class properties
        DepthFirst.__init__(self, protein, writeCsv)
        self.name = "BranchNBound"

        self.protein.findBonds()

    def pruneBranchNBound(self, k):
        """ use pruning

        :param k: the aminoacid currently being placed
        :return: True if branch is pruned, False if not pruned
        """

        # prune
        if self.protein.prune(k, self.bestStability):
            self.pruneCount += 1
            return True
        else:
            return False





