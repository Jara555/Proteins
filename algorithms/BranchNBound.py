from algorithms.DepthFirst import DepthFirst
from algorithms.Randomizer import Randomizer


class BranchNBound(DepthFirst):
    """ Subclass of DepthFirst algorithms:
    Implements Branch 'N Bound algorithms in order to efficiently fold a protein """

    def __init__(self, protein, writeCsv, maxIterations=None):
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

        #initialize input varialbes
        self.maxIterations = maxIterations

        # set class properties
        DepthFirst.__init__(self, protein, writeCsv)
        self.name = "BranchNBound"

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
            return True
        else:
            return False





