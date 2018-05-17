from classes.Algorithm import Algorithm


class BranchNBoundsub(Algorithm):
    """ Implements Branch 'N Bound algorithm in order to efficiently fold a protein """

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
        Algorithm.__init__(self, protein, writeCsv)
        self.name = "BranchNBound"
        self.protein.findbonds()

    def run(self, k):
        """ Recursive search function

        :param k: the aminoacid currently being placed
        :return: the found best folding Pattern and Stability

        """

        for orientation in self.orientations:
            self.iterations += 1

            # terminal output
            self.printProgress()

            # if a max iterations is given, don't exceed this
            if self.maxIterations:
                if self.iterations > self.maxIterations:
                    return

            # if end is reached
            if k == self.protein.length:
                # set last orientation and fold
                self.foldPattern[k - 1] = orientation
                self.protein.fold(self.foldPattern)

                # skip if overlap detected
                if self.skipOverlap(k):
                    continue

                # check for lower stability
                self.checkBest(k)

                # write to csv
                self.writeCsvRow()

            else:
                # set orienation and fold
                self.foldPattern[k - 1] = orientation
                self.protein.fold(self.foldPattern)

                # skip if overlap detected
                if self.skipOverlap(k):
                    continue

                # go to next orientation
                if self.protein.prune(k, self.bestStability):
                    self.pruneCount += 1
                    continue

                # go deeper in recursion
                self.run(k + 1)




