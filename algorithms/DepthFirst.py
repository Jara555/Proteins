from classes.Algorithm import Algorithm


class DepthFirst(Algorithm):
    """ Subclass of Algorithm:
    Implements Depth First algorithms in order to efficiently fold a protein """

    def __init__(self, protein, writeCsv, maxIterations=None):
        """ Set and initiate all properties.

        :param protein: protein to be fold
        :param writeCsv: writes solutions to .csv file when ON
        :param maxIterations: stop after maxIterations
        """

        #initialize input varialbes
        self.maxIterations = maxIterations

        # set class properties
        Algorithm.__init__(self, protein, writeCsv)
        self.name = "DepthFirst"

    def run(self, k):
        """ Recursive search function

        :param k: the aminoacid currently being placed
        :return: the found best folding Pattern and Stability
        """

        # loop over orientatoins
        for orientation in self.orientations:

            # keep track of iterations
            self.iterations += 1

            # print progress in terminal output
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
                # set orientation and fold
                self.foldPattern[k - 1] = orientation
                self.protein.fold(self.foldPattern)

                # if instance of branchNbound: use pruning
                if self.name == "BranchNBound":
                    if self.pruneBranchNBound(k):
                        continue

                # skip if overlap detected
                if self.skipOverlap(k):
                    continue

                # go deeper in recursion
                self.run(k + 1)






