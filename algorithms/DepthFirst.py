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

        # set starting index of iteration
        self.startIndex = 3

    def run(self, k):
        """ Recursive search function

        :param k: the aminoacid currently being placed
        :return: the found best folding Pattern and Stability
        """

        # loop over orientations
        for orientation in self.orientations:

            # if a max iterations is given, don't exceed this
            if self.maxIterations:
                if self.iterations > self.maxIterations:
                    return

            # if an optimal stability is known, stop when reached
            if self.checkOptimum():
                return

            # if end of protein is reached
            if k == self.protein.length:

                # keep track of iterations
                self.iterations += 1

                # print progress in terminal output
                self.printProgress()

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

                # prune on basis of stability
                if self.pruneStability(k):
                    continue

                # skip if overlap detected
                if self.skipOverlap(k):
                    continue

                # go deeper in recursion
                self.run(k + 1)

    def pruneStability(self, k):
        """ Non-functional method for depth first class.
        Overriden in the Branch N bound class in order to prune based on stability """
        return False

    def getParameter(self):
        """ Returns the starting parameter for variable k,
        reflecting the starting position of the aminoacid in the recursive function
        :return k: aminoacid location """
        return self.startIndex

    def checkOptimum(self):
        """ Non-functional method for depth first class.
        Overriden in the Branch N bound class in order to stop when optimum is reached """
        return False








