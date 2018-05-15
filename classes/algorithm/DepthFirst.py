import copy

from classes.Algorithm import Algorithm


class DepthFirst(Algorithm):
    """ Implements Depth First algorithm in order to efficiently fold a protein """

    def __init__(self, protein, writeCsv):
        """
        Set and initiate all properties.

        :param protein: protein to be fold
        :param writeCsv: writes solutions to .csv file when ON

        """

        # set class properties
        Algorithm.__init__(self, protein, writeCsv)
        self.name = "Depth First"
        self.iterations = 0

    def searching(self, k):
        """ Recursive search function

        :param k: the aminoacid currently being placed
        :return: the found best folding Pattern and Stability

        """
        for orientation in self.orientations:
            self.iterations += 1

            # terminal output
            if self.iterations % 100000 == 0:
                print()
                print('BranchNBound combination: ' + str(self.iterations) + '     (stability ' + str(
                    self.bestStability) + ')' + ' (foldpattern ' + str(self.bestPattern) + ')')

            if k == self.protein.length:
                self.foldPattern[k - 1] = orientation
                self.protein.fold(self.foldPattern)

                # skip if overlap detected
                if self.protein.checkOverlap(k):
                    self.overlapCount += 1
                    continue

                # get stability score of input protein
                self.protein.stability(k)

                if self.protein.stabilityScore < self.bestStability:
                    self.bestStability = self.protein.stabilityScore
                    self.bestPattern = copy.copy(self.foldPattern)
                    self.bestRun = self.iterations

                # write to csv
                if self.writeCsv == "ON":
                    self.writer.writerow({'run': self.iterations, 'stability': self.protein.stabilityScore, 'foldingPattern': self.foldPattern})

            else:
                self.foldPattern[k - 1] = orientation
                self.protein.fold(self.foldPattern)

                # skip if overlap detected
                if self.protein.checkOverlap(k):
                    self.overlapCount += 1
                    continue

                self.searching(k + 1)




