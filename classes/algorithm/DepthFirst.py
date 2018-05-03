import copy
import csv
from classes.Algorithms import Algorithms


class DepthFirst(Algorithms):
    """ Implements depth first algorithms in order to most efficiently fold a protein """

    def __init__(self, protein):
        Algorithms.__init__(self, protein)
        self.foldPattern = ['+Y']*self.protein.length
        self.foldPattern[0] = '0'
        self.foldPattern[1] = '+Y'
        self.orientations = ['+Y', '-X', '+X', '-Y']
        self.maxStability = 0
        self.bestPattern = []

    def runDepthFirst(self):

        # open csv file
        write_file = ('results/depthFirst' + str(self.protein.number) + '.csv')
        with open(write_file, 'w') as csvfile:
            fieldnames = ['stability', 'foldPattern']
            self.writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            self.writer.writeheader()

            # recursive function
            k = 3
            self.searching(k)

    def searching(self, k):
        """ Recursive search function """

        for orientation in self.orientations:
            if k == self.protein.length:
                self.foldPattern[k - 1] = orientation
                self.protein.fold(self.foldPattern)

                # skip if overlap detected
                if self.protein.checkOverlap(k):
                    continue

                # get stability score of input protein
                self.protein.stability()

                if self.protein.stabilityScore < self.maxStability:
                    self.maxStability = self.protein.stabilityScore
                    self.bestPattern = copy.copy(self.foldPattern)

                # write to csv
                self.writer.writerow({'stability': self.protein.stabilityScore, 'foldPattern': self.foldPattern})

            else:
                self.foldPattern[k - 1] = orientation
                self.protein.fold(self.foldPattern)

                # skip if overlap detected
                if self.protein.checkOverlap(k):
                    continue

                self.searching(k + 1)

    def printBestDepth(self):
        """ Print and visualise the best foldingPattern found.

        :return: print foldingPattern and associated stability to terminal, visualise in plot
        """

        print()
        print('Depth first maximal stability: ' + str(self.maxStability))
        print(self.bestPattern)
        print()

        # plot protein
        self.protein.fold(self.bestPattern)
        self.protein.visualize(('Best random solution ' + str(self.maxStability)))

