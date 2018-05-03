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
                    self.overlapCount += 1
                    continue

                # get stability score of input protein
                self.protein.stability()

                if self.protein.stabilityScore < self.maxStability:
                    self.maxStability = self.protein.stabilityScore
                    self.bestPattern = copy.copy(self.foldPattern)
                    self.bestProtein = copy.copy(self.protein)

                    # write to csv
                    self.writer.writerow({'stability': self.protein.stabilityScore, 'foldPattern': self.foldPattern})

            else:
                self.foldPattern[k - 1] = orientation
                self.protein.fold(self.foldPattern)

                # skip if overlap detected
                if self.protein.checkOverlap(k):
                    self.overlapCount += 1
                    continue

                self.searching(k + 1)

    def printBest(self):
        """ Prints the best found solution """

        # print info
        print()
        print('DEPTH FIRST')
        print(' Maximal stability: ' + str(self.maxStability))
        print(' First found in run: ' + str(self.bestRun))
        print(' Total overlap: ' + str(self.overlapCount))
        print(' Elapsed time: ' + "{0:.4f}".format(self.elapsed))
        print()

        # print fold pattern
        print('Fold pattern: ')
        print(self.bestPattern)
        print()

        # plot protein
        self.bestProtein.visualize(('Best random solution ' + str(self.maxStability)))

