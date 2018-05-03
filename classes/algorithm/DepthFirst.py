import copy
import csv
import time

from classes.Algorithms import Algorithms


class DepthFirst(Algorithms):
    """ Implements depth first algorithms in order to most efficiently fold a protein """

    def __init__(self, protein):

        Algorithms.__init__(self, protein)

        self.foldPattern = ['+Y']*self.protein.length
        self.foldPattern[0] = '0'
        self.foldPattern[1] = '+Y'

        self.orientations = ['+Y', '-X', '+X', '-Y']
        self.bestPattern = []
        self.writer = None

        self.maxStability = 0
        self.overlapCount = 0
        self.combinations = 0
        self.elapsed = 0

    def runDepthFirst(self):

        print()
        print("------------  Depth first started ----------------")
        print()

        start = time.time()

        # open csv file
        write_file = ('results/depthFirst' + str(self.protein.number) + '.csv')
        with open(write_file, 'w') as csvfile:
            fieldnames = ['stability', 'foldPattern']
            self.writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            self.writer.writeheader()

            # recursive function
            k = 3
            self.searching(k)

        end = time.time()
        self.elapsed = end - start

        print()
        print()
        print("------------  Depth first finished ----------------")
        print()

    def searching(self, k):
        """ Recursive search function """

        for orientation in self.orientations:
            if k == self.protein.length:
                self.foldPattern[k - 1] = orientation
                self.protein.fold(self.foldPattern)

                self.combinations += 1
                if self.combinations % 10000 == 0:
                    print()
                    print('Depth first combination: ' + str(self.combinations) + '     (stability ' + str(self.maxStability) + ')')

                # skip if overlap detected
                if self.protein.checkOverlap(k):
                    self.overlapCount += 1
                    continue

                # get stability score of input protein
                self.protein.stability(k)

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
                    self.overlapCount += 1
                    continue

                self.searching(k + 1)

    def printBest(self):
        """ Print and visualise the best foldingPattern found.

        :return: print foldingPattern and associated stability to terminal, visualise in plot
        """

        # print info
        print()
        print('DEPTH FIRST')
        print(' Maximal stability: ' + str(self.maxStability))
        print(' First found in combination: ' + str(self.combinations))
        print(' Total overlap: ' + str(self.overlapCount))
        print(' Elapsed time: ' + "{0:.4f}".format(self.elapsed))
        print()

        # print fold pattern
        print('Fold pattern: ')
        print(self.bestPattern)
        print()

        # plot protein
        self.protein.fold(self.bestPattern)
        self.protein.visualize(('Best random solution ' + str(self.maxStability)))



