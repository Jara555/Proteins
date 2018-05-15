import copy
import csv
import time

from classes.Algorithms import Algorithms


class DepthFirst(Algorithms):
    """ Implements depth first algorithm in order to most efficiently fold a protein"""

    def __init__(self, protein, dimensions):
        """ Set and initiate all properties.

        :param protein: protein being folded
        :param dimensions: 2 for 2D or 3 for 3D
        """

        Algorithms.__init__(self, protein)

        self.foldPattern = ['+Y']*self.protein.length
        self.foldPattern[0] = '0'
        self.foldPattern[1] = '+Y'
        self.dimensions = dimensions

        if self.dimensions == 2:
            self.orientations = ['+Y', '-X', '+X', '-Y']
        elif self.dimensions == 3:
            self.orientations = ['+Y', '-X', '+X', '-Y', '+Z', '-Z']
        self.bestPattern = []
        self.writer = None

        self.maxStability = 0
        self.overlapCount = 0
        self.combinations = 0
        self.elapsed = 0

    def runDepthFirst(self):
        """ run the Depth First algorithm. Guarantees best solution.

        :return: .csv file with the best folding patterns and associated stability's
        """

        print()
        print("------------  Depth first started ----------------")
        print()

        start = time.time()

        # open csv file
        write_file = ('results/depthFirst' + str(self.protein.number) + '.' + str(self.dimensions) + 'D.csv')
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

    def runFastDepthFirst(self):
        """ Run the Depth First Algortihm with guaranteed best solution. Writes nothing to .csv file.

        :return: nothing
        """

        print()
        print("------------  Depth first started ----------------")
        print()

        start = time.time()

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
        """ Recursive search function

        :param k: the aminoacid currently being placed
        :return: calculates self.bestPattern and self.maxStability
        """
        for orientation in self.orientations:
            if k == self.protein.length:
                self.foldPattern[k - 1] = orientation
                self.protein.fold(self.foldPattern, self.dimensions)

                self.combinations += 1
                if self.combinations % 100000 == 0:
                    print()
                    print('Depth first combination: ' + str(self.combinations) + '     (stability ' + str(self.maxStability) + ')')
                    print(self.bestPattern)

                # skip if overlap detected
                if self.protein.checkOverlap(k):
                    self.overlapCount += 1
                    continue

                # get stability score of input protein
                self.protein.stability(k, self.dimensions)

                if self.protein.stabilityScore < self.maxStability:
                    self.maxStability = self.protein.stabilityScore
                    self.bestPattern = copy.copy(self.foldPattern)

                    # write to csv
                    # self.writer.writerow({'stability': self.protein.stabilityScore, 'foldPattern': self.foldPattern})

            else:
                self.foldPattern[k - 1] = orientation
                self.protein.fold(self.foldPattern, self.dimensions)

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
        self.protein.fold(self.bestPattern, self.dimensions)
        self.protein.visualize(('Best depth first solution ' + str(self.maxStability)), self.dimensions)



