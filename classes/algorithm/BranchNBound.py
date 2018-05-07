import copy
import csv
import time

from classes.Algorithms import Algorithms


class BranchNBound(Algorithms):
    """ Implements branch 'n bound algorithms in order to efficiently fold a protein """

    def __init__(self, protein):
        """ Set and initiate all properties.

        :param protein: protein you want to fold
        """

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

    def runBranchNBound(self):
        """ run the branch 'n bound algorithm

        :return: .csv file with the best folding patterns and associated stabilities
        """

        print()
        print("------------  Branch 'n Bound started ----------------")
        print()

        start = time.time()

        # open csv file
        write_file = ('results/BranchNBound' + str(self.protein.number) + '.csv')
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
        print("------------  Branch 'n Bound finished ----------------")
        print()

    def searching(self, k):
        """ Recursive search function

        :param k: the aminoacid currently being placed
        :return: calculates self.bestPattern and self.maxStability
        """

        for orientation in self.orientations:
            self.combinations += 1
            if self.combinations % 100000 == 0:
                print()
                print('BranchNBound combination: ' + str(self.combinations) + '     (stability ' + str(
                    self.maxStability) + ')' + ' (foldpattern ' + str(self.bestPattern) + ')')

            if k == self.protein.length:
                self.foldPattern[k - 1] = orientation
                self.protein.fold(self.foldPattern)

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

                # pruning: initiate
                self.protein.stability(k)
                bondOptions = copy.copy(self.protein.bondPossibilities)
                newBondOptions = []

                # remove the H-bonds that are in the protein already
                for hBond in self.protein.HBonds:
                    bondOptions.remove(hBond)

                # make a new list of only hBonds with second H after k (still potential H-bonds)
                for hBond in bondOptions:
                    if hBond[1] > k:
                        newBondOptions.append(hBond)

                # prone if potStability > maxStability
                potStability = self.protein.stabilityScore + (-1 * len(newBondOptions))  # calculate potential stability
                if potStability > self.maxStability:
                    continue

                self.searching(k + 1)

    def printBest(self):
        """ Print and visualise the best foldingPattern found.

        :return: print foldingPattern and associated stability to terminal, visualise in plot
        """

        # print info
        print()
        print("BRANCH 'N BOUND")
        print(' Maximal stability: ' + str(self.maxStability))
        print(' Total combinations: ' + str(self.combinations))
        print(' Total overlap: ' + str(self.overlapCount))
        print(' Elapsed time: ' + "{0:.4f}".format(self.elapsed))
        print()

        # print fold pattern
        print('Fold pattern: ')
        print(self.bestPattern)
        print()

        # plot protein
        self.protein.fold(self.bestPattern)
        self.protein.visualize(("Best branch 'n bound solution " + str(self.maxStability)))



