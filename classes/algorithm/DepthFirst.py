import csv
from classes.Algorithms import Algorithms


class DepthFirst(Algorithms):
    """ Implements depth first algorithms in order to most efficiently fold a protein """

    def __init__(self, protein):
        Algorithms.__init__(self, protein)
        self.foldingPattern = ['+Y']*self.protein.length
        self.foldingPattern[0] = '0'
        self.foldingPattern[1] = '+Y'
        self.orientations = ['+Y', '-X', '+X', '-Y']

    def runDepthFirst(self):

        # open csv file
        write_file = ('results/depthFirst' + str(self.protein.number) + '.csv')
        with open(write_file, 'w') as csvfile:
            fieldnames = ['stability', 'foldingPattern']
            self.writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            self.writer.writeheader()

            # recursive function
            k = 3
            self.searching(k)

    def searching(self, k):
        """
        :param foldingPattern: pattern protein follows on 2D grid
        :param length: number of aminoacids in protein
        :param k: aminoacid currently being placed on grid
        :param orientations: directions the aminoacid can go on the grid
        :param writer: defines to which file the results are written
        :param protein: protein being folded
        :result: .csv file with foldingPatterns and associated stability
        """
        for orientation in self.orientations:
            if k == self.protein.length:
                self.foldingPattern[k - 1] = orientation
                self.protein.fold(self.foldingPattern)

                # skip if overlap detected
                if self.protein.checkOverlap():
                    continue

                # get stability score of input protein
                self.protein.stability()

                # write to csv
                self.writer.writerow({'stability': self.protein.stabilityScore, 'foldingPattern': self.foldingPattern})

            else:
                self.foldingPattern[k - 1] = orientation

                self.protein.fold(self.foldingPattern)

                # skip if overlap detected
                if self.protein.checkOverlap():
                    continue

                self.searching(k + 1)

