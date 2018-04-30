import csv
from random import randint
from classes.Algorithms import Algorithms


class Randomizer(Algorithms):
    """ Randomize algorithm: finding best protein solution based on random patterns """

    def __init__(self, protein, iterations, writeOptions):
        Algorithms.__init__(self, protein)
        self.iterations = iterations
        self.writeOptions = writeOptions
        self.bestPattern = []
        self.maxStability = 0
        self.firstHit = 0

    def runRandomizer(self):
        """ Runs the randomizer and finds best pattern with highest stability """

        # create csv file to write output to
        write_file = ('results/random' + str(self.protein.number) + '.csv')
        with open(write_file, 'w') as csvfile:
            fieldnames = ['stability', 'foldingPattern']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            i = 0

            # iterate till the repetition amount is reached
            while i <= self.iterations:
                # get random folding pattern and fold protein according to this pattern
                self.generator()
                self.protein.fold(self.foldPattern)

                # skip if overlap detected
                if self.protein.checkOverlap(self.protein.length):
                    continue

                # get stability score of input protein
                self.protein.stability()

                # if write all is on, write every solution to csv
                if self.writeOptions == 0:
                    writer.writerow(
                        {'stability': self.protein.stabilityScore, 'foldingPattern': self.foldPattern})

                # if stability score is equal or better than max stability save new
                if self.protein.stabilityScore <= self.maxStability:
                    self.maxStability = self.protein.stabilityScore
                    self.bestPattern = self.foldPattern
                    if self.firstHit == 0:
                        self.firstHit = i

                    # if write all is off, write only best solutions to csv
                    if self.writeOptions == 1:
                        writer.writerow(
                            {'stability': self.protein.stabilityScore, 'foldingPattern': self.foldPattern})

                # next iteration
                i += 1

    def generator(self):
        """ Creates random folding pattern """

        # get size protein and create empty folding pattern list
        size = self.protein.length
        self.foldPattern = []

        # index counter for number of items in folding pattern list
        i = 0

        # iterate over required folding pattern length
        while i <= size:

            # get random orientation index
            orientation = randint(1, 4)

            # first element should always start at 0, 0
            if i == 0:
                self.foldPattern.append('0')
                i += 1
            # previous element cannot be inverse orientation
            elif orientation == 1 and self.foldPattern[i - 1] != '-X':
                self.foldPattern.append('+X')
                i += 1
            elif orientation == 2 and self.foldPattern[i - 1] != '+X':
                self.foldPattern.append('-X')
                i += 1
            elif orientation == 3 and self.foldPattern[i - 1] != '-Y':
                self.foldPattern.append('+Y')
                i += 1
            elif orientation == 4 and self.foldPattern[i - 1] != '+Y':
                self.foldPattern.append('-Y')
                i += 1

    def printBestRandom(self):
        """ Prints the best found solution """

        # print folding pattern and stability in terminal
        print()
        print('Randomizer maximal stability: ' + str(self.maxStability))
        print(self.bestPattern)
        print()

        # plot protein
        self.protein.fold(self.bestPattern)
        self.protein.visualize(('Best random solution ' + str(self.maxStability)))



