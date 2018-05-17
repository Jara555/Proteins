import copy
import csv
import time


class Algorithm:
    """ Superclass of all protein folding algorithms:
     Contains all common properties and methods of the following subclasses:
        - Randomizer
        - HillClimber
        - DepthFirst
            - BranchNBound
    """

    def __init__(self, protein, writeCsv="OFF"):
        """
        Set and initiate all properties.

        :param protein: protein to be fold
        :param writeCsv: writes solutions to .csv file when ON

        """

        # initialize input variables
        self.protein = protein
        self.writeCsv = writeCsv

        # set class properties
        self.bestPattern = []
        self.bestStability = 0
        self.bestRun = 0
        self.overlapCount = 0
        self.iterations = 0
        self.pruneCount = 0
        self.elapsed = 0
        self.writer = None
        self.name = "Algorithm"

        # get possible orientations based on dimension
        if self.protein.dimensions == 2:
            self.orientations = ['+X', '-X', '+Y', '-Y']
        elif self.protein.dimensions == 3:
            self.orientations = ['+X', '-X', '+Y', '-Y', '+Z', '-Z']

        # starting fold pattern: first 2 aminoacids have always same coordinates ['0', '+Y', ... ]
        self.foldPattern = ['+Y'] * self.protein.length
        self.foldPattern[0] = '0'

    def runAlgorithm(self):
        """ run the algorithm

        :return: when ON returns .csv file with the best folding patterns
                 and associated stability's

        """

        print()
        print("------------   " + str(self.name) + " started   ----------------")
        print()

        # start timer
        start = time.time()

        # index used for recursive algorithms
        k = 3

        # write to csv file
        if self.writeCsv == "ON":
            write_file = ("results/P" + str(self.protein.number) + "-" + str(self.protein.dimensions) + "D-" + str(self.name) + ".csv")
            with open(write_file, 'w') as csvfile:
                fieldnames = ['run', 'stability', 'foldingPattern']
                self.writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                self.writer.writeheader()

                # start iteration
                self.run(k)

        # do not write to csv file
        elif self.writeCsv == "OFF":
            # start iteration
            self.run(k)

        # end timer
        end = time.time()
        self.elapsed = end - start

        print()
        print("------------   " + str(self.name) + " finished   ----------------")
        print()

    def checkBest(self, length=None):
        """ Checks if stability is lower and saves these new values """

        print("HALLO HIILLCLIMBER")
        print(str(self.protein.stabilityScore))

        if not length:
            length = self.protein.length

        self.protein.stability(length)

        if self.protein.stabilityScore <= self.bestStability:
            self.bestStability = self.protein.stabilityScore
            self.bestPattern = copy.copy(self.foldPattern)
            self.bestRun = self.iterations

    def skipOverlap(self, length=None):
        """ Counts and returns true if overlap """

        if not length:
            length = self.protein.length

        if self.protein.checkOverlap(length):
            self.overlapCount += 1
            return True
        else:
            return False

    def printProgress(self):
        """ Update progress in terminal output """

        if self.maxIterations == None:
            printNow = 10000
        else:
            printNow = (self.maxIterations * 0.05)

        if self.iterations % printNow == 0:
            print(str(self.name) + " iteration: " + str(self.iterations) +
                  "    ----    Stability: " + str(self.protein.stabilityScore) +
                  "    ----    Stability: " + str(self.bestStability) +
                  "    ----    Pattern: " + str(self.bestPattern))

    def writeCsvRow(self):
        """ Writes rows to CSV file """

        # if ON write every pattern to csv
        if self.writeCsv == "ON":
            self.writer.writerow(
                {'run': self.iterations, 'stability': self.protein.stabilityScore, 'foldingPattern': self.foldPattern})

    def printBest(self):
        """ Prints all running info and the best solution found """

        # print running info
        print()
        print(self.name)
        print(' Maximal stability: ' + str(self.bestStability))
        print(' Total iterations: ' + str(self.iterations))
        print(' Found in run: ' + str(self.bestRun))
        print(' Total overlap: ' + str(self.overlapCount))
        if self.name == "BranchNBound":
            print(' Total times pruned ' + str(self.pruneCount))
        print(' Elapsed time: ' + "{0:.4f}".format(self.elapsed))
        print()

        # print best found fold pattern (if found)
        print('Fold pattern: ')
        if self.bestPattern:
            print(self.bestPattern)
            print()

            # visualize protein in plot
            self.protein.fold(self.bestPattern)
            self.protein.stability(self.protein.length)
            self.protein.visualize(('Best ' + str(self.name) + ' solution ' + str(self.bestStability)))
        else:
            print('... No best pattern found...')

        print("-------------------------------------------------------------------")


