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
        self.name = "Algorithm"

        # initialize input variables
        self.protein = protein
        self.writeCsv = writeCsv

        # set 'best' variables
        self.bestPattern = []
        self.bestStability = 0
        self.bestRun = 0

        # set counters
        self.overlapCount = 0
        self.iterations = 0
        self.pruneCount = 0
        self.elapsed = 0

        # set other
        self.writer = None
        self.randStartStability = 0

        # get orientations based on dimension
        if self.protein.dimensions == 2:
            self.orientations = ['+X', '-X', '+Y', '-Y']
        elif self.protein.dimensions == 3:
            self.orientations = ['+X', '-X', '+Y', '-Y', '+Z', '-Z']

        # starting fold pattern: first 2 aminoacids have always same coordinates ['0', '+Y', ... ]
        self.foldPattern = ['+Y'] * self.protein.length
        self.foldPattern[0] = '0'

        # determine print moment for progress printing
        if not self.maxIterations:
            self.printNow = 10000
        else:
            self.printNow = (self.maxIterations * 0.05)
            if self.printNow < 1000:
                self.printNow = 1000

    def runAlgorithm(self):
        """ run the algorithms

        :return: when ON returns .csv file with the best folding patterns
                 and associated stability's

        """

        print()
        print("------------   " + str(self.name) + " started   ----------------")
        print()

        # start timer
        start = time.time()

        # set parameter (if needed) for run method
        parameter = self.setParameter()

        # write to csv file
        if self.writeCsv == "ON":
            write_file = ("results/P" + str(self.protein.number) + "-" + str(self.protein.dimensions) + "D-" + str(self.name) + ".csv")
            with open(write_file, 'w') as csvfile:
                fieldnames = ['run', 'stability', 'foldingPattern']
                self.writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                self.writer.writeheader()

                # start iteration
                self.run(parameter)

        # do not write to csv file
        elif self.writeCsv == "OFF":
            # start iteration
            self.run(parameter)

        # end timer
        end = time.time()
        self.elapsed = end - start

        print()
        print("------------   " + str(self.name) + " finished   ----------------")
        print()

    def setParameter(self):
        """ Sets parameter for the run method.
        As default the run method uses no input parameters """
        return None

    def checkBest(self, length=None):
        """ Checks if stability is lower and saves these new values """

        # calculate stability
        self.protein.stability(length)

        # saves if better stability score
        if self.protein.stabilityScore < self.bestStability:
            self.bestStability = self.protein.stabilityScore
            self.bestPattern = copy.copy(self.foldPattern)
            self.bestRun = self.iterations

    def skipOverlap(self, length=None):
        """ Counts and returns true if overlap detected """

        # checks overlap and count
        if self.protein.checkOverlap(length):
            self.overlapCount += 1
            return True
        else:
            return False

    def printProgress(self):
        """ Update progress in terminal output """

        # check if printing time
        if self.iterations % self.printNow == 0:
            print(str(self.name) + " iteration: " + str(self.iterations) +
                  "    ----    Stability: " + str(self.bestStability))

    def writeCsvRow(self):
        """ Writes rows to CSV file """

        # if ON: write every pattern to csv
        if self.writeCsv == "ON":
            self.writer.writerow(
                {'run': self.iterations, 'stability': self.protein.stabilityScore, 'foldingPattern': self.foldPattern})

    def printBest(self):
        """ Prints all running info and the best solution found """

        # print running info
        print()
        print(self.name)
        if self.name == "SimulatedAnnealing":
            print(' Started with random stablity of: ' + str(self.randStartStability))
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

        # write log to csv
        self.writeLog()

    def writeLog(self):
        """ Writes all running info and the best solution found to .csv file"""

        # write to csv file
        write_file = ("results/P" + str(self.protein.number) + "-"
                      + str(self.protein.dimensions) + "D-" + str(self.name) + ".log")
        with open(write_file, 'w') as logfile:
            writer = csv.writer(logfile)

            # HEADER
            writer.writerow({'_________________________________________________'})
            writer.writerow({})
            writer.writerow({str(self.name)})
            writer.writerow({'_________________________________________________'})
            writer.writerow({})
            writer.writerow({'Protein ' + str(self.protein.number) + ': ' + str(self.protein.string) })
            writer.writerow({})

            # PROTEIN INFO
            writer.writerow({'Best Stability Found: ' + str(self.bestStability)})
            writer.writerow({})
            writer.writerow({'Best Folding Pattern Found:'})
            writer.writerow({str(self.bestPattern)})
            writer.writerow({})

            # RUNNING INFO
            if self.name == "SimulatedAnnealing":
                writer.writerow({'Started with stability of:  ' + str(self.randStartStability)})
            writer.writerow({'Total Iterations:  ' + str(self.iterations)})
            writer.writerow({'Total Overlap:     ' + str(self.overlapCount)})
            if self.name == "BranchNBound":
                writer.writerow({'Total pruned:     ' + str(self.pruneCount)})
            writer.writerow({'Found in Run:      ' + str(self.bestRun)})

            writer.writerow({'Elapsed Time:      ' + "{0:.4f}".format(self.elapsed)})






