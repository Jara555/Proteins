import csv
import time


class Algorithm:
    """ Contains all protein folding algorithms """

    def __init__(self, protein, writeCsv="OFF"):

        # initialize input variables
        self.protein = protein
        self.writeCsv = writeCsv

        # set class properties
        self.bestPattern = []
        self.bestStability = 0
        self.bestRun = 0
        self.overlapCount = 0
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
        """
        run the algorithm

        :return: when ON returns .csv file with the best folding patterns
                 and associated stability's

        """

        print()
        print("------------   " + str(self.name) + " started   ----------------")
        print()

        # start timer
        start = time.time()

        # write to csv file
        if self.writeCsv == "ON":
            write_file = ('results/' + str(self.name) + str(self.protein.number) + '-' + str(self.protein.dimensions) + 'D.csv')
            with open(write_file, 'w') as csvfile:
                fieldnames = ['run', 'stability', 'foldingPattern']
                self.writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                self.writer.writeheader()

                # start iteration
                self.pickAlgorithm()

        # do not write to csv file
        elif self.writeCsv == "OFF":
            # start iteration
            self.pickAlgorithm()

        # end timer
        end = time.time()
        self.elapsed = end - start

        print()
        print("------------   " + str(self.name) + " finished   ----------------")
        print()

    def pickAlgorithm(self):

        if self.name == "Randomizer":
            self.runRandomizer()
        elif self.name == "Branch N Bound" \
                or self.name == "Depth First":
            k = 3
            self.searching(k)

    def printBest(self):
        """ Prints all running info and the best solution found """

        # print running info
        print()
        print(self.name)
        print(' Maximal stability: ' + str(self.bestStability))
        print(' Total iterations: ' + str(self.iterations))
        print(' Found in run: ' + str(self.bestRun))
        print(' Total overlap: ' + str(self.overlapCount))
        if self.name == "Branch N Bound":
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
            self.protein.visualize(('Best ' + str(self.name) + ' solution ' + str(self.bestStability)))
        else:
            print('... No best pattern found...')


