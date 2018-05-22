import csv
import numpy as np
import matplotlib.pyplot as plt
import os, sys
from os import listdir, path


class Experiment(object):
    """ Implements an experiment where statistics are visualized. """

    def __init__(self, protein, dimensions, algorithm):
        self.protein = protein
        self.dimensions = dimensions
        self.algorithm = algorithm

        self.filename = ""
        self.bestStability = []
        self.iterations = []
        self.overlap = []
        self.foundRun = []
        self.timeElapsed = []

    def setFile(self):
        """ Sets filename. """

        self.protein = "P" + str(self.protein)
        self.dimensions = str(self.dimensions) + "D"

        if self.algorithm == "R":
            self.algorithm = "Randomizer"
        elif self.algorithm == "HC":
            self.algorithm = "HillClimber"
        elif self.algorithm == "DF":
            self.algorithm = "DepthFirst"
        elif self.algorithm == "BB":
            self.algorithm = "BranchNBound"
        else:
            self.algorithm = "SimulatedAnnealing"

        self.filename = self.protein + '-' + self.dimensions + '-' + self.algorithm + '.log'

    def readFile(self):

        algorithms = ["R", "HC", "DF", "BB", "SA"]

        with open("results/" + self.filename, 'r') as logfile:
            logfile = logfile.readlines()

            for line in logfile:
                if "Stability" in line:
                    self.bestStability.append([float(s) for s in line.split() if "." in s])
                if "Iterations" in line:
                    self.iterations.append([int(s) for s in line.split() if s.isdigit()])
                if "Overlap" in line:
                    self.overlap.append([int(s) for s in line.split() if s.isdigit()])
                if "Run" in line:
                    self.foundRun.append([int(s) for s in line.split() if s.isdigit()])
                if "Time" in line:
                    self.timeElapsed.append([float(s) for s in line.split() if "." in s])

            print(self.bestStability)
            print(self.iterations)
            print(self.overlap)
            print(self.foundRun)
            print(self.timeElapsed)

    def visualize(self):
        """ Visualizes analyzed data. """

        # plt.bar(self.bestStability[0], )
        # plt.show()

        # def print
        # def stability(self, log=None):
        #     """ Analyzes stability. """
        #
        #     self.readCSVfiles(log)
        #
        #     low = 0
        #     score = []
        #
        #     with open("results/" + self.filename + ".csv", 'r') as csvfile:
        #         # determine lowest stability and adds stability to list
        #         next(csvfile)
        #         for row in csv.reader(csvfile):
        #             if row:
        #                 stability = int(abs(float(row[1])))
        #                 score.append(stability)
        #                 if stability > low:
        #                     low = stability
        #
        #     # stores occurrences of stability scores
        #     low = low + 1
        #     y = np.zeros(low)
        #     for item in score:
        #         y[item] = y[item] + 1
        #
        #     x = []
        #     start = 0
        #
        #     # create array for x axis
        #     for i in range(low):
        #         x.append(start)
        #         start = start + 1
        #
        #     # creates plot
        #     plt.bar(x, y)
        #     plt.xticks(np.arange(min(x), max(x) + 1, 1.0))
        #     plt.title(self.algorithm)
        #     plt.xlabel("Stability score")
        #     plt.ylabel("Count")
        #
        #     # add labels
        #     label = []
        #     for item in range(low):
        #         label.append(str(int(abs(float(y[item])))))
        #
        #     count = 0
        #     for bar in range(low):
        #         plt.text(x=count, y=3, s=label[count])
        #         count = count + 1
        #
        #     plt.show()









