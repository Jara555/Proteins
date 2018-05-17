import csv
import numpy as np
import matplotlib.pyplot as plt
import os, sys
from os import listdir, path


class Experiment(object):
    """ Analyzes different aspects of the folding algorithms. """

    def __init__(self, branchBound, depthFirst, hillClimber, randomizer):
        self.branchBound = branchBound
        self.depthFirst = depthFirst
        self.hillClimber = hillClimber
        self.randomizer = randomizer

    def readCSVfiles(self):
        """ Read CSV files. """

        # list all algorithms
        path = "classes/algorithms/"
        dirs = os.listdir(path)

    def stability(self):
        """ Analyzes stability. """

        low = 0
        score = []

        with open("results/" + self.hillClimber + ".csv", 'r') as csvfile:
            # determine lowest stability and adds stability to list
            next(csvfile)
            for row in csv.reader(csvfile):
                if row:
                    stability = int(abs(float(row[1])))
                    score.append(stability)
                    if stability > low:
                        low = stability

        # stores occurrences of stability scores
        low = low + 1
        y = np.zeros(low)
        for item in score:
            y[item] = y[item] + 1

        x = []
        start = 0

        # create array for x axis
        for i in range(low):
            x.append(start)
            start = start + 1

        # creates plot
        plt.bar(x, y)
        plt.xticks(np.arange(min(x), max(x) + 1, 1.0))
        plt.title(self.hillClimber)
        plt.xlabel("Stability score")
        plt.ylabel("Count")

        # add labels
        label = []
        for item in range(low):
            label.append(str(int(abs(float(y[item])))))

        count = 0
        for bar in range(low):
            plt.text(x=count, y=y[count] + 5, s=label[count])
            count = count + 1

        plt.show()

    def overlap(self):
        """ Analyzes overlap. """

    def runningTime(self):
        """ Analyzes running time. """

    def visualize(self):
        """ Visualizes analyzed data. """









