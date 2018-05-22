import csv
import numpy as np
import matplotlib.pyplot as plt
import os, sys
from os import listdir, path
import getopt
import sys


def main(argv):
    """ Implements an experiment where statistics are visualized.
    """

    # START ERROR CHECKING

    if len(sys.argv) < 3 or sys.argv[1] == "help":
        usage()
        sys.exit(1)

    # END ERROR CHECKING

    protein = 'P' + sys.argv[1][2]
    dimensions = sys.argv[2][2] + 'D'
    algorithms = ["Randomizer", "HillClimber", "DepthFirst", "BranchNBound"]

    bestStability = []
    iterations = []
    overlap = []
    foundRun = []
    timeElapsed = []
    params = [bestStability, iterations, overlap, foundRun, timeElapsed]

    # START READ FROM FILE

    for i in range(len(algorithms)):
        filename = protein + '-' + dimensions + '-' + algorithms[i] + '.log'

        with open("results/" + filename, 'r') as logfile:
            logfile = logfile.readlines()

            for line in logfile:
                if "Stability" in line:
                    temp = [float(s) for s in line.split() if "." in s]
                    bestStability.append(abs(temp[0]))
                if "Iterations" in line:
                    temp = [int(s) for s in line.split() if s.isdigit()]
                    iterations.append(temp[0])
                if "Overlap" in line:
                    temp = [int(s) for s in line.split() if s.isdigit()]
                    overlap.append(temp[0])
                if "Run" in line:
                    temp = [int(s) for s in line.split() if s.isdigit()]
                    foundRun.append(temp[0])
                if "Time" in line:
                    temp = [float(s) for s in line.split() if "." in s]
                    timeElapsed.append(temp[0])

    # END READ FROM FILE

    # START VISUALISATION

    x = [1, 2, 3, 4]
    fig = plt.figure()
    fig.suptitle('Protein ' + protein[1] + ' ' + dimensions, fontsize=18)
    plotNum = 231

    # loops over params and plots data
    for i in range(4):
        ax = fig.add_subplot(plotNum)
        ax.bar(x, params[i], width=0.5, color="green")
        ax.set_xticks(x)
        ax.set_xticklabels(algorithms)
        ax.set_ylabel('Best stability')
        plotNum += 1

    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')
    plt.show()

    # END VISUALISATION


def usage():
    """Prints the usage of the command line arguments in the terminal """
    print()
    print("usage: python3 Experiment.py "
          " -p <protein> -d <dimensions>")
    print()


if __name__ == "__main__":
    main(sys.argv)








