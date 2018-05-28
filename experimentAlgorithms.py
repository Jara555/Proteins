import getopt
import sys
import numpy as np
import math

import matplotlib.pyplot as plt


def main(argv):
    """ Implements an experiment where statistics are visualized.
    usage:
        python Experiment.py -p <protein> -d <dimensions>

    :argument -p <protein> Protein to fold (1, 2, 3, ...)
    :argument -d <dimensions> Dimension, 2 or 3
        2: 2D
        3: 3D
    """

    # START ERROR CHECKING

    if len(sys.argv) < 3 or sys.argv[1] == "help":
        usage()
        sys.exit(1)

    # try to catch the parsers for command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hp:d:", ["help"])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(1)

    # find arguments for the input options
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt == "-p":
            # check if integer is entered or string
            try:
                protein = int(arg)
            except ValueError:
                usage()
                sys.exit()
        elif opt == "-d":
            # check if integer is entered or string
            try:
                dimensions = int(arg)
            except ValueError:
                usage()
                sys.exit()

    # check if protein number is correct
    while protein < 1 or protein > 9:
        usage()
        print("Enter a protein number (integer 1-9)")
        protein = int(input("    Protein: "))

    # check if dimension is correct
    while dimensions < 2 or dimensions > 3:
        usage()
        print("Enter dimension (2 or 3)")
        dimensions = int(input("    Dimensions: "))

    # END ERROR CHECKING

    # set list of input arguments
    algorithmsFound = []
    bestStability = []
    iterations = []
    overlap = []
    foundRun = []
    timeElapsed = []
    params = [bestStability, iterations, overlap, foundRun, timeElapsed]
    algorithms = ["Randomizer", "HillClimber", "DepthFirst", "BranchNBound", "SimulatedAnnealing"]
    xlabels = ["R", "HC", "DF", "BNB", "SA"]
    ylabels = ['Best stability (* -1)', 'Number of iterations', 'Number  of overlap', 'Found in run...',
               'Elapsed time (in s)']
    plotTitles = ['Stability', 'Number of iterations', 'Overlap', 'Found in run', 'Elapsed time']

    # START READ FROM FILE

    for i in range(len(algorithms)):
        filename = 'P' + str(protein) + '-' + str(dimensions) + 'D-' + algorithms[i] + '.log'

        try:
            with open("results/" + filename, 'r') as logfile:
                logfile = logfile.readlines()

                algorithmsFound.append(xlabels[i])

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

        except IOError:
            print(algorithms[i] + ' does not exist...')

    # END READ FROM FILE

    # START VISUALISATION

    print()
    print('----- Start comparing algorithms...')
    print()

    x = range(1, len(algorithmsFound) + 1)
    fig = plt.figure()
    fig.suptitle('Protein ' + str(protein) + ' in ' + str(dimensions) + 'D', fontsize=22)
    plotNum = 231

    # loops over params and plots data
    for i in range(len(params)):
        ax = fig.add_subplot(plotNum)
        ax.bar(x, params[i], width=0.5, color="green")
        ax.set_xticks(x)
        ax.set_xticklabels(algorithmsFound)
        ax.set_ylabel(ylabels[i])
        ax.set_yscale("log")
        ax.set_title(plotTitles[i])
        plotNum += 1

    # add legend of xlabels
    ax = fig.add_subplot(plotNum)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    ax.text(0, 0.9, 'R = Randomizer', fontsize=12)
    ax.text(0, 0.7, 'HC = Hill Climber*', fontsize=12)
    ax.text(0, 0.5, 'DF = Depth First', fontsize=12)
    ax.text(0, 0.3, 'BB = Branch \'n Bound', fontsize=12)
    ax.text(0, 0.1, 'SA = Simulated Annealing*', fontsize=12)
    ax.text(0, -0.1, '* starts with pattern extracted from randomizer')

    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')
    plt.show()

    print()
    print('----- End comparing algorithms...')
    print()

    # END VISUALISATION

    # START VISUALISATION HILL CLIMBER VERSUS SIMULATED ANNEALING

    print()
    print('----- Start Hill Climber vs Simulated Annealing')
    print()

    # initializes variables
    iteration = 100
    alg = ["HillClimber", "SimulatedAnnealing"]
    fig = plt.figure()
    fig.suptitle('Protein ' + str(protein) + ' in ' + str(dimensions) + 'D', fontsize=22)
    k = 121
    csvCount = 0

    for i in range(len(alg)):

        stability = []
        lineCount = 0
        stabilityMean = []
        count = 0

        # opens and reads csv file
        try:
            with open('results/P' + str(protein) + '-' + str(dimensions) + 'D-' + alg[i] + '.csv', 'r') as csvfile:
                csvfile = csvfile.readlines()
                csvCount += 1

                for line in csvfile:
                    lineCount += 1

                    line = line.rstrip()

                    # skip first two lines
                    if lineCount == 0 or lineCount == 1:
                        continue
                    if line:
                        line.split(',')
                        stability.append(abs(float(line.split(',')[1])))
        except IOError:
            print('csv file for ' + alg[i] + ' does not exist, craete one with \'proteins.py ... -cON\'')
            continue

        # calculates average stability per n iterations
        stabilityNum = len(stability)
        while stabilityNum > iteration:
            start = count * iteration
            end = start + iteration - 1
            stabilityMean.append(math.ceil(np.mean(stability[start:end])))

            stabilityNum -= iteration
            count += 1

        # plots average stability over time
        ax = fig.add_subplot(k)
        x = range(0, len(stabilityMean))
        ax.plot(x, stabilityMean)
        ax.set_xlabel('Iteration')
        ax.set_ylabel('Stability (* -1)')
        ax.set_title(alg[i])
        plt.ylim((0, np.amax(bestStability)))
        k += 1

    if csvCount > 0:
        plt.show()

    print()
    print('----- End Hill Climber vs Simulated Annealing')
    print()

    # START VISUALISATION HILL CLIMBER VERSUS SIMULATED ANNEALING

def usage():
    """Prints the usage of the command line arguments in the terminal """
    print()
    print("usage: python Experiment.py -p <protein> -d <dimensions>")
    print()


if __name__ == "__main__":
    main(sys.argv)








