import csv
import sys
import getopt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from experiments.ProteinRandomizer import ProteinRandomizer
from algorithms.BranchNBound import BranchNBound
from classes.Protein import Protein


def main(argv):
    """ Implements an experiment where protein properties are evaluated with the Branch 'n Bound algorithm

    usage:
            python proteins.py  -d <dimensions> -i <iterations> -n <number>
        :argument -d <dimensions> : dimensions to be fold in
                2 = 2D
                3 = 3D (default)
        :argument -n <number> : Number of proteins to be created
        :argument -l <length> : Lenght of the protein to be crated (default: 14)
        :argument -f <fixed h-number> : Fixed number of H's in the protein (optional)
    """

    # default values
    dimensions = 3
    number = 0
    length = 14
    fixedHNumber = None

    # ERROR CHECKING:

    # output usage in terminal
    if len(argv) < 2:
        usage()
        sys.exit(1)
    elif sys.argv[1] == "help":
        usage()
        sys.exit(1)

    # try to catch the parsers for command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hd:n:l:f:", ["help"])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(1)

    # find arguments for the input options
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt == "-d":
            dimensions = int(arg)
        elif opt == "-n":
            number = int(arg)
        elif opt == "-l":
            length = int(arg)
        elif opt == "-f":
            fixedHNumber = int(arg)

    # check if user gave the number of proteins
    if number == 0:
        usage()
        print("Enter the number of proteins to be generated: ")
        number = int(input("    Number: "))

    # END ERROR CHECKING

    # create the proteins and run them
    createProteins(number, length, fixedHNumber, dimensions)
    runAlgorithm(number, length, dimensions, fixedHNumber)

    # create the stats and visualise
    if fixedHNumber:
        createStatsListsCluster(number, length, dimensions, fixedHNumber)
    else:
        createStatsListHNumber(number, length, dimensions, fixedHNumber)


def createProteins(number, length, fixedHNumber, dimensions):
    """ create proteinStrings in textfiles

    :param length: the length of the proteins
    :param number: the amount of proteins
    :param fixedHNumber: fixed number of H's
    :param dimensions: 3 for 3D, 2 for 2D
    :return: .txt fils with proteinStrings
    """

    proteinRandomizer = ProteinRandomizer(length, number, fixedHNumber, dimensions)
    proteinRandomizer.run()


def runAlgorithm(number, length, dimensions, fixedH):
    """ runs the Branch 'n Bound algorithm on the proteins

    :param algorithmName: which algorithm should be runned
    :param number: the number of proteins to be runned
    :param length: the length of the proteins to be runned
    :param dimensions: 3 for 3D, 2 for 2D
    :param fixedH: fixed number of H's
    :return: addition of upperbound stability and found stability to .csv file (experimentProteins.csv)
    """

    writeCSV = "OFF"
    maxIterations = None

    # startnumber for saving the proteins (do not overwrite default proteins of project)
    saveNumber = 100

    # loop through the (number of) proteins created
    for i in range(number):
        proteinNumber = i + saveNumber
        resultsList = []

        # create protein from proteinstrings
        protein = Protein(dimensions, proteinNumber)

        # calculate an upperbound minStability
        protein.findBonds()
        maxStability = -1 * len(protein.bondPossibilities[0])
        counter = [0] * protein.length
        bondsMin = protein.updateBondOptions(counter, protein.bondPossibilities[0])
        minStability = -1 * len(bondsMin)

        writeResults = ("results/experimentProteins-l" + str(length) + "-d" + str(dimensions) + "-f" +
                        str(fixedH) + ".csv")

        # read results file
        with open(writeResults, 'r') as resultsReadfile:
            reader = csv.reader(resultsReadfile)
            resultsList.extend(reader)
            resultsReadfile.close()

            # add minStability to resultsList
            resultsList[i][3] = minStability

            # run BranchNBound algorithm
            algorithm = BranchNBound(protein, writeCSV, maxIterations)
            algorithm.runAlgorithm()

            # add stability to resultsList
            resultsList[i][4] = algorithm.bestStability

            # write results to csv file
            with open(writeResults, 'w', newline='') as resultsWritefile:
                writer = csv.writer(resultsWritefile)
                writer.writerows(resultsList)


def createStatsListsCluster(number, length, dimensions, fixedH):
    """ create lists of the H-count, max (H) cluster length and number of (H) clusters with associated statistics

    :param number: the number of proteins being analysed
    :param length: the length of the proteins
    :param dimensions: 2 for 2d, 3 for 3d
    :param fixedH: fixed number of H's
    :return: lists of the statistics
    """

    results = ("results/experimentProteins-l" + str(length) + "-d" + str(dimensions) + "-f" + str(fixedH) + ".csv")

    # create empty results and property lists
    resultsList = []
    stabilityHcount = [[] for _ in range(length+1)]
    stabilityMaxClusterLength = [[] for _ in range(length+1)]
    stabilityNClusters = [[] for _ in range(length+1)]

    # read results file
    with open(results, 'r') as resultsReadfile:
        reader = csv.reader(resultsReadfile)
        resultsList.extend(reader)
        resultsReadfile.close()

        # loop through the result rows and save in separate property lists
        for i in range(number):
            for j in range(length+1):
                stability = resultsList[i][4]

                if int(resultsList[i][0]) == j:
                    stabilityHcount[j].append(float(stability)*-1)

                # max cluster length
                if int(resultsList[i][1]) == j:
                    stabilityMaxClusterLength[j].append(float(stability)*-1)

                # number of clusters
                if int(resultsList[i][2]) == j:
                    stabilityNClusters[j].append(float(stability)*-1)

    visualiseStatsCluster(stabilityMaxClusterLength, stabilityNClusters, number, dimensions, length, fixedH)

def createStatsListHNumber(number, length, dimensions, fixedH):
    """ create lists of the H-count, max (H) cluster length and number of (H) clusters with associated statistics

    :param number: the number of proteins being analysed
    :param length: the length of the proteins
    :param dimensions: 2 for 2d, 3 for 3d
    :param fixedH: fixed number of H's
    :return: lists of the statistics
    """

    results = ("results/experimentProteins-l" + str(length) + "-d" + str(dimensions) + "-f" + str(fixedH) + ".csv")

    # create empty results and property lists
    resultsList = []
    stabilityHcount = [[] for _ in range(length+1)]

    # read results file
    with open(results, 'r') as resultsReadfile:
        reader = csv.reader(resultsReadfile)
        resultsList.extend(reader)
        resultsReadfile.close()

        # loop through the result rows and save in separate property lists
        for i in range(number):
            for j in range(length+1):
                stability = resultsList[i][4]

                if int(resultsList[i][0]) == j:
                    stabilityHcount[j].append(float(stability)*-1)

    visualiseStatsHCount(stabilityHcount, number, dimensions, length)

def visualiseStatsCluster(clusterLengthStatistics, clusterCountStatistics, number, dimensions, length, fixedH):
    """ visualise statistics related to H-clusters

        :param clusterLengthStatistics: list of lists of specific max cluster lengths with associated stability's
        :param clusterCountStatistics: list of lists of specific cluster counts with associated stability's
        :param number: number of proteins
        :param dimensions: 2 for 2D, 3 for 3D
        :param length: length of the proteins
        :param fixedH: fixed number of H's
        :return:
        """

    # create groups
    n = len(clusterLengthStatistics)
    ind = np.arange(n)  # the x locations for the groups

    # empty lists
    minClusterLength = []
    maxClusterLength = []
    meanClusterLength = []
    minClusterCount = []
    meanClusterCount = []
    maxClusterCount = []
    minList = [minClusterLength, minClusterCount]
    meanList = [meanClusterLength, meanClusterCount]
    maxList = [maxClusterLength, maxClusterCount]

    # plot characteristics
    ylabel = "stability (* -1)"
    xlabels = ["Longest cluster of H's", "Number of H clusters"]
    plotTitles = ["Influence of cluster length on stability", "Influence of number of H clusters on stability"]

    # calculate mean values clusterLength
    for clusterLength in clusterLengthStatistics:
        if clusterLength == []:
            minClusterLength.append(0)
            meanClusterLength.append(0)
            maxClusterLength.append(0)
        else:
            # calculate min, mean and max stability
            min = np.min(clusterLength)
            mean = np.mean(clusterLength)
            max = np.max(clusterLength)
            # set in the respective lists
            minClusterLength.append(min)
            meanClusterLength.append(mean)
            maxClusterLength.append(max)

    # calculate mean values cluster count
    for clusterCount in clusterCountStatistics:
        if clusterCount == []:
            minClusterCount.append(0)
            meanClusterCount.append(0)
            maxClusterCount.append(0)
        else:
            # calculate min, mean and max stability
            min = np.min(clusterCount)
            mean = np.mean(clusterCount)
            max = np.max(clusterCount)
            # set in the respective lists
            minClusterCount.append(min)
            meanClusterCount.append(mean)
            maxClusterCount.append(max)

    # open figure
    fig = plt.figure()
    fig.suptitle('Statistical results of ' + str(number) + ' randomly generated proteins in ' + str(dimensions) + 'D' +
                 ' with length ' + str(length) + ' and ' + str(fixedH) + 'H\'s', fontsize=14)

    # create plots
    gs1 = gridspec.GridSpec(1, 2)

    # loops over params and plots data
    for i in range(len(meanList)):
        ax = fig.add_subplot(gs1[i])
        ax.bar(ind, meanList[i], width=0.5, color="green", label="Mean stability")
        ax.set_xticks(ind)
        ax.set_xlabel(xlabels[i], fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(plotTitles[i], fontsize=12)
        ax.vlines(ind, ymin=minList[i], ymax=maxList[i], label="Range")

    plt.legend()
    plt.show()
    gs1.tight_layout(fig)

def visualiseStatsHCount(HCountStatistics, number, dimensions, length):
    """ visualise statistics

    :param HCountStatistics: list of lists of specific H-count with associated stability's
    :param number: number of proteins
    :param dimensions: 2 for 2D, 3 for 3D
    :param length: length of the proteins
    :return:
    """

    # create groups
    n = len(HCountStatistics)
    ind = np.arange(n)  # the x locations for the groups
    minHcount = []
    meanHcount = []
    maxHcount =[]

    ylabel = "stability (* -1)"
    xlabel = "Number of H's"
    plotTitle = "Influence of number of H's on stability"

    # calculate mean values H count
    for Hcount in HCountStatistics:
        if Hcount == []:
            minHcount.append(0)
            meanHcount.append(0)
            maxHcount.append(0)
        else:
            # calculate min, mean and max stability
            min = np.min(Hcount)
            mean = np.mean(Hcount)
            max = np.max(Hcount)
            # set in the respective lists
            minHcount.append(min)
            meanHcount.append(mean)
            maxHcount.append(max)

    # open figure
    fig, ax = plt.subplots()
    fig.suptitle('Statistical results of ' + str(number) + ' randomly generated proteins in ' + str(dimensions) + 'D' +
                 ' with length ' + str(length), fontsize=14)

    # plot data
    ax.bar(ind, meanHcount, width=0.5, color="green", label="Mean stability")
    ax.set_xticks(ind)
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(plotTitle, fontsize=12)
    ax.vlines(ind, ymin=minHcount, ymax=maxHcount, label="Range")

    plt.legend()
    plt.show()


def usage():
    """Prints the usage of the command line arguments in the terminal """
    print()
    print("usage: python3 proteins.py "
          "-n <number> -d <dimensions> -l <length> -f <fixed H-number>")
    print()


def printGiveNumber():
    """Tells the user to give the number of proteins to be generated """

    print("Enter the number of proteins to be generated: ")


if __name__ == "__main__":
    main(sys.argv)
