import csv
import sys
import getopt
from experiment.ProteinRandomizer import ProteinRandomizer
from algorithms.BranchNBound import BranchNBound
from classes.Protein import Protein


def main(argv):
    """ Implements an experiment where protein properties are evaluated.

    usage:
            python proteins.py -a <algorithm> -d <dimensions> -i <iterations> -n <number>

        :argument -a <algorithm> : First letters of algorithm names
                R = Randomizer
                HC = HillClimber
                SA = SimulatedAnnealing
                DF = DepthFirst
                BB = BranchNBound (default)
        :argument -d <dimensions> : dimensions to be fold in
                2 = 2D
                3 = 3D (default)
        :argument -i <iterations> : Maximal iterations to be run (required for R and HC, optional for DF and BB)
        :argument -n <number> : Number of proteins to be created (default: 100)
        :argument -l <length> : Lenght of the protein to be crated (default: 14)
    """

    # set lists for input arguments
    algorithmNames = ["R", "HC", "DF", "BB", "SA"]

    # default values
    dimensions = 3
    maxIterations = None
    randIterations = 1000
    number = 100
    length = 14
    algorithmName = "BB"
    algorithm = None

    # ERROR CHECKING:

    # output usage in terminal
    if len(argv) < 1:
        usage()
        sys.exit(1)
    elif sys.argv[1] == "help":
        usage()
        sys.exit(1)

    # try to catch the parsers for command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ha:d:i:n:l:", ["help"])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(1)

    # find arguments for the input options
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt == "-a":
            algorithmName = arg.upper()
        elif opt == "-d":
            dimensions = int(arg)
        elif opt == "-i":
            maxIterations = int(arg)
        elif opt == "-n":
            number = int(arg)
        elif opt == "-l":
            length = int(arg)

    # check if algorithm name is correct
    while algorithmName not in algorithmNames:
        printAlgorithmOptions()
        algorithmName = input("    Algorithm: ").upper()

    # check if max iterations is entered and set if needed
    if not maxIterations and (algorithmName == "R" or algorithmName == "HC"):
        # go with default value
        maxIterations = 1000

    # set amount of randomizer iterations to generate a starting pattern for HC and SA
    # if algorithmName == "SA" or algorithmName == "HC":
    #     if proteinNumber in [1, 2, 3]:
    #         randIterations = 10000
    #     elif proteinNumber in [4, 6, 7]:
    #         randIterations = 100000
    #     elif proteinNumber in [5, 8, 9]:
    #         randIterations = 1000000

    # END ERROR CHECKING

    createProteins(length, number)
    runAlgorithm(algorithmName, number, dimensions, maxIterations)


def createProteins(length, number):
    """ create proteinStrings in textfiles

    :param length: the length of the proteins
    :param number: the amount of proteins
    :return: .txt fils with proteinStrings
    """
    proteinRandomizer = ProteinRandomizer(length, number)
    proteinRandomizer.run()


def runAlgorithm(algorithmName, number, dimensions, maxIterations):
    """

    :param algorithmName: which algorithm should be runned
    :param number: the number of proteins to be runned
    :param dimensions: 3 for 3D, 2 for 2D
    :param maxIterations: how many iterations should the algorithm do
    :return: addition of upperbound stability and found stability to .csv file (experimentProteins.csv)
    """

    writeCSV = "OFF"

    # loop through the (number of) proteins created
    for i in range(number):
        proteinNumber = i + 100
        resultsList = []

        # create protein from proteinstrings
        protein = Protein(dimensions, proteinNumber)

        # calculate an upperbound (maxStability high upperbound, minStability lower (upper?)bound) of stability
        protein.findBonds()
        maxStability = -1 * len(protein.bondPossibilities[0])
        counter = [0] * protein.length
        bondsMin = protein.updateBondOptions(counter, protein.bondPossibilities[0])
        minStability = -1 * len(bondsMin)

        write_results = ("results/experimentProteins" + ".csv")

        # read results file
        with open(write_results, 'r') as resultsReadfile:
            reader = csv.reader(resultsReadfile)
            resultsList.extend(reader)
            resultsReadfile.close()

            # add minStability to resultsList
            resultsList[i][3] = minStability

            # run BranchNBound algorithm
            if algorithmName == "BB":
                algorithm = BranchNBound(protein, writeCSV, maxIterations)
                algorithm.runAlgorithm()

                # add stability to resultsList
                resultsList[i][4] = algorithm.bestStability

                # write results to csv file
                with open(write_results, 'w', newline='') as resultsWritefile:
                    writer = csv.writer(resultsWritefile)
                    writer.writerows(resultsList)


def usage():
    """Prints the usage of the command line arguments in the terminal """
    print()
    print("usage: python3 proteins.py "
          "-a <algorithm> -p <protein> -d <dimensions> -i <iterations> -n <number>")
    print()

def printAlgorithmOptions():
    """Prints the algorithm input options in the terminal """

    usage()
    print("Algorithm input options:")
    print("     R  - Randomizer")
    print("     HC - HillClimber")
    print("     SA - SimulatedAnnealing")
    print("     DF - DepthFirst")
    print("     BB - BranchNBound")
    print()
    print("Enter an algorithm")

if __name__ == "__main__":
    main(sys.argv)
