#!/usr/bin/python3
import getopt
import sys

from algorithms.SimulatedAnnealing import SimulatedAnnealing
from classes.Protein import Protein
from algorithms.BranchNBound import BranchNBound
from algorithms.DepthFirst import DepthFirst
from algorithms.HillClimber import HillClimber
from algorithms.Randomizer import Randomizer


def main(argv):
    """ Implements random algorithms in order to most efficiently fold a protein
        usage:
            python3 proteins.py -a <algorithm> -p <protein> -d <dimensions> -i <iterations> -c <csv>

        :argument -a <algorithm> : First letters of algorithm names
                R = Randomizer
                HC = HillClimber
                DF = DepthFirst
                BB = BranchNBound
                SA = SimulatedAnnealing
        :argument -p <protein> : Protein to be fold, can be a number (1/2/3/ .. 9) or a string (types H / P / C)
        :argument -d <dimensions> : dimensions to be fold in
                2 = 2D
                3 = 3D (default)
        :argument -i <iterations> : Maximal iterations to be run (required for R and HC, optional for DF and BB)
        :argument -c <csv> : Write solutions to .csv file
                ON = write
                OFF = not write (default)
    """

    # set variables
    algorithmNames = ["R", "HC", "DF", "BB", "SA"]
    aminoAcidTypes = ["H", "P", "C"]

    # initial values
    algorithmName = None
    proteinNumber = 0
    proteinString = None

    # default values
    dimensions = 3
    maxIterations = None
    writeCsv = "OFF"

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
        opts, args = getopt.getopt(sys.argv[1:], "ha:p:d:i:c:", ["help"])
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
            algorithmName = arg
        elif opt == "-p":
            # check if integer is entered or string
            try:
                proteinNumber = int(arg)
            except ValueError:
                proteinString = arg.upper()
        elif opt == "-d":
            dimensions = int(arg)
        elif opt == "-i":
            maxIterations = int(arg)
        elif opt == "-c":
            writeCsv = arg.capitalize()

    # check if algorithm name is correct
    while algorithmName not in algorithmNames:
        usage()
        printAlgorithmOptions()
        algorithmName = input("    Algorithm: ").upper()

    # check if protein is correct
    if not proteinString:
        # check if protein number is correct
        while proteinNumber < 1 or proteinNumber > 9:
            usage()
            print("Enter a protein number (integer 1-9)")
            proteinNumber = int(input("    Protein: "))
    else:
        # check if valid amino acid types in protein string
        for i in proteinString:
            if i not in aminoAcidTypes:
                printAminoAcidOptions()

    # check if max iterations is entered (if needed)
    if not maxIterations and (algorithmName == "R" or algorithmName == "HC" or algorithmName == "SA"):
        # go with default value
        if proteinNumber < 3:
            maxIterations = 1000
        else:
            maxIterations = 100000

    # END ERROR CHECKING

    # initialize protein
    protein = Protein(dimensions, proteinNumber, proteinString)

    # Randomizer
    if algorithmName == "R":
        algorithm = Randomizer(protein, writeCsv, maxIterations)

    # HillClimber
    elif algorithmName == "HC":
        # random start pattern
        randomAlgorithm = Randomizer(protein, writeCsv, maxIterations=100)
        randomAlgorithm.runAlgorithm()
        startPattern = randomAlgorithm.bestPattern
        algorithm = HillClimber(protein, writeCsv, maxIterations, startPattern)

    # DepthFirst
    elif algorithmName == "DF":
        algorithm = DepthFirst(protein, writeCsv, maxIterations)

    # BranchNBound
    elif algorithmName == "BB":
        algorithm = BranchNBound(protein, writeCsv, maxIterations)

    elif algorithmName == "SA":
        randomAlgorithm = Randomizer(protein, writeCsv, maxIterations=100)
        randomAlgorithm.runAlgorithm()
        startPattern = randomAlgorithm.bestPattern
        algorithm = SimulatedAnnealing(protein, writeCsv, maxIterations, startPattern)

    else:
        print("Error: Input algorithm could not be found")
        sys.exit(2)

    # run the created instance of the Algorithm class and print the best solution
    algorithm.runAlgorithm()
    algorithm.printBest()


def usage():
    print()
    print("usage: python3 proteins.py "
          "-a <algorithm> -p <protein> -d <dimensions> -i <iterations> -c <csv>")
    print()

def printAlgorithmOptions():
    print("Algorithm input options:")
    print("     R  - Randomizer")
    print("     HC - HillClimber")
    print("     DF - DepthFirst")
    print("     BB - BranchNBound")
    print()
    print("Enter an algorithm")

def printAminoAcidOptions():
    usage()
    print("Possible amino acid types for the protein string are:")
    print("     H / P / C")
    print()
    sys.exit(4)


if __name__ == "__main__":
    main(sys.argv)
