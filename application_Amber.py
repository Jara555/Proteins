from classes.Protein import Protein
from classes.algorithm.HillClimber import HillClimber


def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    number = 2
    iterations = 10
    writeOptions = 1

    # run random algorithm
    protein = Protein(number)
    bestHillClimber = HillClimber(protein, iterations, writeOptions)
    bestHillClimber.runHillClimber()
    bestHillClimber.printBestHill()


if __name__ == "__main__":
    main()