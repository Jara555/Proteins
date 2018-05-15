from classes.Protein import Protein
from classes.algorithm.Randomizer import Randomizer
from classes.algorithm.DepthFirst import DepthFirst


def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    # TODO: Change these numbers per protein / run !!
    number = 1
    randomIterations = 100
    dimensions = 2
    writeCsv = "OFF"

    # run random algorithm
    protein = Protein(number, dimensions)
    # protein.fold(['+Y', '+Y', '-X', '-X', '+Y', '+Z', '+Z', '-X'])
    # protein.visualize("test")
    randomAlgorithm = Randomizer(protein, writeCsv, randomIterations)
    randomAlgorithm.runAlgorithm()

    # run depth first algorithm
    # protein = Protein(number)
    # depthFirstAlgorithm = DepthFirst(protein, dimensions)
    # depthFirstAlgorithm.runDepthFirst()

    # print solutions
    randomAlgorithm.printBest()
    # depthFirstAlgorithm.printBest()

if __name__ == "__main__":
    main()