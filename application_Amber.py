from classes.Protein import Protein
from classes.algorithm.BranchNBound import BranchNBound
from classes.algorithm.Randomizer import Randomizer
from classes.algorithm.DepthFirst import DepthFirst
#from classes.algorithm.HillClimber import HillClimber
from Experiment import Experiment


def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    # TODO: Change these numbers per protein / run !!
    number = 3
    iterations = 1000
    dimensions = 2
    writeCsv = "OFF"

    # run random algorithm
    protein = Protein(number, dimensions)
    randomAlgorithm = Randomizer(protein, writeCsv, iterations)
    randomAlgorithm.runAlgorithm()

    # # run depth first algorithm
    # protein = Protein(number, dimensions)
    # depthFirstAlgorithm = DepthFirst(protein, writeCsv)
    # depthFirstAlgorithm.runAlgorithm()
    #
    # # run branch n bound algorithm
    # protein = Protein(number, dimensions)
    # branchNBoundAlgorithm = BranchNBound(protein, writeCsv)
    # branchNBoundAlgorithm.runAlgorithm()

    # run hill climber algorith
    # protein = Protein(number, dimensions)
    # hillclimberAlgirhtm = HillClimber(protein, writeCsv)
    # hillclimberAlgirhtm.runHillClimber()
    #
    # # print solutions
    randomAlgorithm.printBest()
    # depthFirstAlgorithm.printBest()
    # branchNBoundAlgorithm.printBest()
    # hillclimberAlgirhtm.printBest()

    #analyzer = Experiment("random1.2", "random1.3", "run2", "run1")
    #analyzer.stability()

    # bonds = [(13, 18), (15, 18), (17, 20), (18, 13), (18, 15), (20, 17), (3, 8), (4, 7), (7, 4), (8, 3), (4, 7), (7, 4), (30, 33), (33, 30), (4, 7), (7, 4), (20, 23), (23, 20), (16, 19), (19, 16)]
    #
    # noDouble = []
    #
    # bondsSet = list(set(bonds))
    # print(bondsSet)
    #
    # for bond in bondsSet:
    #     bond = bond[1], bond[0]
    #     print(bond)
    #     if bond in bondsSet:
    #         bondsSet.remove(bond)
    #
    #
    #
    #
    #
    # # print(bondsSet[0])
    # #
    # # print([(t[1], t[0]) for t in bondsSet])
    # #
    # # print(bondsSet[0])
    # #
    # # print(bondsSet[0])
    # #
    # # #for bond in bondsSet:
    #
    # bondsSet[0] = bondsSet[0][1], bondsSet[0][0]
    #
    # print(bondsSet)





if __name__ == "__main__":
    main()