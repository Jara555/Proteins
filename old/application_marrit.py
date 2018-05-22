from classes.Protein import Protein
from algorithms.BranchNBound import BranchNBound
from algorithms.Randomizer import Randomizer
from algorithms.DepthFirst import DepthFirst
from experiment.ProteinRandomizer import ProteinRandomizer

def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    # TODO: Change these numbers per protein run !!
    maxIterations = 2
    length = 14

    randomProtein = ProteinRandomizer(length, maxIterations)
    randomProtein.run()


if __name__ == "__main__":
    main()