from old.old_classes import Protein
from functions.findRandom import findRandom


def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    protein_nr = 2

    # open protein text file
    data_file = ('data/protein' + str(protein_nr) + '.txt')
    with open(data_file, 'r') as file:
        protein_string = file.read()

    # print protein
    print()
    print(protein_string)
    print()

    # create protein of class Protein
    protein = Protein(protein_string)

    # find best random solution
    # input: protein, number of repetitions, protein number, csv write all (0) or write best only (1)

    findRandom(protein, 1000, protein_nr, 0)


if __name__ == "__main__":
    main()
