from new_classes.Protein import Protein
from functions.findRandom import findRandom


def main():
    """ Implements random algorithms in order to most efficiently fold a protein """

    protein_nr = 1

    # open protein text file
    data_file = ('data/protein' + str(protein_nr) + '.txt')
    with open(data_file, 'r') as file:
        protein_string = file.read()

    # create protein of class Protein
    protein = Protein(protein_string)

    # print protein string
    print(protein)

    # fold protein according to pattern
    protein.fold(['+Y', '+Y', '+Y', '+Y', '+Y', '+Y', '+Y', '+Y'])

    # check for overlap
    if protein.checkOverlap():
        print("Invalid")
    else:
        print("Valid")

    # visualize protein
    protein.visualize("test")

    # randomizer!
    findRandom(protein, 1000, protein_nr, 0)


if __name__ == "__main__":
    main()