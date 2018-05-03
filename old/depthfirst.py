import csv
from classes.Protein import Protein


def main():
    """ Implements depth first algorithms in order to most efficiently fold a protein """

    # create protein of class Protein
    protein = Protein(1)
    length = protein.length

    # declare variables
    foldingPattern = ['+Y']*length
    foldingPattern[0] = '0'
    foldingPattern[1] = '+Y'

    # starting amino acid (i.e first 2 are fixed)
    k = 3
    orientations = ['+Y', '-X', '+X', '-Y']

    # open csv file
    write_file = ('results/depthFirst' + str(protein.number) + '.csv')
    with open(write_file, 'w') as csvfile:
        fieldnames = ['stability', 'foldingPattern']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        searching(foldingPattern, length, k, orientations, writer, protein)


def searching(foldingPattern, length, k, orientations, writer, protein):
    """
    :param foldingPattern: pattern protein follows on 2D grid
    :param length: number of aminoacids in protein
    :param k: aminoacid currently being placed on grid
    :param orientations: directions the aminoacid can go on the grid
    :param writer: defines to which file the results are written
    :param protein: protein being folded
    :result: .csv file with foldingPatterns and associated stability
    """
    for orientation in orientations:
        if k == length:
            foldingPattern[k - 1] = orientation
            protein.fold(foldingPattern)

            # skip if overlap detected
            if protein.checkOverlap():
                continue

            # get stability score of input protein
            protein.stability()

            # write to csv
            writer.writerow({'stability': protein.StabilityScore, 'foldingPattern': foldingPattern})

        else:
            foldingPattern[k - 1] = orientation

            proteinTemp = Protein(protein.list[0:k-1].type)
            proteinTemp.fold(foldingPattern[0:k-1])

            # skip if overlap detected
            if proteinTemp[0:k-1].checkOverlap():
                continue

            searching(foldingPattern, length, k + 1, orientations, writer, protein)


if __name__ == "__main__":
    main()
