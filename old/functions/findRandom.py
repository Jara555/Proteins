from old.functions.randomizer import randomizer
from old.functions.foldprotein import foldprotein
from old.functions.overlap import overlap
from old.functions import visualize
from old.functions import stability
import csv


def findRandom(protein, repeats, protein_nr, write):
    """
    Finds best random folded protein
            protein = list of aminoacids
            repeats = times to try random solutions
            protein_nr = number of the protein as saved in data directory
            write =
                0: write all solutions to csv output file
                1: write only the best solutions to csv output file
    """

    # get protein list and string
    protein_list = protein.list
    protein_string = protein.string

    # declare variables
    i = 0
    max_stability = 0
    best_pattern = []

    # create csv file to write output to
    write_file = ('results/run' + str(protein_nr) + '.csv')
    with open(write_file, 'w') as csvfile:
        fieldnames = ['number', 'stability', 'foldingPattern']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        # iterate till the repetition amount is reached
        while i <= repeats:
            # get random folding pattern and fold protein according to this pattern
            folding_pattern = randomizer(protein_string)
            folded_protein = foldprotein(protein_list, folding_pattern)

            # skip if overlap detected
            if overlap(folded_protein):
                continue

            # get stability score of input protein
            stability_score = stability(folded_protein)

            # if write all is on, write every solution to csv
            if write == 0:
                writer.writerow({'number': i, 'stability': stability_score, 'foldingPattern': folding_pattern})

            # if stability score is equal or better than max stability save new
            if stability_score <= max_stability:
                max_stability = stability_score
                best_pattern = folding_pattern
                best_i = i

                # if write all is off, write only best solutions to csv
                if write == 1:
                    writer.writerow({'number': i, 'stability': stability_score, 'foldingPattern': folding_pattern})

            # next iteration
            i += 1

    # print folding pattern and stability in terminal
    print()
    print('The highest stability found was: ' + str(max_stability) + ' (run ' + str(best_i) + ')')
    print()
    print(best_pattern)
    print()

    # plot protein
    best_protein = foldprotein(protein_list, best_pattern)
    visualize(best_protein, ('Best random solution ' + str(max_stability)))


