# the Protein class reflects a list of aminoacids
from classes.AminoAcid import AminoAcid


class Protein(object):

    # input should be a string reflecting the protein (e.g. "HPHHHPHH")
    def __init__(self, protein_string):

        self.protein_string = protein_string

        # create empty protein list to be filled with aminoacids
        self.protein_list = []

        # determine length of protein and set starting coordinates at 0
        length_protein = len(self.protein_string)
        x, y = 0, 0

        # loop over chars/aminoacids in protein string
        for aa_index in range(length_protein):

            # create aminoacid of class AminoAcid and append to protein list
            aminoacid = AminoAcid(self.protein_string[aa_index], x, y)
            self.protein_list.append(aminoacid)

            # increment y coordinate with 1 for next aminoacid
            y += 1

    # print method
    def __str__(self):
        return self.protein_string
