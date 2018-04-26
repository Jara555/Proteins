import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

from new_classes.AminoAcid import AminoAcid


class Protein(object):
    """ Contains all protein info and methods  """

    def __init__(self, protein_string):
        """ Initializes all aminoacids in the protein based on input string """

        # store protein string and determine length
        self.protein_string = protein_string
        self.length_protein = len(self.protein_string)

        # create empty protein list to be filled with aminoacids
        self.protein_list = []

        # loop over chars/aminoacids in protein string
        for aa_index in range(self.length_protein):
            # create aminoacid and append to protein list
            aminoacid = AminoAcid(self.protein_string[aa_index])
            self.protein_list.append(aminoacid)

    def fold(self, folding_pattern):
        """ Folds protein according to input pattern """

        # let 1st aminoacid start at coordinates (0, 0)
        x, y = 0, 0

        self.protein_list[0].setCoordinates(x, y)

        # set all possible orientations
        plusY = [0, 1]
        minY = [0, -1]
        plusX = [1, 0]
        minX = [-1, 0]

        # set default orientation
        orientation = plusY

        # iterate over aminoacids in protein (skipping the 1st)
        for index in range(1, self.length_protein):
            # adapting orientation to folding pattern
            if folding_pattern[index] == '+X':
                orientation = plusX
            elif folding_pattern[index] == '-X':
                orientation = minX
            elif folding_pattern[index] == '+Y':
                orientation = plusY
            elif folding_pattern[index] == '-Y':
                orientation = minY

            # set new index based on (new) orientation
            x = x + orientation[0]
            y = y + orientation[1]

            # set new coordinates of aminoacid
            self.protein_list[index].setCoordinates(x, y)

    def checkOverlap(self):
        """ Checks the protein for overlap"""

        coords_list = []

        # put x and y coordinates of aminoacid in coord lists
        for aminoacid in self.protein_list:
            coords = (aminoacid.x, aminoacid.y)
            if coords in coords_list:
                self.overlap = True
            coords_list.append(coords)

        self.overlap = False

    def visualize(self, name):
        """ Prints the protein in scatter plot with lines"""

        print()
        print('Protein coordinates:')

        x = []
        y = []
        color = []

        # put x and y coordinates of aminoacid in x and y lists
        for aminoacid in self.protein_list:
            x.append(aminoacid.x)
            y.append(aminoacid.y)

            # creates color list for H = red and P = blue
            if aminoacid.type == 'H':
                color.append('red')
            else:
                color.append('blue')

        # print coordinates in terminal
        print(x)
        print(y)
        print()

        # scatter plot with line
        plt.plot(x, y, 'C3', zorder=1, lw=2, color='black')
        plt.scatter(x, y, s=50, zorder=2, color=color)
        plt.title('Protein: ' + name)
        plt.tight_layout()
        plt.axis('scaled')
        plt.ylim(min(y) - 1, max(y) + 1)
        plt.xlim(min(x) - 1, max(x) + 1)
        plt.xticks(np.arange(min(x), max(x) + 1, 1.0))
        plt.yticks(np.arange(min(y), max(y) + 1, 1.0))

        # plt.axhline(0, linestyle='--', color='gray', linewidth=0.5)  # horizontal lines
        # plt.axvline(0, linestyle='--', color='gray', linewidth=0.5)  # vertical lines

        hydrofoob = mpatches.Patch(color='red', label='H')
        polair = mpatches.Patch(color='blue', label='P')
        plt.legend(handles=[hydrofoob, polair])

        plt.show()

    def makeCombinations(self):
        """ Checks which H's can make a H-bond """

        # empty lists
        oddH = []
        evenH = []
        firstH = []
        secondH = []
        possibilities = []
        combinations = []

        # iterate through protein to find all H's and append the location to the oddH or evenH list
        for i in range(self.length_protein):
            if self.protein_list[i].type == 'H':
                if i % 2 == 0:                  # i.e. even
                    evenH.append(i)
                else:                           # i.e. odd
                    oddH.append(i)

        # rule: the second H has to be at a later position in the protein than the first H
        # make a list of first and second H's
        for i in range(len(evenH)):
            for j in range(len(oddH)):
                if evenH[i] < oddH[j]:
                    firstH.append(evenH[i])
                    secondH.append(oddH[j])
                else:
                    firstH.append(oddH[j])
                    secondH.append(evenH[i])

        # rule: the distance between two H's must be at least 2 to be able to fold
        # make a list of tuples representing the first and second H of a possible folding
        for i in range(len(firstH)):
            if (secondH[i] - firstH[i]) > 2:
                possibilities.append((firstH[i], secondH[i]))

        # get the first and the second H's of a possibility
        possibilitiesSeconds = [x[1] for x in possibilities]  # Last H of possible fold
        possibilitieFirsts = [x[0] for x in possibilities]  # First H of possible fold

        # rule: you can make a combination if the second possibilities H's come both after the first possibility H's
        for i in range(len(possibilitiesSeconds)):
            for j in range(len(possibilitieFirsts)):
                if possibilitieFirsts[j] >= possibilitiesSeconds[i]:
                    combinations.append((possibilities[i], possibilities[j]))

        # rule: you can make a combination of possibilities if the first H is smaller and the second H is larger than
        # H's of the other possible fold
        for i in range(len(possibilities)):
            for j in range(len(possibilities)):
                if possibilities[j][0] > possibilities[i][0] and possibilities[j][1] < possibilities[i][1]:
                    combinations.append((possibilities[i], possibilities[j]))

        return combinations

    def __str__(self):
        """ Prints the protein as a string """

        return self.protein_string
