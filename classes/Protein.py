import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import copy
from classes.AminoAcid import AminoAcid
from mpl_toolkits.mplot3d import Axes3D


class Protein(object):
    """ Contains all protein properties and methods """

    def __init__(self, number, dimensions):
        """ Set properties and initialize all aminoacids

        :param number: protein file number
        """

        self.number = number

        # open protein text file
        with open('data/protein' + str(number) + '.txt', 'r') as file:
            self.string = file.read()

        # initiate properties
        self.length = len(self.string)
        self.list = []
        self.listH = []
        self.HBonds = []
        self.stabilityScore = 0
        self.bondPossibilities = []
        self.dimensions  = dimensions

        # append list with aminoacids
        for aa_index in range(self.length):
            aminoacid = AminoAcid(self.string[aa_index])
            self.list.append(aminoacid)

    def fold(self, folding_pattern):
        """ Folds protein according to input pattern 2D

        :param folding_pattern: pattern followed to fold protein
        :param dimensions: 2 for 2D or 3 for 3D
        :return: coordinates of aminoacids set in the self.list
        """

        # starting coordinates
        x, y, z = 0, 0, 0

        # iterate over aminoacids in protein
        for index in range(self.length):
            # set orientation to folding pattern
            if folding_pattern[index] == '+X':
                x += 1
            elif folding_pattern[index] == '-X':
                x -= 1
            elif folding_pattern[index] == '+Y':
                y += 1
            elif folding_pattern[index] == '-Y':
                y -= 1
            elif folding_pattern[index] == '+Z':
                z += 1
            elif folding_pattern[index] == '-Z':
                z -= 1

            # set new coordinates to aminoacid
            self.list[index].setCoordinates(x, y, z)

    def checkOverlap(self, maxLength):
        """ Checks for overlap in the folded protein, both 2D and 3D

        :param maxLength: until where should the overlap be checked
        :return: True if overlap is found
        """
        coordinatesList = []

        # iterate over aminoacids
        for aminoacid in self.list[0:maxLength]:
            coordinates = (aminoacid.x, aminoacid.y, aminoacid.z)

            # if coordinates already exist: overlap detected
            if coordinates in coordinatesList:
                return True

            # else put in the list
            coordinatesList.append(coordinates)

        # if iteration finished there was no overlap
        return False

    def stability(self, maxLength):

        # set variables
        aminoType = ["H", "C"]
        stabilityEffect = [-1, -5]
        bonds = []
        #noDoubles = []
        score = 0

        # checks stability score for C-bonds and H-bonds
        for type in range(len(aminoType)):

            x = []
            y = []
            z = []
            a = []
            i = -1

            # set orientations
            orientation = self.setOrientations()

            # stores x and y coordinates of amino acids with either type "H" or "C"
            for aminoacid in self.list[0:maxLength]:
                i += 1
                if aminoacid.type == aminoType[type]:
                    x.append(aminoacid.x)
                    y.append(aminoacid.y)
                    z.append(aminoacid.z)
                    a.append(i)

            # loops over amino acids with type "H" or "C" and determines number of bonds
            for i in range(len(x)):

                # determine index of current "C" or "H" in protein
                counter = 0
                for index in range(self.length):
                    if self.string[index] == aminoType[type]:
                        if counter == i:
                            currentType = index
                            break
                        counter = counter + 1

                # iterates over orientations and checks for H- or C-bonds
                for k in range(6):
                    xbond = x[i] + orientation[k][0]
                    ybond = y[i] + orientation[k][1]
                    zbond = z[i] + orientation[k][2]

                    for n in range(len(x)):
                        if i == 0:
                            if (x[n] == xbond and y[n] == ybond and z[n] == zbond) and (
                                    self.list[currentType + 1].x != xbond or self.list[
                                    currentType + 1].y != ybond or self.list[currentType + 1].z != zbond):
                                # do the following
                                bonds.append((currentType, a[n]))
                                score = score - stabilityEffect[i]
                        elif i == len(x) - 1:
                            if (x[n] == xbond and y[n] == ybond and z[n] == zbond) and (
                                    self.list[currentType - 1].x != xbond or self.list[
                                    currentType - 1].y != ybond or self.list[currentType - 1].z != zbond):
                                # do the following
                                bonds.append((currentType, a[n]))
                                score = score - stabilityEffect[i]
                        else:
                            if (x[n] == xbond and y[n] == ybond and z[n] == zbond) and \
                                    (self.list[currentType - 1].x != xbond or self.list[
                                        currentType - 1].y != ybond or
                                        self.list[currentType - 1].z != zbond) and (
                                        self.list[currentType + 1].x != xbond or
                                        self.list[currentType + 1].y != ybond or self.list[currentType + 1].z != zbond):
                                # do the following
                                bonds.append((currentType, a[n]))
                                score = score - stabilityEffect[i]

        self.stabilityScore = score / 2

    def setOrientations(self):
        # set all possible orientations
        plusY = [0, 1, 0]
        minY = [0, -1, 0]
        plusX = [1, 0, 0]
        minX = [-1, 0, 0]
        if self.dimensions == 2:
            plusZ = [0, 0, 0]
            minZ = [0, 0, 0]
        else:
            plusZ = [0, 0, 1]
            minZ = [0, 0, -1]

        orientation = [plusX, minX, plusY, minY, plusZ, minZ]
        return orientation

    def visualize(self, name):
        """ Prints the protein in scatter plot with lines in 3D

        :param name: title of the plot
        :param dimensions: 2 for 2D or 3 for 3D
        :return: a plot
        """

        # empty lists for coordinates and color
        x, y, z = [], [], []
        color = []

        # put x and y coordinates of aminoacid in x and y lists
        for aminoacid in self.list:
            x.append(aminoacid.x)
            y.append(aminoacid.y)
            z.append(aminoacid.z)

            # creates color list for H = red, P = blue, C = orange
            if aminoacid.type == 'H':
                color.append('red')
            elif aminoacid.type == 'P':
                color.append('blue')
            elif aminoacid.type == 'C':
                color.append('orange')

        # print coordinates in terminal
        print()
        print('Protein coordinates:')
        print(x)
        print(y)
        print(z)
        print()

        fig = plt.figure()
        ax = fig.gca(projection='3d')

        # scatter plot with line
        ax.plot(x, y, z, 'C3', zorder=1, lw=2, color='black')
        ax.scatter(x, y, z, s=50, zorder=2, color=color)

        # layout
        plt.title(name)

        # limits
        ax.set_ylim(min(y) - 1, max(y) + 1)
        ax.set_xlim(min(x) - 1, max(x) + 1)
        ax.set_zlim(min(z) - 1, max(z) + 1)

        # ticks
        ax.set_xticks(np.arange(min(x), max(x) + 1, 1.0))
        ax.set_yticks(np.arange(min(y), max(y) + 1, 1.0))
        ax.set_zticks(np.arange(min(z), max(z) + 1, 1.0))

        # legend
        hydrofoob = mpatches.Patch(color='red', label='H')
        polair = mpatches.Patch(color='blue', label='P')
        cysteine = mpatches.Patch(color='orange', label='C')
        plt.legend(handles=[hydrofoob, polair, cysteine])

        plt.show()

    def findHs(self):
        """ Find the indexes of H's in the protein, both for 2D and 3D"""

        # remember H-indices
        for i in range(self.length):
            if self.list[i].type == "H":
                self.listH.append(i)

    def findHbonds(self):
        """ Checks which H's can make a H-bond, both for 2D and 3D

        :return: self.bondPossibilities a list with tuples of H-indices
        """

        # empty lists
        oddH = []
        evenH = []
        firstH = []
        secondH = []
        possibilities = []
        combinations = []

        # iterate through protein to find all H's and append the location to the oddH or evenH list
        for i in range(self.length):
            if self.list[i].type == 'H':
                if i % 2 == 0:  # i.e. even
                    evenH.append(i)
                else:  # i.e. odd
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

        self.bondPossibilities = possibilities
        self.listH = evenH + oddH

        # # get the first and the second H's of a possibility
        # possibilitiesSeconds = [x[1] for x in possibilities]  # Last H of possible fold
        # possibilitieFirsts = [x[0] for x in possibilities]  # First H of possible fold
        #
        # # rule: you can make a combination if the second possibilities H's come both after the first possibility H's
        # for i in range(len(possibilitiesSeconds)):
        #     for j in range(len(possibilitieFirsts)):
        #         if possibilitieFirsts[j] >= possibilitiesSeconds[i]:
        #             combinations.append((possibilities[i], possibilities[j]))
        #
        # # rule: you can make a combination of possibilities if the first H is smaller and the second H is larger than
        # # H's of the other possible fold
        # for i in range(len(possibilities)):
        #     for j in range(len(possibilities)):
        #         if possibilities[j][0] > possibilities[i][0] and possibilities[j][1] < possibilities[i][1]:
        #             combinations.append((possibilities[i], possibilities[j]))



    def __str__(self):
        """ Prints the protein as a string """

        return self.string

    def prune(self, maxLength, maxStability, dimensions):
        """ Check if the protein can be pruned, both for 2D and 3D

        :param maxLength: the aminoacid you reached in the protein and after which you might want to prune
        :param maxStability: the max stability found for this protein so far
        :param dimensions: 2 for 2D or 3 for 3D
        :return: True if potential stability > maxStability (attention: negative values!)
        """

        # pruning: initiate
        self.stability(maxLength, self.dimensions)
        bondOptions = copy.copy(self.bondPossibilities)
        newBondOptions = []

        # remove the H-bonds that are in the protein already
        for hBond in self.HBonds:
            bondOptions.remove(hBond)

        # make a new list of only hBonds with second H after k (still potential H-bonds)
        for hBond in bondOptions:
            if hBond[1] >= maxLength:
                newBondOptions.append(hBond)

        # for index in self.listH:
        #     w = 0
        #     q = 0
        #     for hBond in reversed(newBondOptions):
        #         # while (q != len(newBondOptions) | w != 2):
        #         if index in hBond:  # 10 should be index
        #             print(hBond)
        #             # w += 1                  # make new list with max 2 times one number (3 times for first and last)
        #             q += 1

        # prune if potential Stability >= maxStability
        potStability = self.stabilityScore + (-1 * len(newBondOptions))  # calculate potential stability
        if potStability >= maxStability:
            return True

        return False

