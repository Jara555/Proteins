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
        :param dimensions: 2 for 2D or 3 for 3D
        """

        # initialize input variables
        self.number = number
        self.dimensions = dimensions

        # open protein text file
        with open('data/protein' + str(number) + '.txt', 'r') as file:
            self.string = file.read()

        # set class properties
        self.length = len(self.string)
        self.list = []
        self.listH = []
        self.bonds = []
        self.stabilityScore = 0
        self.bondPossibilities = []

        # append list with aminoacids
        for aa_index in range(self.length):
            aminoacid = AminoAcid(self.string[aa_index])
            self.list.append(aminoacid)

        # first 2 aminoacids have always same coordinates ['0', 'Y', ... ]
        self.list[0].setCoordinates(0, 0, 0)
        self.list[1].setCoordinates(0, 1, 0)

    def fold(self, folding_pattern):
        """ Folds protein according to input pattern 2D

        :param folding_pattern: pattern followed to fold protein
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
        bondsH = []
        bondsC = []
        self.bonds = [bondsH, bondsC]
        noDoublesH = []
        noDoublesC = []
        noDoubles = [noDoublesH, noDoublesC]
        score = 0

        # checks stability score for C-bonds and H-bonds
        for type in range(2):

            x = []
            y = []
            z = []
            a = []
            i = -1
            currentType = 0

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
                            if (x[n] == xbond and y[n] == ybond and z[n] == zbond) and (self.list[currentType + 1].x != xbond or self.list[currentType + 1].y != ybond or self.list[currentType + 1].z != zbond):
                                self.bonds[type].append((currentType, a[n]))
                                score = score + stabilityEffect[type]
                        elif i == len(x) - 1:
                            if (x[n] == xbond and y[n] == ybond and z[n] == zbond) and (
                                    self.list[currentType - 1].x != xbond or self.list[
                                    currentType - 1].y != ybond or self.list[currentType - 1].z != zbond):
                                # do the following
                                self.bonds[type].append((currentType, a[n]))
                                score = score + stabilityEffect[type]
                        else:
                            if (x[n] == xbond and y[n] == ybond and z[n] == zbond) and \
                                    (self.list[currentType - 1].x != xbond or self.list[currentType - 1].y != ybond or
                                     self.list[currentType - 1].z != zbond) and (self.list[currentType + 1].x != xbond or
                                                        self.list[currentType + 1].y != ybond or
                                                                              self.list[currentType + 1].z != zbond):
                                self.bonds[type].append((currentType, a[n]))
                                score = score + stabilityEffect[type]

        # remove double counted bonds
        for type in range(len(aminoType)):
            print(type)
            for bond in self.bonds[type]:
                print(bond)
                if bond[0] < bond[1]:
                    noDoubles[type].append((bond[0], bond[1]))

        self.bonds = noDoubles
        self.stabilityScore = score / 2

    def setOrientations(self):
        # set all possible orientations
        plusY = [0, 1, 0]
        minY = [0, -1, 0]
        plusX = [1, 0, 0]
        minX = [-1, 0, 0]
        plusZ = [0, 0, 1]
        minZ = [0, 0, -1]

        orientation = [plusX, minX, plusY, minY, plusZ, minZ]

        return orientation

    def visualize(self, name):
        """ Prints the protein in scatter plot with lines in 3D

        :param name: title of the plot
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

    def findbonds(self):
        """ Checks which H/C's can make a HH-bond and CC-bond, both for 2D and 3D

        :return: self.bondPossibilities a list with 2 lists of tuples of H-indices en C-indices
        """

        aminoType = ["H", "C"]

        # empty lists combining for two aminotypes
        odd = [[], []]
        even = [[], []]
        first = [[], []]
        second = [[], []]
        possibilities = [[], []]
        self.possibilities = [[], []]

        for type in range(len(aminoType)):

            # iterate through protein to find all H's and append the location to the oddH or evenH list
            for i in range(self.length):
                if self.list[i].type == aminoType[type]:
                    if i % 2 == 0:  # i.e. even
                        even[type].append(i)
                    else:  # i.e. odd
                        odd[type].append(i)

            # rule: the second H/C of a bond has to be at a later position in the protein than the first H/C
            # make a list of first and second H/C's
            for i in range(len(even[type])):
                for j in range(len(odd[type])):
                    if even[type][i] < odd[type][j]:
                        first[type].append(even[type][i])
                        second[type].append(odd[type][j])
                    else:
                        first[type].append(odd[type][j])
                        second[type].append(even[type][i])

            # rule: the distance between two H/C's must be at least 2 to be able to fold
            # make a list of tuples representing the first and second H/C of a possible folding
            for i in range(len(first[type])):
                if (second[type][i] - first[type][i]) > 2:
                    possibilities[type].append((first[type][i], second[type][i]))

            self.bondPossibilities = possibilities

    def __str__(self):
        """ Prints the protein as a string """

        return self.string

    def prune(self, maxLength, maxStability):
        """ Check if the protein can be pruned, both for 2D and 3D

        :param maxLength: the aminoacid you reached in the protein and after which you might want to prune
        :param maxStability: the max stability found for this protein so far
        :return: True if potential stability > maxStability (attention: negative values!)
        """

        # pruning: initiate
        self.stability(maxLength)
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

