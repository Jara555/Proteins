import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import copy
from classes.AminoAcid import AminoAcid
from mpl_toolkits.mplot3d import Axes3D


class Protein(object):
    """ Contains all protein properties and methods """

    def __init__(self, dimensions, number, string=None):
        """ Set properties and initialize all aminoacids
        :param number: protein file number
        :param dimensions: 2 for 2D or 3 for 3D
        """

        # initialize input variables
        self.dimensions = dimensions

        # if a number was entered open the protein file of that number
        if number:
            self.number = number
            # open protein text file
            with open('data/protein' + str(number) + '.txt', 'r') as file:
                self.string = file.read()

        # if a string was entered use the string
        elif string:
            self.string = string

        # set class properties
        self.length = len(self.string)
        self.list = []
        self.listH = []
        self.bonds = []
        self.Hbonds = []
        self.Cbonds = []
        self.stabilityScore = 0
        self.bondPossibilities = [[], []]

        # append list with aminoacids
        for aa_index in range(self.length):
            aminoacid = AminoAcid(self.string[aa_index])
            self.list.append(aminoacid)

        # first 2 aminoacids have always same coordinates ['0', 'Y', ... ]
        self.list[0].setCoordinates(0, 0, 0)
        self.list[1].setCoordinates(0, 1, 0)

        # set optimum (only available for first 5 benchmark proteins)
        self.optimum = self.setOptimum()

    def setOptimum(self):
        """ The first 5 proteins of the benchmark proteins in /data have
        according to the literature specific optimal stabilities.
        This method sets the optimum for the current protein. """

        # 2D and 3D optima for proteins 1 - 5 according to literature
        optima2D = [-3, -6, -9, -14, -21]
        optima3D = [-3, -7, -11, -18, -34]

        # only for benchmark proteins 1 - 5
        if self.number in [1, 2, 3, 4, 5]:
            if self.dimensions == 2:
                return optima2D[self.number - 1]
            elif self.dimensions == 3:
                return optima3D[self.number - 1]
        else:
            return None

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

    def checkOverlap(self, maxLength=None):
        """ Checks for overlap in the folded protein, both 2D and 3D
        :param maxLength: until where should the overlap be checked
        :return: True if overlap is found
        """

        if not maxLength:
            maxLength = self.length

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

    def stability(self, maxLength=None):
        """ Determines stability score of protein based on the bonds between "C" and "H" amino types.
        :param: length (i.e. number of aminoacid) of which stability will be determined.
        :return: stability score based on C-C, C-H, and H-H bonds.
        """

        if not maxLength:
            maxLength = self.length

        self.bonds, self.Hbonds, self.Cbonds, aminoType, x, y, z, a = [], [], [], [], [], [], [], []
        currentType, score = 0, 0
        i = -1

        # set orientations
        orientation = self.setOrientationsStability()

        # stores x, y and z coordinates of amino acids with either type "H" or "C"
        x, y, z, a, aminoType = self.storeCoordinatesStability(i, x, y, z, a, aminoType, maxLength)

        # loops over amino acids with type "H" or "C" and determines number of bonds
        for i in range(len(x)):

            # determine index of current "C" or "H" in protein
            counter = 0
            for index in range(self.length):
                if self.string[index] == "C" or self.string[index] == "H":
                    if counter == i:
                        currentType = index
                        break
                    counter = counter + 1

            # iterates over orientations and checks for H- or C-bonds
            for k in range(2*self.dimensions):
                xbond = x[i] + orientation[k][0]
                ybond = y[i] + orientation[k][1]
                zbond = z[i] + orientation[k][2]

                for n in range(len(x)):
                    if i == 0:
                        if (x[n] == xbond and y[n] == ybond and z[n] == zbond) and (
                                self.list[currentType + 1].x != xbond or self.list[currentType + 1].y != ybond or
                                self.list[currentType + 1].z != zbond):
                            # update stability score
                            score = self.updateStability(i, n, a, aminoType, currentType, score)
                    elif i == len(x) - 1:
                        if (x[n] == xbond and y[n] == ybond and z[n] == zbond) and (
                                self.list[currentType - 1].x != xbond or self.list[
                                currentType - 1].y != ybond or self.list[currentType - 1].z != zbond):
                            # update stability score
                            score = self.updateStability(i, n, a, aminoType, currentType, score)
                    else:
                        if (x[n] == xbond and y[n] == ybond and z[n] == zbond) and \
                                (self.list[currentType - 1].x != xbond or
                                 self.list[currentType - 1].y != ybond or
                                 self.list[currentType - 1].z != zbond) and (
                                 self.list[currentType + 1].x != xbond or
                                 self.list[currentType + 1].y != ybond or
                                 self.list[currentType + 1].z != zbond):
                            # update stability score
                            score = self.updateStability(i, n, a, aminoType, currentType, score)

        self.removeDoubleStability()
        self.stabilityScore = score / 2

    def setOrientationsStability(self):
        """" Sets different orientations.
        :return: different orientations
        """

        plusY = [0, 1, 0]
        minY = [0, -1, 0]
        plusX = [1, 0, 0]
        minX = [-1, 0, 0]
        plusZ = [0, 0, 1]
        minZ = [0, 0, -1]

        orientation = [plusX, minX, plusY, minY, plusZ, minZ]

        return orientation

    def storeCoordinatesStability(self, i, x, y, z, a, aminoType, maxLength):
        """ Stores coordinates and type of aminoacids "H" and "C".
        :param: arrays containing x, y and z coordinates and amino type for "C" and "H"
        :return: arrays containing x, y and z coordinates and amino type for "C" and "H"
        """

        for aminoacid in self.list[0:maxLength]:
            i += 1
            if aminoacid.type == "C" or aminoacid.type == "H":
                x.append(aminoacid.x)
                y.append(aminoacid.y)
                z.append(aminoacid.z)
                aminoType.append(aminoacid.type)
                a.append(i)

        return x, y, z, a, aminoType

    def updateStability(self, i, n, a, aminoType, currentType, score):
        """ Updates stability score and bonds.
        :return: stability score
        """

        if (aminoType[i] == "H" and aminoType[n] == "H") or (aminoType[i] == "H" and aminoType[n] == "C") or (
                aminoType[i] == "C" and aminoType[n] == "H"):
            score = score - 1
            self.Hbonds.append((currentType, a[n]))
        else:
            score = score - 5
            self.Cbonds.append((currentType, a[n]))

        return score

    def removeDoubleStability(self):
        """ Removes double bonds.
        :return: array of bond indexes
        """

        self.bonds = [self.Hbonds, self.Cbonds]
        copyBonds = copy.deepcopy(self.bonds)

        for i in range(len(self.bonds)):
            for bond in self.bonds[i]:
                bond2 = bond[1], bond[0]
                if bond2 in copyBonds[i]:
                    self.bonds[i].remove(bond2)
                    copyBonds[i].remove(bond2)
                else:
                    if bond2[1] > bond2[0]:
                        copyBonds[i].append(bond2)
                        copyBonds[i].remove(bond)

        self.bonds = copyBonds

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

        # visualizes CC/CH/HH bonds
        xbond, ybond, zbond = [], [], []
        colorBonds = ["red", "orange"]

        for i in range(len(self.bonds)):
            for j in range(len(self.bonds[i])):
                xbond.append(self.list[self.bonds[i][j][0]].x)
                xbond.append(self.list[self.bonds[i][j][1]].x)
                ybond.append(self.list[self.bonds[i][j][0]].y)
                ybond.append(self.list[self.bonds[i][j][1]].y)
                zbond.append(self.list[self.bonds[i][j][0]].z)
                zbond.append(self.list[self.bonds[i][j][1]].z)
                ax.plot(xbond, ybond, zbond, linestyle=':', color=colorBonds[i])
                xbond, ybond, zbond = [], [], []

        # scatter plot with line
        ax.plot(x, y, z, 'C3', zorder=1, lw=2, color='black')
        ax.scatter(x, y, z, s=100, zorder=2, color=color)

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

        # turn off tick labels
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_zticklabels([])

        # legend
        hydrofoob = mpatches.Patch(color='red', label='H')
        polair = mpatches.Patch(color='blue', label='P')
        cysteine = mpatches.Patch(color='orange', label='C')
        Hbond, = ax.plot(xbond, ybond, zbond, lw=1, color="red", linestyle=':', label="stability -1")
        Cbond, = ax.plot(xbond, ybond, zbond, lw=1, color="orange", linestyle=':', label="stability - 5")
        plt.legend(handles=[hydrofoob, polair, cysteine, Hbond, Cbond])

        plt.show()

    def findHs(self):
        """ Find the indexes of H's in the protein, both for 2D and 3D"""

        # remember H-indices
        for i in range(self.length):
            if self.list[i].type == "H":
                self.listH.append(i)

    def findBonds(self):
        """ Checks which H/C's can make a HH-bond and CC-bond, both for 2D and 3D
        :return: self.bondPossibilities a list with 2 lists of tuples of HBond-indices en CBond-indices
        """

        # empty lists combining for two aminotypes
        odd = []
        even = []
        first = []
        second = []
        i = -1

        # iterate through protein to find all H/C's and append the location to the odd or even list
        for aminoacid in self.list:
            i += 1
            if aminoacid.type == "C" or aminoacid.type == "H":
                if i % 2 == 0:  # i.e. even
                    even.append(i)
                else:  # i.e. odd
                    odd.append(i)

        # rule: the second H/C of a bond has to be at a later position in the protein than the first H/C
        # make a list of first and second H/C's
        for i in range(len(even)):
            for j in range(len(odd)):
                if even[i] < odd[j]:
                    first.append(even[i])
                    second.append(odd[j])
                else:
                    first.append(odd[j])
                    second.append(even[i])

        # rule: the distance between two H/C's must be at least 2 to be able to fold
        # make a list of tuples representing the first and second H/C of a possible folding
        for i in range(len(first)):
            if (second[i] - first[i]) > 2:
                self.detectBondKind(first[i], second[i])

    def detectBondKind(self, first, second):
        """ Detect which kind of bond is found

        :param first: the index of the first aminoacid of the bond
        :param second: the index of the second aminoacid of the bond
        """

        # if both first and second aminoacid are C: CC-bond, otherwise H-bond
        if self.list[first].type == "C" and self.list[second].type == "C":
            self.bondPossibilities[1].append((first, second))
        else:
            self.bondPossibilities[0].append((first, second))

    def __str__(self):
        """ Prints the protein as a string """

        return self.string

    def prune(self, maxLength, maxStability):
        """ Check if the protein can be pruned, both for 2D and 3D
        :param maxLength: the aminoacid you reached in the protein and after which you might want to prune
        :param maxStability: the max stability found for this protein so far
        :return: True if potential stability worse than maxStability
        """

        # pruning: initiate
        self.stability(maxLength)
        bondOptions = copy.deepcopy(self.bondPossibilities)
        newBondOptions = [[], []]
        stabilityEffect = [-1, -5]

        # count down from 1 to 0, to start with C-bonds then H-bonds
        for type in range(1, -1, -1):
            counter = [0] * self.length

            # remove the bonds that are in the protein already, count bonds per aminoacid in counter
            if self.bonds[type] != []:
                for bond in self.bonds[type]:
                    bondOptions[type].remove(bond)
                    counter[bond[0]] += 1
                    counter[bond[1]] += 1

            # make a new list of only bonds with second aminoacid after k (still potential bonds)
            for bond in bondOptions[type]:
                if bond[1] >= maxLength:
                    newBondOptions[type].append(bond)

            # correct newbondOptions for max bindings per aminoacid
            newBondOptions[type] = self.updateBondOptions(counter, newBondOptions[type])

        # calculate potential stability
        potStability = self.stabilityScore + (stabilityEffect[0] * len(newBondOptions[0]) +
                (stabilityEffect[1] * len(newBondOptions[1])))

        # prune if potential Stability is worse than maxStability
        if potStability >= maxStability:
            return True

        return False

    def updateBondOptions(self, counter, bondOptions):
        """ Correct new bondOptions for max bindings per aminoacid.

        :param counter: counter for bonds per aminoacid
        :param bondOptions: list of possibilities to make bonds (index of aminoacids)
        :return: list of possible bonds, corrected for max bonds per
        """

        bondOptionsMax = []

        if bondOptions == []:
            return bondOptions

        # set max for inner and outer aminoacids depending on dimensions
        if self.dimensions == 2:
            outerMax = 3
            innerMax = 2
        else:
            outerMax = 5
            innerMax = 4

        # check if aminoacid doesn't have to many bonds
        for bond in reversed(bondOptions):
            if bond[0] == 0 or bond[1] == self.length:
                if counter[bond[0]] < outerMax and counter[bond[1]] < outerMax:
                    bondOptionsMax.append(bond)
                    counter[bond[0]] += 1
                    counter[bond[1]] += 1
            else:
                if counter[bond[0]] < innerMax and counter[bond[1]] < innerMax:
                    bondOptionsMax.append(bond)
                    counter[bond[0]] += 1
                    counter[bond[1]] += 1

        return bondOptionsMax