import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

from classes.AminoAcid import AminoAcid


class Protein(object):
    """ Contains all protein properties and methods """

    def __init__(self, number):
        """ Set properties and initializes all aminoacids """

        self.number = number

        # open protein text file
        with open('data/protein' + str(number) + '.txt', 'r') as file:
            self.string = file.read()

        # set properties
        self.length = len(self.string)
        self.list = []

        # append list with aminoacids
        for aa_index in range(self.length):
            aminoacid = AminoAcid(self.string[aa_index])
            self.list.append(aminoacid)

    def fold(self, folding_pattern):
        """ Folds protein according to input pattern """

        # let 1st aminoacid start at coordinates (0, 0)
        x, y = 0, 0
        self.list[0].setCoordinates(x, y)

        # set all possible orientations
        plusY = [0, 1]
        minY = [0, -1]
        plusX = [1, 0]
        minX = [-1, 0]

        # set default orientation
        orientation = plusY

        # iterate over aminoacids in protein (skipping the 1st)
        for index in range(1, self.length):
            # set orientation to folding pattern
            if folding_pattern[index] == '+X':
                orientation = plusX
            elif folding_pattern[index] == '-X':
                orientation = minX
            elif folding_pattern[index] == '+Y':
                orientation = plusY
            elif folding_pattern[index] == '-Y':
                orientation = minY

            # set new coordinates based on orientation
            x = x + orientation[0]
            y = y + orientation[1]

            # set new coordinates to aminoacid
            self.list[index].setCoordinates(x, y)

    def checkOverlap(self):
        """ Checks the protein for overlap"""

        coordinatesList = []

        # iterate over aminoacids
        for aminoacid in self.list:
            coordinates = (aminoacid.x, aminoacid.y)

            # if coordinates already exist: overlap detected
            if coordinates in coordinatesList:
                self.overlap = True

            # else put in the list
            coordinatesList.append(coordinates)

        # if iteration finished there was no overlap
        self.overlap = False

    def stability(self):
        """ Checks the stability of the protein """

        x = []
        y = []
        score = 0
        orientation = [1, 0, -1, 0, 0, 1, 0, -1]
        currentH = 0

        # stores x and y coordinates of aminoacids with type "H"
        for i in range(self.length):
            if self.list[i].type == "H":
                x.append(self.list[i].x)
                y.append(self.list[i].y)

        # loops over aminoacids with type "H" and determines number of H-bonds
        for i in range(len(x)):
            # determine index of H in protein
            counter = 0

            for index in range(self.length):
                if self.string[index] == "H":
                    if counter == i:
                        currentH = index
                        break
                    counter = counter + 1

            # iterates over orientations and checks for H-bonds
            for k in range(4):
                xbond = x[i] + orientation[(k * 2)]
                ybond = y[i] + orientation[(k * 2) + 1]

                for n in range(len(x)):
                    if i == 0:
                        if (x[n] == xbond and y[n] == ybond) and (
                           self.list[currentH + 1].x != xbond or self.list[currentH + 1].y != ybond):
                            score = score - 1
                    elif i == len(x) - 1:
                        if (x[n] == xbond and y[n] == ybond) and (
                           self.list[currentH - 1].x != xbond or self.list[currentH - 1].y != ybond):
                            score = score - 1
                    else:
                        if (x[n] == xbond and y[n] == ybond) and \
                                (self.list[currentH - 1].x != xbond or self.list[currentH - 1].y != ybond) and \
                                (self.list[currentH + 1].x != xbond or self.list[
                                    currentH + 1].y != ybond):
                            score = score - 1

    def visualize(self, name):
        """ Prints the protein in scatter plot with lines"""

        x = []
        y = []
        color = []

        # put x and y coordinates of aminoacid in x and y lists
        for aminoacid in self.list:
            x.append(aminoacid.x)
            y.append(aminoacid.y)

            # creates color list for H = red and P = blue
            if aminoacid.type == 'H':
                color.append('red')
            else:
                color.append('blue')

        # print coordinates in terminal
        print()
        print('Protein coordinates:')
        print(x)
        print(y)
        print()

        # scatter plot with line
        plt.plot(x, y, 'C3', zorder=1, lw=2, color='black')
        plt.scatter(x, y, s=50, zorder=2, color=color)

        # layout
        plt.title('Protein: ' + name), plt.tight_layout(), plt.axis('scaled')
        plt.ylim(min(y) - 1, max(y) + 1), plt.xlim(min(x) - 1, max(x) + 1)
        plt.xticks(np.arange(min(x), max(x) + 1, 1.0)), plt.yticks(np.arange(min(y), max(y) + 1, 1.0))

        # legend
        hydrofoob = mpatches.Patch(color='red', label='H')
        polair = mpatches.Patch(color='blue', label='P')
        plt.legend(handles=[hydrofoob, polair])

        plt.show()

    def findHbonds(self):
        """ Checks which H's can make a H-bond """

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

        self.combinations = combinations

    def __str__(self):
        """ Prints the protein as a string """

        return self.string
