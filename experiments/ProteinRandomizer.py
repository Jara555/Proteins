import csv
import random


class ProteinRandomizer():
    """ Creates random Protein"""

    def __init__(self, length, iterations, fixedHNumber, dimensions):
        """ Set and initiate all properties.

        :param length: length of the protein to be created
        :param iterations: number of random proteins to be generated
        """

        # initialize input variables
        self.length = length
        self.iterations = iterations
        self.fixedHNumber = fixedHNumber
        self.dimensions = dimensions
        self.name = "ProteinRandomizer"
        self.proteinString = ""
        self.countH = 0
        self.countC = 0
        self.clusters = 0
        self.maxClusterLength = 0
        self.clusterLength = 0



    def run(self):
        """ Runs a maximal amount of iterations in which random proteins be created

        :return: .txt file with the protein and number starting at 100
        """

        minStability = ''
        Stability = ''

        writeResults = ("results/experimentProteins-l" + str(self.length) + "-d" + str(self.dimensions) + "-f" +
                        str(self.fixedHNumber) + ".csv")
        with open(writeResults, 'w', newline='') as resultsfile:
            self.writer = csv.writer(resultsfile)


            # loop over max iterations
            for iteration in range(self.iterations):

                if self.fixedHNumber is None:
                    self.generator()
                elif self.fixedHNumber == 0:
                    self.generatorFixedHall()
                    break
                else:
                    self.generatorFixedH(self.fixedHNumber)

                # write to text file
                writeFile = open(("data/protein" + str(iteration + 10000) + ".txt"), "w")
                writeFile.write(self.proteinString)

                # write to csv file
                self.writer.writerow([self.countH, self.maxClusterLength, self.clusters, minStability, Stability])


    def generator(self):
        """ Creates random protein string"""

        # set all properties to null
        self.proteinString = ""
        self.countH = 0
        self.countC = 0
        self.clusters = 0
        self.maxClusterLength = 0

        aminoTypes = ["P", "H"]

        # create random protein
        for i in range(self.length):
            aminoType = aminoTypes[random.randrange(len(aminoTypes))]
            self.proteinString += aminoType

            self.saveProperties(aminoType, i)

    def generatorFixedH(self, fixedHNumber):
        """ creates random protein string with fixed number of H's

        :param fixedHNumber: the number of H in the protein to be generated
        """

        # set all properties to null
        self.proteinString = ""
        self.countH = 0
        self.countC = 0
        self.clusters = 0
        self.maxClusterLength = 0

        population = []

        # make a population list for random sampling without replacement
        for i in range(self.length):
            population.append(i)

        # create a list of random H locations
        HList = random.sample(population, fixedHNumber)

        # create protein
        for i in range(self.length):
            if i in HList:
                aminoType = "H"
                self.proteinString += aminoType
            else:
                aminoType = "P"
                self.proteinString += aminoType

            self.saveProperties(aminoType, i)


    def saveProperties(self, aminoType, k):
        """ save properties of the protein

        :param aminoType: the type of the aminoacid just added
        :param k: the index of the aminoacid just added
        """

        # count cluster properties if the added aminoacid is a H
        if aminoType == "H":
            self.countH += 1
            self.clusterLength += 1

            if k == self.length-1:
                if self.clusterLength > self.maxClusterLength:
                    self.maxClusterLength = self.clusterLength
                    self.clusters += 1

        # save cluster and start with new cluster if aminoacid is a P
        elif aminoType == "P":

            # save if the longest cluster of H's
            if self.clusterLength > self.maxClusterLength:
                self.maxClusterLength = self.clusterLength

            if self.clusterLength > 0:
                self.clusters += 1

            self.clusterLength = 0



