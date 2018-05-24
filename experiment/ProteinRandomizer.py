import csv
import random


class ProteinRandomizer():
    """ Creates random Protein"""

    def __init__(self, length, iterations, fixedHNumber):
        """ Set and initiate all properties.

        :param length: length of the protein to be created
        :param iterations: number of random proteins to be generated
        """

        # initialize input variables
        self.length = length
        self.iterations = iterations
        self.fixedHNumber = fixedHNumber
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

        write_results = ("results/experimentProteins" + ".csv")
        with open(write_results, 'w', newline='') as resultsfile:
            self.writer = csv.writer(resultsfile)

            if self.fixedHNumber is None:
                # loop over max iterations
                for iteration in range(self.iterations):

                    self.generator()

                    # write to text file
                    write_file = open(("data/protein" + str(iteration + 10000) + ".txt"), "w")
                    write_file.write(self.proteinString)

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

    def generatorFixedH(self, n):

        # create random protein with n H's


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



