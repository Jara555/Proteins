import random


class ProteinRandomizer():
    """ Creates random Protein"""

    def __init__(self, length, Iterations):
        """ Set and initiate all properties.

        :param length: length of the protein to be created
        :param Iterations: number of random proteins to be generated
        """

        # initialize input variables
        self.length = length
        self.maxIterations = Iterations
        self.name = "ProteinRandomizer"
        self.proteinString = ""

    def run(self):
        """ Runs a maximal amount of iterations in which random proteins be created

        :return: .txt file with the protein and number starting at 100
        """

        # loop over max iterations
        for iteration in range(self.Iterations):

            self.generator()

            # write to text file
            write_file = open(("data/protein" + str(iteration + 100) + ".txt"), "w")
            write_file.write(self.proteinString)

    def generator(self):
        """ Creates random protein string"""

        self.proteinString = ""
        aminoTypes = ["P", "H", "C"]

        for i in range(self.length):
            aminoType = aminoTypes[random.randrange(len(aminoTypes))]
            self.proteinString += aminoType




