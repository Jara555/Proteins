class AminoAcid(object):
    """ Contains aminoacid info """

    def __init__(self, type):
        """ Initializes type of the aminoacid """

        self.type = type
        self.x = 0
        self.y = 0
        self.z = 0

    def setCoordinates(self, x, y, z=0):
        """ Initializes coordinates of the aminoacid """

        self.x = x
        self.y = y
        self.z = z




