class AminoAcid(object):
    """ Contains aminoacid info """

    def __init__(self, type):
        """ Initializes type of the aminoacid """

        self.type = type

    def setCoordinates(self, x, y):
        """ Initializes 2D coordinates of the aminoacid """

        self.x = x
        self.y = y

    def set3Dcoordinates(self, x, y, z):
        """ Initializes 3D coordinates of the aminoacid """

        self.x = x
        self.y = y
        self.z = z



