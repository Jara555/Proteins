# the AminoAcid class sets up the properties of each aminoacid
class AminoAcid(object):

    # uses a char (H/P) and coordinates as input
    def __init__(self, type, x, y):

        # set type and type specific code
        self.type = type
        if self.type == 'H':
            self.code = -1
        elif self.type == 'P':
            self.code = 1
        elif self.type == 'C':
            self.code = -5

        # set coordinates
        self.x = x
        self.y = y

        # link to previous and next aminoacid
        if self.y != 0:
            self.prev = y - 1
        else:
            self.prev = None

        self.nex = y + 1
