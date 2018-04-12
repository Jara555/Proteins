
class aminoAcid(object):

    def __init__(self, type, x, y):

        self.type = type
        self.x = x
        self.y = y

        self.prev = y - 1
        self.nex = y + 1

    def __str__(self):
        return "type = yes gllafd"

def main():

    H1 = aminoAcid('H', 0, 0)

    print(H1)







if __name__ == "__main__":
    main()