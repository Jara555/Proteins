
class AminoAcid(object):

    def __init__(self, type, x, y):

        self.type = type
        self.x = x
        self.y = y

        self.prev = y - 1
        self.nex = y + 1

    def __str__(self):
        return (self.type + " -")

def main():

    h4 = AminoAcid('H', 0, 3)
    h1 = AminoAcid('H', 0, 0)
    h2 = AminoAcid('H', 0, 1)
    p3 = AminoAcid('P', 0, 2)
    h4 = AminoAcid('H', 0, 3)
    h1 = AminoAcid('H', 0, 0)
    h2 = AminoAcid('H', 0, 1)
    p3 = AminoAcid('P', 0, 2)

    print()
    print(h1, h2, p3, h4)
    print()


if __name__ == "__main__":
    main()