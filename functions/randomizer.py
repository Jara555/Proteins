from random import randint


def randomizer(protein):
    """ Creates random folding patterns """

    # get size protein and create empty folding pattern list
    size = len(protein)
    patternFold = []

    # index counter for number of items in folding pattern list
    i = 0

    # iterate over required folding pattern length
    while i <= size:

        # get random orientation index
        orientation = randint(0, 4)

        # first element should always start at 0, 0
        if i == 0:
            patternFold.append('0')
            i += 1
        elif orientation == 0:
            patternFold.append('0')
            i += 1
        # previous element cannot be inverse orientation
        elif orientation == 1 and patternFold[i - 1] != '-X':
            patternFold.append('+X')
            i += 1
        elif orientation == 2 and patternFold[i - 1] != '+X':
            patternFold.append('-X')
            i += 1
        elif orientation == 3 and patternFold[i - 1] != '-Y':
            patternFold.append('+Y')
            i += 1
        elif orientation == 4 and patternFold[i - 1] != '+Y':
            patternFold.append('-Y')
            i += 1

    return patternFold