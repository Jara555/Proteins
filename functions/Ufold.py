def Ufold(protein, noFold):
    """ Creates folding pattern """

    size = len(protein)
    patternFold = []

    if size % 2 == 0:
        fold1 = size/2
        fold2 = fold1 + 1
    else:
        fold1 = (size + 1) / 2
        fold2 = fold1 + 1

    print(fold1)
    print(fold2)

    if noFold == 0:

        for i in range(size):
            # starting direction is +Y
            if i == 0:
                patternFold.append('+Y')
            elif i == fold1:
                patternFold.append('+X')
            elif i == fold2:
                patternFold.append('-Y')
            else:
                direction = patternFold[i - 1]
                patternFold.append(direction)

    elif noFold == 1:
        patternFold.append('0')

        for i in range(size - 1):
            # starting direction
            if i == 0:
                patternFold.append('+Y')
            elif i == fold1:
                patternFold.append('-X')
            elif i == fold2:
                patternFold.append('+Y')
            else:
                direction = patternFold[i - 1]
                patternFold.append(direction)

    return patternFold