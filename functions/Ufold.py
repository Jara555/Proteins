def Ufold(protein, noFold):
    """ Creates folding pattern """

    size = len(protein)
    patternFold = []

    if size % 2 == 0:
        fold1 = (size/2) - 1
        fold2 = size/2
    else:
        fold2 = (size + 1) / 2
        fold1 = fold2 - 1

    # odd = []
    # even = []
    #
    # for i in range(size):
    #     if protein[i] == 1:
    #         if i % 2 == 0:
    #             even.append(i)
    #         else:
    #             odd.append(i)

    if noFold == 0:
        patternFold.append('0')

        for i in range(size - 1):
            if i == fold1:
                patternFold.append('+X')
            elif i == fold2:
                patternFold.append('-Y')
            else:
                direction = patternFold[i - 1]
                patternFold.append(direction)
    elif noFold == 1:
        patternFold.append('0')

        for i in range(size - 1):
            if i == fold1:
                patternFold.append('-X')
            elif i == fold2:
                patternFold.append('+Y')
            else:
                direction = patternFold[i - 1]
                patternFold.append(direction)


    return patternFold