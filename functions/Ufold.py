def Ufold(protein, direction):
    """ Creates folding pattern """

    size = len(protein)
    patternFold = []

    if size % 2 == 0:
        fold1 = size/2
        fold2 = fold1 + 1
    else:
        fold1 = (size + 1) / 2
        fold2 = fold1 + 1

    # Make a Ufold starting from a certain direction
    for i in range(size):
        if direction == '+Y':
            if i == 0:
                patternFold.append('+Y')
            elif i == fold1:
                patternFold.append('+X')
            elif i == fold2:
                patternFold.append('-Y')
            else:
                direction = patternFold[i - 1]
                patternFold.append(direction)
        elif direction == '-Y':
            if i == 0:
                patternFold.append('-Y')
            elif i == fold1:
                patternFold.append('-X')
            elif i == fold2:
                patternFold.append('+Y')
            else:
                direction = patternFold[i - 1]
                patternFold.append(direction)
        elif direction == '+X':
            if i == 0:
                patternFold.append('+X')
            elif i == fold1:
                patternFold.append('-Y')
            elif i == fold2:
                patternFold.append('-X')
            else:
                direction = patternFold[i - 1]
                patternFold.append(direction)
        elif direction == '-X':
            if i == 0:
                patternFold.append('-X')
            elif i == fold1:
                patternFold.append('+Y')
            elif i == fold2:
                patternFold.append('+X')
            else:
                direction = patternFold[i - 1]
                patternFold.append(direction)

    return patternFold