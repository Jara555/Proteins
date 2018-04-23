from functions.Ufold import Ufold

def multipleUfold(protein):
    """ Creates U-fold on multiple locations
        input is: protein-string, output is: folding patterns"""

    # Variables
    size = len(protein)
    patternFold = []
    oddH = []
    evenH = []
    firstH = []
    secondH = []
    pattern = []
    possibilities = []
    combinations = []
    patterns = []

    # Iterate through protein to find al H's append the location to the oddH or evenH list
    for i in range(size):
        if protein[i] == 'H':
            if i % 2 == 0:                  # i.e. even
                evenH.append(i)
            else:                           # i.e. odd
                oddH.append(i)

    # reality check
    print("oddH " + str(oddH))
    print("evenH" + str(evenH))





    # Find the first and second H                       # Todo: not all folding points are found yet
    for i in range(len(evenH)):
        for j in range(len(oddH)):
            if evenH[i] < oddH[j]:
                firstH.append(evenH[i])
                secondH.append(oddH[j])
            else:
                firstH.append(oddH[j])
                secondH.append(evenH[i])

    print("firstH:" + str(firstH))
    print("secondH:" + str(secondH))

    # Find folding possibilities within firstH, secondH
    for i in range(len(firstH)):
        if (secondH[i] - firstH[i]) > 2:
            #pattern.append(Ufold(protein[firstH:secondH+1])
            print("found folding possibility: " + str(firstH[i]) + str(secondH[i]))
            possibilities.append((firstH[i], secondH[i]))


    # Find combinations within the folding possibilities
    # TODO: now only in combination with the fist option. Need another looping strategy
    # TODO: now only combination of two possibilities?
    print(possibilities)
    possibilitiesSeconds = [x[1] for x in possibilities]
    possibilitieFirsts = [x[0] for x in possibilities]
    for i in range(len(possibilitieFirsts)):
        for j in range(len(possibilitieFirsts)):
            if (possibilitieFirsts[i] >= possibilitiesSeconds[j]):
                combinations.append((possibilities[j], possibilities[i]))       # Tuple of tuples?

    # make Ufold on right points in different combinations
    print("combinations is:" + str(combinations))
    print(len(combinations))
    # for i in range(len(combinations)):
    #     for j in range(len(combinations)):
    #         #pattern += Ufold(protein[combinations[i][j][0]:combinations[i][j][1]+1], j)
    #         print("combinations[i][j] is: " + str(combinations[i][j]) + " i = " + str(i) + " j = " + str(j))
    #         print("combinations[i][j][0] is: " + str(combinations[i][j][0]))
    #     print(pattern)
    #     patterns.append(pattern)
    #     print("patterns is: " + str(patterns))
    #     pattern = []
    # print("patterns is: " + str(patterns))
    #
    # return patterns










    #
    #     print(pattern)
    #     pattern.append('0')
    #     pattern.append('0')
    #     pattern.append('0')
    #     pattern.append('0')
    #     print(pattern)
    #






