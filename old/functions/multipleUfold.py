from old.functions import Ufold

def multipleUfold(protein):
    """ Creates U-fold on multiple locations
        input is: protein-string, output is: folding patterns"""

    # variables
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
    directions = ['+Y', '+X', '-Y', '-X']

    # iterate through protein to find al H's append the location to the oddH or evenH list
    for i in range(size):
        if protein[i] == 'H':
            if i % 2 == 0:                  # i.e. even
                evenH.append(i)
            else:                           # i.e. odd
                oddH.append(i)

    # reality check
    print("oddH " + str(oddH))
    print("evenH" + str(evenH))





    # find the first and second H                       # Todo: not all folding points are found yet
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

    # find folding possibilities within firstH, secondH
    for i in range(len(firstH)):
        if (secondH[i] - firstH[i]) > 2:
            #pattern.append(Ufold(protein[firstH:secondH+1])
            print("found folding possibility: " + str(firstH[i]) + str(secondH[i]))
            possibilities.append((firstH[i], secondH[i]))
    print(possibilities)

    # find combinations of possibilities
    possibilitiesSeconds = [x[1] for x in possibilities]            # Last H of possible fold
    possibilitieFirsts = [x[0] for x in possibilities]              # First H of possible fold

    # make combi with first fold is first of possbilitiesSeconds
    for i in range(len(possibilitiesSeconds)):
        for j in range(len(possibilitieFirsts)):
            if possibilitieFirsts[j] >= possibilitiesSeconds[i]:
                combinations.append((possibilities[i], possibilities[j]))
                print(combinations)

    # check if there is one inbetween others
    for i in range(len(possibilities)):
        for j in range(len(possibilities)):
            if possibilities[j][0] > possibilities[i][0] and possibilities[j][1] < possibilities[i][1]:
                combinations.append((possibilities[i], possibilities[j]))
                print(combinations)

    # make Ufolds
    for i in range(len(combinations)):
        for j in range(len(combinations[i])):
            subProtein = protein[combinations[i][j][0]:combinations[i][j][1]+1]
            print(subProtein)
            pattern[combinations[i][j][0]:combinations[i][j][1]+1] = Ufold(subProtein, directions[j])
            print(pattern)
        length = len(pattern)
        if len(pattern) != len(protein):
            print(pattern[length-1])
            pattern += pattern[length-1]
            print(pattern)
        patterns.append(pattern)
        pattern = []
    print("patterns is: " + str(patterns))


    return patterns










    #
    #     print(pattern)
    #     pattern.append('0')
    #     pattern.append('0')
    #     pattern.append('0')
    #     pattern.append('0')
    #     print(pattern)
    #






