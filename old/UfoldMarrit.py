# def Ufold(protein):


def Ufold(proteinList):

    # Determine length
    length = len(proteinList)

    # Fold at half of the length
    if length % 2 == 0:         # even
        firstFold = length/2
    else:                       #odd
        length-1
        firtsFold = length/2

    secondFold = firstFold + 1

    # fold at first fold to right
    print(proteinList[int(firstFold)].x)
    proteinList[int(firstFold)].x += 1
    print(proteinList[int(firstFold)].x)

    # fold at second fold to right again
    proteinList[int(secondFold)].y -= 1

    return proteinList










    # Determine how many even H

    # Determine how many odds H

    # Last of evens

    # First of odds
