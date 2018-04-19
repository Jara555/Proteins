def stability(protein):
    """ Determines stability (number of H-bonds) of folded protein """

    # initializes params
    x = []
    y = []
    coordinates = []
    score = 0
    orientation = [0, 1, 0, -1, 1, 0, -1, 0]

    # stores x and y coordinates of aminoacids with type "H"

    # ... as integer
    for i in range(len(protein)):
        if protein[i].type == "H":
            x.append(protein[i].x)
            y.append(protein[i].y)

    # ... as string
    for i in range(len(x)):
        coordinates.append(str(x[i]) + str(y[i]))

    # loops over aminoacids with type "H" and determines number of H-bonds
    for i in range(len(x)):

        for k in range(4):
            ybond = int(coordinates[i][0]) + orientation[(k * 2) + 1]
            xbond = int(coordinates[i][1]) + orientation[k * 2]
            xy = str(ybond) + str(xbond)

            if i == 0:
                if xy in coordinates and coordinates[i + 1] != xy:
                    score = score - 1
            if i == len(x) - 1:
                if xy in coordinates and coordinates[i - 1] != xy:
                    score = score - 1
            else:
                if xy in coordinates and coordinates[i - 1] != xy and coordinates[i + 1] != xy:
                    score = score - 1

    # returns stability
    return score / 2
