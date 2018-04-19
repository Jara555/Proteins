def stability(protein):

    x = []
    y = []
    score = 0
    orientation = [0, 1, 0, -1, 1, 0, -1, 0]

    for i in range(len(protein)):
        if protein[i].type == "H":
            x.append(protein[i].x)
            y.append(protein[i].y)

    coordinates = []

    for i in range(len(x)):
        coordinates.append(str(x[i]) + str(y[i]))

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

    return score / 2
