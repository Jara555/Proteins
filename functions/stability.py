def stability(protein):
    """ Determines stability (number of H-bonds) of folded protein """

    x = []
    y = []
    score = 0
    orientation = [1, 0, -1, 0, 0, 1, 0, -1]
    currentH = 0

    # stores x and y coordinates of aminoacids with type "H"
    for i in range(len(protein)):
        if protein[i].type == "H":
            x.append(protein[i].x)
            y.append(protein[i].y)

    print("x")
    print(x)
    print("y")
    print(y)

    # loops over aminoacids with type "H" and determines number of H-bonds
    for i in range(len(x)):

        # determine index of H in protein
        counter = 0

        for index in range(len(protein)):
            if protein[index].type == "H":
                if counter == i:
                    currentH = index
                    break
                counter = counter + 1

        # iterates over orientations and checks for H-bonds
        for k in range(4):

            xbond = x[i] + orientation[(k * 2)]
            ybond = y[i] + orientation[(k * 2) + 1]

            for n in range(len(x)):
                if i == 0:
                    if (x[n] == xbond and y[n] == ybond) and (protein[currentH + 1].x != xbond or protein[currentH + 1].y != ybond):
                        score = score - 1
                elif i == len(x) - 1:
                    if (x[n] == xbond and y[n] == ybond) and (protein[currentH - 1].x != xbond or protein[currentH - 1].y != ybond):
                        score = score - 1
                else:
                    if (x[n] == xbond and y[n] == ybond) and \
                            (protein[currentH - 1].x != xbond or protein[currentH - 1].y != ybond) and \
                            (protein[currentH + 1].x != xbond or protein[currentH + 1].y != ybond):
                        score = score - 1

    # returns stability
    return score / 2
