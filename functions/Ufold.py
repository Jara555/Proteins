def Ufold(protein):
    """ Creates folding pattern """

    size = len(protein)
    patternFold = []

    if size % 2 == 0:
        fold1 = (size/2) - 1
        fold2 = size/2

    # odd = []
    # even = []
    #
    # for i in range(size):
    #     if protein[i] == 1:
    #         if i % 2 == 0:
    #             even.append(i)
    #         else:
    #             odd.append(i)

    patternFold.append("0")

    for i in range(size - 1):
        if i == fold1:
            patternFold.append("-Y")
        elif i == fold2:
            patternFold.append("-X")
        else:
            patternFold.append("0")

    print(patternFold)

if _name_ == "_main_":
    main()