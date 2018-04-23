
def overlap(protein):
    """ Checks the protein for overlap"""

    coords_list = []

    # put x and y coordinates of aminoacid in x and y lists
    for aminoacid in protein:
        coords = (aminoacid.x, aminoacid.y)
        if coords in coords_list:
            return True
        coords_list.append(coords)

    return False

