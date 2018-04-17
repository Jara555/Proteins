
def foldprotein(protein, folding_pattern):
    """ Folds protein according to input pattern """

    # let 1st aminoacid start at coordinates (0, 0)
    x_index = 0
    y_index = 0
    protein[0].x = x_index
    protein[0].y = y_index

    # declare orientations
    right = [0, 1]
    left = [0, -1]
    down = [1, 0]
    up = [-1, 0]

    # set default orientation
    orientation = right

    # iterate over aminoacids in protein (skipping the 1st)
    for index in range(1, len(protein)):

        # adapting orientation to folding pattern
        if folding_pattern[index] == '+X':
            orientation = down
        elif folding_pattern[index] == '-X':
            orientation = up
        elif folding_pattern[index] == '+Y':
            orientation = right
        elif folding_pattern[index] == '-Y':
            orientation = left

        # set new index based on (new) orientation
        x_index = x_index + orientation[0]
        y_index = y_index + orientation[1]

        # set new coordinates of aminoacid
        protein[index].x = x_index
        protein[index].y = y_index

    # return rearranged protein
    return protein

