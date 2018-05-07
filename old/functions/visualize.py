import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


def visualize(protein, name):
    """ Prints the protein in scatter plot with lines"""

    size = len(protein)

    print()
    print('Protein coordinates:')

    x = []
    y = []
    color = []

    # put x and y coordinates of aminoacid in x and y lists
    for aminoacid in protein:
        x.append(aminoacid.x)
        y.append(aminoacid.y)

        # creates color list for H = red and P = blue
        if aminoacid.type == 'H':
            color.append('red')
        else:
            color.append('blue')

    # print coordinates in terminal
    print(x)
    print(y)
    print()

    # scatter plot with line
    plt.plot(x, y, 'C3', zorder=1, lw=2, color='black')
    plt.scatter(x, y, s=50, zorder=2, color=color)
    plt.title('Protein: ' + name)
    plt.tight_layout()
    plt.axis('scaled')
    #
    plt.ylim(min(y) - 1, max(y) + 1)
    plt.xlim(min(x) - 1, max(x) + 1)

    plt.xticks(np.arange(min(x), max(x) + 1, 1.0))
    plt.yticks(np.arange(min(y), max(y) + 1, 1.0))

    # plt.axhline(0, linestyle='--', color='gray', linewidth=0.5)  # horizontal lines
    # plt.axvline(0, linestyle='--', color='gray', linewidth=0.5)  # vertical lines

    hydrofoob = mpatches.Patch(color='red', label='H')
    polair = mpatches.Patch(color='blue', label='P')
    plt.legend(handles=[hydrofoob, polair])

    plt.show()
