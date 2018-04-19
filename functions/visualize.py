import matplotlib.pyplot as plt


def visualize(protein):
    """ Prints the protein in scatter plot with lines"""

    print()
    print('Protein coordinates:')

    x = []
    y = []

    for aminoacid in protein:
        x.append(aminoacid.x)
        y.append(aminoacid.y)

    print(x)
    print(y)
    print()

    plt.plot(x, y, 'C3', zorder=1, lw=3)
    plt.scatter(x, y, s=120, zorder=2)
    plt.title('Protein')
    plt.tight_layout()

    plt.show()
