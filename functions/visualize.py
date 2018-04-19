import matplotlib.pyplot as plt


def visualize(protein):
    """ Prints the protein in scatter plot with lines"""

    print()
    print('Protein coordinates:')

    x = []
    y = []
    color = []

    for aminoacid in protein:
        x.append(aminoacid.x)
        y.append(aminoacid.y)
        if aminoacid.type == 'H':
            color.append('red')
        else:
            color.append('blue')

    print(x)
    print(y)
    print()

    plt.plot(x, y, 'C3', zorder=1, lw=3, color='black')
    plt.scatter(x, y, s=120, zorder=2, color=color)
    plt.title('Protein')
    plt.tight_layout()
    #plt.legend()

    plt.show()
