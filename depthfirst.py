import csv


def main():
    """ Implements depth first algorithms in order to most efficiently fold a protein """

    length = 8

    foldingPattern = [None]*length
    foldingPattern[0] = '0'
    foldingPattern[1] = '+Y'

    k = 3
    options = ['+Y', '-X', '+X', '-Y']

    print(foldingPattern)

    write_file = ('results/run' + 'TEST' + '.csv')
    with open(write_file, 'w') as csvfile:
        fieldnames = ['foldingPattern']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        searching(foldingPattern, length, k, options, writer)


def searching(foldingPattern, length, k, options, writer):

        for option in options:
            if k == length:
                foldingPattern[k - 1] = option
                writer.writerow({'foldingPattern': foldingPattern})

            else:
                foldingPattern[k - 1] = option
                searching(foldingPattern, length, k + 1, options, writer)



if __name__ == "__main__":
    main()
