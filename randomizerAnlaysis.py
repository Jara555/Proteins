import csv
import numpy as np
import matplotlib.pyplot as plt


def main():

    low = 0
    score = []

    with open('results/run3.csv', 'r') as csvfile:
        # determine lowest stability and adds stability to list
        next(csvfile)
        for row in csv.reader(csvfile):
            if row:
                score.append(int(abs(float(row[1]))))
                if int(abs(float(row[1]))) > low:
                    low = int(abs(float(row[1])))

    # stores occurrences of stability scores
    low = low + 1
    y = np.zeros(low)
    for item in score:
        y[item] = y[item] + 1

    x = []
    start = 0

    # create array for x axis
    for i in range(low):
        x.append(start)
        start = start - 1

    # creates plot
    plt.bar(x, y)
    plt.xticks(np.arange(min(x), max(x) + 1, 1.0))
    plt.title("Stability scores")
    plt.xlabel("Stability score")
    plt.ylabel("Count")
    plt.show()


if __name__ == "__main__":
    main()
