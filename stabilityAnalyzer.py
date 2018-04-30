import csv
import numpy as np
import matplotlib.pyplot as plt


def stabilityAnalyzer(csvname):
    low = 0
    score = []
    filename = "results/" + csvname + ".csv"

    with open(filename, 'r') as csvfile:
        # determine lowest stability and adds stability to list
        next(csvfile)
        for row in csv.reader(csvfile):
            if row:
                stability = int(abs(float(row[0])))
                score.append(stability)
                if stability > low:
                    low = stability

    # stores occurrences of stability scores
    low = low + 1
    y = np.zeros(low)
    for item in score:
        y[item] = y[item] + 1

    print(y)

    x = []
    start = 0

    # create array for x axis
    for i in range(low):
        x.append(start)
        start = start + 1

    # creates plot
    plt.bar(x, y)
    plt.xticks(np.arange(min(x), max(x) + 1, 1.0))
    plt.title("Stability scores")
    plt.xlabel("Stability score")
    plt.ylabel("Count")

    # add labels
    label = []
    for item in range(low):
        label.append(str(int(abs(float(y[item])))))

    count = 0
    for bar in range(low):
        plt.text(x=count, y=y[count] + 5, s=label[count])
        count = count + 1

    plt.show()
