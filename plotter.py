
import argparse
import random
import csv
import bisect
import time
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from collections import Counter

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--x',
                        required=True,
                        help="Name of x")
    parser.add_argument("--y", required=True, help="Name of y")
    parser.add_argument("--filename", required=True, help="help file csv to open")
    parser.add_argument("--label", nargs=1, help="help file csv to open")
    parser.add_argument("--x_name", default="", nargs=1, help="help file csv to open")
    parser.add_argument("--y_name", default="", nargs=1, help="help file csv to open")

    args = parser.parse_args()
    with open(args.filename, 'r') as csvFile:
        reader = csv.DictReader(csvFile)
        all_rows = [row for row in reader]
        x_col = [float(row[args.x]) for row in all_rows]
        y_col = [float(row[args.y]) for row in all_rows]


        plt.hold(True)
        ax = plt.subplot(1,1,1)
        p1, = ax.plot(x_col, y_col, color='g', label=args.label)
        #p2, = ax.plot(r_fs, time_fs, color='b', label='SubSumScan')
        #ax.legend(loc='upper left')
        ax.set_xlabel("Number of Trajectories in Sample")
        ax.set_ylabel("Time(sec)")
        plt.show()
