#!/usr/bin/python
import argparse
import random
import csv
import bisect
import time
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import plotting_tools

from collections import Counter

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--x',
                        required=True,
                        help="Name of x")
    parser.add_argument('--x_diff',
                        help="Name of x")
    parser.add_argument('--y_diff',
                        help="Name of x")
    parser.add_argument("--y", required=True, help="Name of y")
    parser.add_argument("--filenames", required=True, nargs="+", help="help file csv to open")
    parser.add_argument("--labels", nargs="+", help="help file csv to open")
    parser.add_argument("--x_name", default="", help="help file csv to open")
    parser.add_argument("--y_name", default="", help="help file csv to open")
    parser.add_argument("--logx", action="store_true")
    parser.add_argument("--logy", action="store_true")
    parser.add_argument("--smooth", type=float)
    parser.add_argument("--save", help="Name of the saved file")

    args = parser.parse_args()
    ax = plt.subplot(1, 1, 1)

    if args.logx:
        ax.set_xscale("log")
    if args.logy:
        ax.set_yscale("log")
    for filename, label in zip(args.filenames, args.labels):
        with open(filename, 'r') as csvFile:
            reader = csv.DictReader(csvFile)
            all_rows = [row for row in reader]
            if args.x_diff is not None:
                x_col = [abs(float(row[args.x_diff]) - float(row[args.x])) for row in all_rows]
            else:
                x_col = [float(row[args.x]) for row in all_rows]
            if args.y_diff is not None:
                y_col = [abs(float(row[args.y_diff]) - float(row[args.y])) for row in all_rows]
            else:
                y_col = [float(row[args.y]) for row in all_rows]

            #plt.hold(True)
            print(label)

            ax.scatter(x_col, y_col, label=label)
            if args.smooth is not None:
                plotting_tools.plot_interp(ax, x_col, y_col, sig=args.smooth, logsc=args.logx)
            #p2, = ax.plot(r_fs, time_fs, color='b', label='SubSumScan')

    ax.legend()
    ax.set_xlabel(args.x_name)
    ax.set_ylabel(args.y_name)
    plt.tight_layout()
    if args.save is not None:
        plt.savefig(args.save, format="pdf")
    plt.show()
