import random
import bisect
import math
import itertools
import csv
import time
import numpy as np

import numpy.random as npr

import pyscan
from collections import deque
from geometric import approx_eq, Segment, Line, to_line


def horizontal_split_vertices(points, segment):
    up = []
    down = []
    for p in points:
        if segment.pt_eq_below(p):
            down.append(p)
        else:
            up.append(p)
    return up, down


def score(end_point, w_lines):
    b_score = 0
    r_score = 0
    for line, rw, bw in w_lines:
        if line.pt_eq_below(end_point):
            r_score += rw
            b_score += bw
    return b_score, r_score




def simple_generate_blue_and_red(full_points, q, r, eps=.01, labels=[], region="line"):

    red_points = []
    blue_points = []
    red_labels = []
    blue_labels = []


    label_set = set(labels)
    red_label_set = set()

    for l in label_set:
        if .5 < random.random():
            red_label_set.add(l)


    for p, label in zip(full_points, labels):
        if label in red_label_set:
            red_points.append(p)
            red_labels.append(label)
        else:
            blue_points.append(p)
            blue_labels.append(label)

    return red_points, red_labels, blue_points, blue_labels


def generate_blue_and_red(full_points, q, r, eps=.01, labels=[], region="line"):
    #724 931 2294

    if region == "line":
        l = line_planted_test(int(2/eps + 1), full_points, r, labels)
    elif region == "disk":
        l = disk_planted_test(int(2/eps + 1), full_points, r, labels)
    else:
        return
    red_points = []
    blue_points = []
    red_labels = []
    blue_labels = []
    red_below = 0
    blue_below = 0
    labels_below = set()
    labels_above = set()

    if not labels:
        labels = range(len(full_points))

    if region == "line":
        for p, label in zip(full_points, labels):

            if l.above_closed(p):
                labels_below.add(label)
            else:
                labels_above.add(label)
    elif region == "disk":
        for p, label in zip(full_points, labels):
            if l.contains(p):
                labels_below.add(label)
            else:
                labels_above.add(label)

    labels_part_below = labels_below
    labels_com_above = labels_above - labels_below
    red_label_set = set()

    for l in labels_part_below:
            if .5 + q < random.random():
                red_below += 1
                red_label_set.add(l)
            else:
                blue_below += 1
    for l in labels_com_above:
            if .5 < random.random():
                red_label_set.add(l)

    for p, label in zip(full_points, labels):
        if label in red_label_set:
            red_points.append(p)
            red_labels.append(label)
        else:
            blue_points.append(p)
            blue_labels.append(label)

    print(abs(red_below / len(red_points) - blue_below / len(blue_points)))
    return red_points, red_labels, blue_points, blue_labels


if __name__ == "__main__":
    print("running experiment")
    import partitioning_tests
    #pts = partitioning_tests.upload_crimes("crimes.csv")
    pts = [pyscan.Point(random.random(), random.random(), 1.0) for _ in range(1000)]

    red_points, blue_points = generate_blue_and_red(pts, .2, .4)
    # red_points = [(random.random(), random.random()) for _ in range(1000000)]
    # blue_points = [(random.random(), random.random()) for _ in range(1000000)]

    # red_points = random.sample(red_points, 100000)
    # blue_points = random.sample(blue_points, 100000)

    # import matplotlib.pyplot as plt
    # f, ax = plt.subplots()
    # x, y = zip(*red_points)
    # ax.scatter(x, y, color="red")
    # x, y = zip(*blue_points)
    # ax.scatter(x, y, color="blue")
    #
    # plt.show()

    # for s_alg in ["naive", "quad", "mat", "chan", "chans"]:
    #     for sc_alg in ["fast"]:
    #         print(s_alg + " " + sc_alg)
    #         testing_framework(red_points, blue_points, -.5, -2.5, 80, input_size=min(len(red_points), len(blue_points)), sampling_alg=s_alg, scan_alg=sc_alg)
