import pyscan
import matplotlib.pyplot as plt
import random
import math

if __name__ == "__main__":
    st_pt = (0,0)
    pts = [st_pt]
    for i in range(100):
        pts.append((random.random() / 100.0 + pts[-1][0], random.random() / 100.0 + pts[-1][1]))

    alpha = .01
    chord_l = 2 * math.sqrt(2 * alpha * .1- alpha * alpha)
    pts = pyscan.approx_traj_cells(pts, chord_l, alpha)
    print(pts)
