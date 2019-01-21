import pyscan
import matplotlib.pyplot as plt
import random
import math

import statistics

def plot_points(ax, pts, c):
    xs = []
    ys = []
    for pt in pts:
        xs.append(pt[0] )
        ys.append(pt[1])
    ax.scatter(xs, ys, color=c)

def plot_points_traj(ax, pts, c):
    xs = []
    ys = []
    for pt in pts:
        xs.append(pt.getX())
        ys.append(pt.getY())
    ax.plot(xs, ys, color=c)

def plot_tuples(ax, pts, c):
    xs = []
    ys = []
    for pt in pts:
        xs.append(pt[0])
        ys.append(pt[1])
    ax.plot(xs, ys, color=c)
    #ax.scatter(xs, ys, color=c)


def boxed_trajectory(count):
    """
    Creates a trajectory constrained to a box
    :return:
    """
    #random.seed(2)
    st_pt = (0.0,0.0)
    pts = [st_pt]
    for i in range(count):
        pts.append((random.random() - .5 + pts[-1][0], random.random() -.5 + pts[-1][1]))


    mnx = min(pts, key= lambda x: x[0])[0]
    mxx = max(pts, key= lambda x: x[0])[0]
    mny = min(pts, key=lambda x: x[1])[1]
    mxy = max(pts, key=lambda x: x[1])[1]

    mdx = statistics.median([pt[0] for pt in pts])
    mdy = statistics.median([pt[1] for pt in pts])

    pts = [((pt[0] - mnx) / (mxx - mnx), (pt[1] - mny) / (mxy - mny)) for pt in pts]
    mdx = statistics.median([pt[0] for pt in pts])
    mdy = statistics.median([pt[1] for pt in pts])
    pts = [(pt[0] - mdx , pt[1] - mdy) for pt in pts]
    return pts


def full_coreset_example(point_count, alpha, method, min_r=None):
    if min_r is None:
        min_r = alpha

    pts = boxed_trajectory(point_count)
    if method == "lifting":
        core_set_pts = pyscan.lifting_kernel(pts, math.sqrt(alpha))
    elif method == "grid":
        core_set_pts = pyscan.grid_kernel(pts, alpha)
    elif method == "even":
        core_set_pts = pyscan.even_sample_error(pts, alpha, False)
    elif method =="grid_alt":
        core_set_pts = pyscan.grid_trajectory(pts, alpha)
    elif method == "grid_direc":
        chord_l = math.sqrt(4 * alpha * min_r - 2 * alpha * alpha)
        core_set_pts = pyscan.grid_direc_kernel(pts, chord_l, alpha)
    elif method == "kernel":
        core_set_pts = pyscan.halfplane_kernel(pts, alpha)
    elif method == "graham":
        core_set_pts = pyscan.hull(pts)
    elif method == "dp":
        core_set_pts = pyscan.dp_compress(pts, alpha)
    else:
        return

    ax = plt.subplot()

    plot_tuples(ax, pts, "g")
    plot_points(ax, core_set_pts, "b")



def partial_coreset_example(traj_count, point_count, s, method):
    trajectories = [boxed_trajectory(point_count) for i in range(traj_count)]

    if method == "even":
        core_set_pts = pyscan.even_sample(trajectories, s, False)
    elif method == "uniform":
        core_set_pts = pyscan.uniform_sample(trajectories, s, False)
    elif method == "block":
        core_set_pts = pyscan.block_sample(trajectories, s, False)
    ax = plt.subplot()

    for traj in trajectories:
        plot_tuples(ax, traj, "g")

    plot_points(ax, core_set_pts, "b")


for method in ["even", "block", "uniform"]:
    partial_coreset_example(1, 20, 3, method)
    plt.savefig(method + "_partial.pdf")
    plt.clf()

for method in ["lifting", "grid", "even", "grid_alt", "grid_direc", "kernel", "dp", "graham"]:
    full_coreset_example(20, .1, method, min_r=.1)
    plt.savefig(method + "_full.pdf")
    plt.clf()
