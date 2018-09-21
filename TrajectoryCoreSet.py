import pyscan
import matplotlib.pyplot as plt
import random
import math

def plot_points(ax, pts, c):
    xs = []
    ys = []
    for pt in pts:
        xs.append(pt.getX())
        ys.append(pt.getY())
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

def halfplane_approx_example():
    st_pt = (0.0,0.0)
    pts = [st_pt]
    for i in range(100):
        pts.append(((random.random() - .5) / 10.0 + pts[-1][0], (random.random() - .5) / 10.0 + pts[-1][1]))

    new_pts = pyscan.approximate_hull(.00001, pts)

    ax = plt.subplot()
    plot_tuples(ax, pts, "b")
    plot_points(ax, new_pts, "r")
    plt.show()



def grid_approx_example():
    st_pt = (0.0,0.0)
    pts = [st_pt]
    for i in range(100):
        pts.append(((random.random() - .5) / 10.0 + pts[-1][0], (random.random() - .5) / 10.0 + pts[-1][1]))

    alpha = .001
    min_r = 10000
    chord_l = math.sqrt(2 * alpha * (2 * min_r - alpha))
    print(chord_l)
    new_pts = pyscan.grid_traj(pts, chord_l)
    print(len(new_pts))
    ax = plt.subplot()
    plot_tuples(ax, pts, "b")
    new_pt_list = []
    for key in new_pts:
        new_pt_list.extend(new_pts[key])
    print(new_pts)
    plot_points(ax, new_pt_list, "r")
    plt.show()


def grid_coreset_example():
    st_pt = (0.0,0.0)
    pts = [st_pt]
    for i in range(100):
        pts.append(((random.random() - .5) / 10.0 + pts[-1][0], (random.random() - .5) / 10.0 + pts[-1][1]))

    alpha = .01
    min_r = 1.0
    chord_l = math.sqrt(2 * alpha * (2 * min_r - alpha))
    print(chord_l)
    grid_pts = pyscan.grid_traj(pts, chord_l)
    grid_pts_list = []
    for key in grid_pts:
        grid_pts_list.extend(grid_pts[key])

    new_pts = pyscan.approx_traj_cells(pts, chord_l, alpha)

    ax = plt.subplot()

    new_pt_list = []
    for key in new_pts:
        new_pt_list.extend(new_pts[key])

    plot_tuples(ax, pts, "g")
    plot_points(ax, grid_pts_list, "b")
    plot_points(ax, new_pt_list, "r")
    plt.show()


def coreset_3d_example():
    st_pt = (0.0,0.0)
    pts = [st_pt]
    for i in range(100):
        pts.append(((random.random() - .5) / 10.0 + pts[-1][0], (random.random() - .5) / 10.0 + pts[-1][1]))

    alpha = .01

    grid_pts = pyscan.grid_traj(pts, .005)
    grid_pts_list = []
    for key in grid_pts:
        grid_pts_list.extend(grid_pts[key])
    core_set_pts = pyscan.core_set_3d_traj(grid_pts_list, alpha)
    ax = plt.subplot()

    plot_tuples(ax, pts, "g")
    plot_points(ax, core_set_pts, "r")
    plt.show()


coreset_3d_example()