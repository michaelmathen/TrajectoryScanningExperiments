import matplotlib.pyplot as plt
import matplotlib.patches as patch
import pyscan
import utils
import numpy as np
import paths



def plot_tuples(ax, pts, c):
    xs = []
    ys = []
    for pt in pts:
        xs.append(pt[0])
        ys.append(pt[1])
    ax.plot(xs, ys, color=c)


def plot_partial_trajectories(trajectories, r, p, q, eps_r):
    disc = utils.disc_to_func("disc")
    red, blue, mx_disk, _ = pyscan.plant_partial_disk(trajectories, r, p, q, eps_r, disc)


    ax = plt.subplot()
    for traj in red:
        plot_tuples(ax, traj.get_pts(), "r")

    for traj in blue:
        plot_tuples(ax, traj.get_pts(), "b")

    actor = plt.Circle((mx_disk.get_origin()[0], mx_disk.get_origin()[1]), mx_disk.get_radius())
    ax.add_artist(actor)
    plt.show()


def plot_full_trajectories(trajectories, r, p, q):
    disc = utils.disc_to_func("disc")
    red, blue, mx_disk, _ = pyscan.plant_full_disk(trajectories, r, p, q, disc)


    ax = plt.subplot()
    for traj in blue:
        plot_tuples(ax, traj, "b")
    for traj in red:
        plot_tuples(ax, traj, "r")



    actor = plt.Circle((mx_disk.get_origin()[0], mx_disk.get_origin()[1]), mx_disk.get_radius())
    print(mx_disk.get_origin(), mx_disk.get_radius())
    ax.add_artist(actor)
    plt.show()



def plot_plane_partial_trajectories(trajectories, r, p, q, eps_r):
    disc = utils.disc_to_func("disc")
    red, blue, mx_plane, _ = pyscan.plant_partial_halfplane(trajectories, r, p, q, eps_r, disc)


    ax = plt.subplot()
    for traj in red:
        plot_tuples(ax, traj.get_pts(), "r")

    for traj in blue:
        plot_tuples(ax, traj.get_pts(), "b")
    xs = np.arange(0, 1, .01)


    ys = (-1 - mx_plane.get_coords()[0] * xs) * 1 / mx_plane.get_coords()[1]
    ax.plot(xs, ys, color="g")
    plt.show()


def plot_plane_full_trajectories(trajectories, r, p, q):
    disc = utils.disc_to_func("disc")
    red, blue, mx_plane, _ = pyscan.plant_full_halfplane(trajectories, r, p, q, disc)


    ax = plt.subplot()
    for traj in blue:
        plot_tuples(ax, traj, "b")
    for traj in red:
        plot_tuples(ax, traj, "r")

    xs = np.arange(0, 1, .01)
    ys = (-1 - mx_plane.get_coords()[0] * xs) * 1/ mx_plane.get_coords()[1]
    ax.plot(xs, ys, color="g")
    plt.show()


trajectories = paths.read_geolife_files(1000)
#trajectories = clean(trajectories)
#print(len(trajectories))
plot_plane_full_trajectories(trajectories, .25, 1.0, 0.0)

