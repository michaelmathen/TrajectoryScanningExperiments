import matplotlib.pyplot as plt
import matplotlib.patches as patch
import pyscan
import utils
import numpy as np
import paths
import itertools
import math
import random

def plot_tuples(ax, pts, c):
    xs = []
    ys = []
    for pt in pts:
        xs.append(pt[0])
        ys.append(pt[1])
    ax.scatter(xs, ys, color=c, marker='.')


def plot_line(ax, pts, c):
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
    rpts = list(itertools.chain.from_iterable([traj.get_pts() for traj in red]))
    rpts = random.sample(rpts, 500)

    plot_tuples(ax, list(rpts), "r")


    bpts = list(itertools.chain.from_iterable([traj.get_pts() for traj in blue]))
    bpts = random.sample(bpts, 500)
    plot_tuples(ax, bpts, "b")

    actor = plt.Circle((mx_disk.get_origin()[0], mx_disk.get_origin()[1]), mx_disk.get_radius(), edgecolor="k", linewidth=2, fill=False)
    ax.add_artist(actor)
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 1.0)
    plt.tight_layout()
    plt.axis('off')
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


def plot_full_trajectories_intersection_check(trajectories, r):
    disc = utils.disc_to_func("disc")
    red, blue, mx_disk, _ = pyscan.plant_full_disk(trajectories, r, .5, .5, disc)

    #mx_disk = pyscan.Disk(mx_disk.get_origin()[0], mx_disk.get_origin()[1],  .01)
    ax = plt.subplot()

    blue_c = 0
    for traj in trajectories:
        if pyscan.Trajectory(traj).intersects_disk(mx_disk):
            blue_c += 1
            plot_line(ax, traj, "b")
        else:
             plot_line(ax, traj, "r")
    print(blue_c / len(trajectories), r)

    actor = plt.Circle((mx_disk.get_origin()[0], mx_disk.get_origin()[1]), mx_disk.get_radius())
    print(mx_disk.get_origin(), mx_disk.get_radius())
    ax.plot(mx_disk.get_origin()[0], mx_disk.get_origin()[1], marker='o')
    ax.add_artist(actor)
    plt.show()


def plot_full_trajectories_intersection_rect(trajectories, r):
    disc = utils.disc_to_func("disc")
    red, blue, rect, _ = pyscan.plant_full_square(trajectories, r, 0.5, 0.5, disc)

    ax = plt.subplot()

    blue_c = 0
    for traj in trajectories:
        if rect.intersects_trajectory(traj):
            blue_c += 1
            plot_line(ax, traj, "b")
        else:
            plot_line(ax, traj, "r")
    print(blue_c / len(trajectories), r)

    actor = plt.Rectangle((rect.upX(), rect.upY()), rect.upX() - rect.lowX(), rect.upY() - rect.lowY())
    ax.add_artist(actor)
    plt.show()


def plot_plane_partial_trajectories(trajectories, r, p, q, eps_r):
    disc = utils.disc_to_func("disc")
    red, blue, mx_plane, _ = pyscan.plant_partial_halfplane(trajectories, r, p, q, eps_r, disc)


    print(len(red), len(blue))
    ax = plt.subplot()
    for traj in red:
        plot_tuples(ax, traj.get_pts(), "r")

    for traj in blue:
        plot_tuples(ax, traj.get_pts(), "b")

    xs = np.arange(0, 1, .01)
    ys = (-1 - mx_plane.get_coords()[0] * xs) * 1 / mx_plane.get_coords()[1]
    ax.plot(xs, ys, color="g")
    plt.show()


def plot_plane_flux_trajectories(trajectories, r, q, eps_r):

    st_pts, end_pts = pyscan.trajectories_to_flux(trajectories)

    st_pts = [pyscan.WPoint(1.0, float(p[0]), float(p[1]), 1.0) for p in st_pts]
    end_pts = [pyscan.WPoint(1.0, float(p[0]), float(p[1]), 1.0) for p in end_pts]

    scan = utils.range_to_func("halfplane")
    red, blue, mx_plane = pyscan.paired_plant_region(st_pts, end_pts, r, q, eps_r, scan)
    ax = plt.subplot()
    ax.scatter([x[0] for x in red], [x[1] for x in red], c="r")
    ax.scatter([x[0] for x in blue], [x[1] for x in blue], c="b")

    xs = np.arange(0, 1, .01)


    ys = (-1 - mx_plane.get_coords()[0] * xs) * 1 / mx_plane.get_coords()[1]
    ax.plot(xs, ys, color="g")
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 1.0)
    plt.show()

def plot_disk_flux_trajectories(trajectories, r, q, eps_r):

    st_pts, end_pts = pyscan.trajectories_to_flux(trajectories)

    st_pts = [pyscan.WPoint(1.0, float(p[0]), float(p[1]), 1.0) for p in st_pts]
    end_pts = [pyscan.WPoint(1.0, float(p[0]), float(p[1]), 1.0) for p in end_pts]

    scan = utils.range_to_func("disk")
    red, blue, mx_disk = pyscan.paired_plant_region(st_pts, end_pts, r, q, eps_r, scan)
    ax = plt.subplot()
    actor = plt.Circle((mx_disk.get_origin()[0], mx_disk.get_origin()[1]), mx_disk.get_radius(), fill=False)
    ax.add_artist(actor)

    ax.scatter([x[0] for x in red], [x[1] for x in red], c="r")
    ax.scatter([x[0] for x in blue], [x[1] for x in blue], c="b")

    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 1.0)
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


def testing(trajectories, alpha, max_r):

    curr_r = alpha
    while curr_r <= max_r:
        chord_l = math.sqrt(4 * alpha * curr_r -  2 * alpha * alpha)
        sample = [pyscan.grid_direc_kernel(pyscan.dp_compress(traj, alpha), chord_l, alpha) for traj in trajectories]
        pts = list(itertools.chain.from_iterable(sample))
        print("Grid Directional radius = {0:.4f} : {1:.4f} ".format(curr_r * 3000, len(pts) / len(trajectories)))
        curr_r *= 2

    sample = [pyscan.grid_kernel(pyscan.dp_compress(traj, alpha), alpha) for traj in trajectories]
    pts = list(itertools.chain.from_iterable(sample))
    print("Grid : {0:.2f}".format(len(pts) / len(trajectories)))

    sample = [pyscan.halfplane_kernel(pyscan.dp_compress(traj, alpha), alpha) for traj in trajectories]
    pts = list(itertools.chain.from_iterable(sample))
    print("Halfplane : {0:.2f}".format(len(pts) / len(trajectories)))

    sample = [pyscan.dp_compress(traj, alpha) for traj in trajectories]
    pts = list(itertools.chain.from_iterable(sample))
    print("DP : {0:.2f}".format(len(pts) / len(trajectories)))

    # sample = [pyscan.lifting_kernel(pyscan.dp_compress(traj, alpha), .01) for traj in trajectories]
    # pts = list(itertools.chain.from_iterable(sample))
    # print("Lifting : {0:.2f}".format(len(pts) / len(trajectories)))

    sample = [pyscan.convex_hull([pyscan.Point(pt[0], pt[1], 1.0) for pt in traj]) for traj in trajectories]
    pts = list(itertools.chain.from_iterable(sample))
    print("Hull : {0:.2f}".format(len(pts) / len(trajectories)))

    sample = [pyscan.even_sample_error([pyscan.Point(pt[0], pt[1], 1.0) for pt in traj], alpha, False) for traj in trajectories]
    pts = list(itertools.chain.from_iterable(sample))
    print("Even : {0:.2f}".format(len(pts) / len(trajectories)))

    # sample = [pyscan.grid_direc_kernel(pyscan.dp_compress(traj, alpha), chord_l, alpha) for traj in trajectories]
    # sample = [pyscan.grid_direc_kernel(pyscan.dp_compress(traj, alpha), chord_l, alpha) for traj in trajectories]

    # ax = plt.subplot()
    # ax.scatter([x[0] for x in pts], [x[1] for x in pts], c="r")
    # ax.set_xlim(0, 1.0)
    # ax.set_ylim(0, 1.0)
    # plt.show()

def test_halfspace_error(trajectories, core_sets):

    for (traj, traj_approx) in zip(trajectories, core_sets):
        (reg, curr_error) = pyscan.coreset_error_halfplane(traj, traj_approx)
        # if math.isnan(curr_error) or math.isinf(curr_error) or curr_error > 1.0:
        #     print(curr_error)
        #     print(traj)

        yield curr_error

def post_process_error(point_set):
    lset = list(point_set)
    return sorted(lset)[int(.9 * len(lset))]

def test_disk_error(trajectories, core_sets, min_radius, max_radius):
    for (traj, traj_approx) in zip(trajectories, core_sets):
        (reg, curr_error) = pyscan.coreset_error_disk(traj, min_radius, max_radius, traj_approx)
        print(curr_error)
        yield curr_error

# def test_halfspace_error(trajectories, core_sets):
#     max_error = 0
#     for (traj, traj_approx) in zip(trajectories, core_sets):
#         max_error = max(pyscan.coreset_error_halfplane(traj, traj_approx), max_error)
#     return max_error


def testing_geometric_error(trajectories, alpha, max_r, count):
    random.shuffle(trajectories)
    trajectories = trajectories[:30]
    print("got here")
    curr_r = alpha
    while curr_r < max_r:
        chord_l = math.sqrt(4 * alpha * curr_r -  2 * alpha * alpha)
        sample = [pyscan.grid_direc_kernel(pyscan.dp_compress(traj, alpha), chord_l, alpha) for traj in trajectories]
        for error, traj in zip(test_halfspace_error(trajectories, sample), trajectories):
            if error > alpha:
                print(error)
                print(traj)
        #print("Grid Direc Error {} {}".format(i, post_process_error(test_halfspace_error(trajectories, sample))))


    sample = [pyscan.grid_kernel(pyscan.dp_compress(traj, alpha), alpha) for traj in trajectories]
    print("Grid : {}".format(post_process_error(test_halfspace_error(trajectories, sample))))

    sample = [pyscan.halfplane_kernel(pyscan.dp_compress(traj, alpha), alpha) for traj in trajectories]
    print("Halfplane : {}".format( post_process_error(test_halfspace_error(trajectories, sample))))

    sample = [pyscan.dp_compress(traj, alpha) for traj in trajectories]
    print("DP Error: {}".format(post_process_error(test_halfspace_error(trajectories, sample))))

    sample = [pyscan.convex_hull([pyscan.Point(pt[0], pt[1], 1.0) for pt in traj]) for traj in trajectories]
    print("Hull Error: {}".format( post_process_error(test_halfspace_error(trajectories, sample))))

    sample = [pyscan.even_sample_error(pyscan.dp_compress(traj, alpha), alpha) for traj in trajectories]
    print("Even : {}".format(post_process_error(test_halfspace_error(trajectories, sample))))

    # sample = [pyscan.lifting_kernel(pyscan.dp_compress(traj, alpha), .01) for traj in trajectories]
    # pts = list(itertools.chain.from_iterable(sample))
    # print("Lifting : {0:.2f}".format(len(pts) / len(trajectories)))

def testing_disk_geometric_error(trajectories, alpha, max_r, count):
    random.shuffle(trajectories)
    trajectories = trajectories[:30]
    for i in np.linspace(alpha, max_r, count):
        chord_l = math.sqrt(4 * alpha * i -  2 * alpha * alpha)
        sample = [pyscan.grid_direc_kernel(pyscan.dp_compress(traj, alpha), chord_l, alpha) for traj in trajectories]
        print("Grid Direc Error {} {}".format(50 * i, 50 * post_process_error(test_disk_error(trajectories, sample, alpha, max_r))))


    sample = [pyscan.grid_kernel(pyscan.dp_compress(traj, alpha), alpha) for traj in trajectories]
    print("Grid : {}".format(50 * post_process_error(test_disk_error(trajectories, sample, alpha, max_r))))

    sample = [pyscan.halfplane_kernel(pyscan.dp_compress(traj, alpha), alpha) for traj in trajectories]
    print("Halfplane : {}".format(50 * post_process_error(test_disk_error(trajectories, sample, alpha, max_r))))

    sample = [pyscan.dp_compress(traj, alpha) for traj in trajectories]
    print("DP Error: {}".format(50 * post_process_error(test_disk_error(trajectories, sample, alpha, max_r))))

#trajectories = paths.read_geolife_files(1000)
trajectories = paths.read_dong_csv("/data/Dong_sets/Trajectory_Sets/samples/osm_eu_sample_10k_nw.tsv")

#trajectories = clean(trajectories)
#print(len(trajectories))
#plot_plane_full_trajectories(trajectories, .25, 1.0, 0.0)

r=.02
p =.5
q= 1.0
eps_r=.01

#plot_full_trajectories_intersection_rect(trajectories, r)
# OSM EU 1/30000, 1/300
# BJTAXI 1/500 1/5
testing(trajectories, 1/30000, 1/300)
#testing_geometric_error(trajectories, 1/500, 1/10, 3)



#testing_disk_geometric_error(trajectories,.1, .01, 5)
#plot_partial_trajectories(trajectories, r, p, q, eps_r)
#plot_plane_partial_trajectories(trajectories, r, p, q, eps_r)

