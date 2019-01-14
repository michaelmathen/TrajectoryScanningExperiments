import csv
import numpy as np
import utils
import paths
import matplotlib.pyplot as plt
import time
import pyscan
import itertools
import math


def multiscale_disk(min_disk_r, max_disk_r, alpha, red_sample, blue_sample, net, disc, fast_disk):

    mx = -1
    curr_disk_r = min_disk_r
    reg = None
    while curr_disk_r < max_disk_r:

        chord_l = math.sqrt(4 * alpha * curr_disk_r - 2 * alpha * alpha)
        m_sample = [pyscan.grid_direc_kernel(pyscan.dp_compress(traj, alpha), chord_l, alpha) for traj in red_sample]
        b_sample = [pyscan.grid_direc_kernel(pyscan.dp_compress(traj, alpha), chord_l, alpha) for traj in blue_sample]
        pt_net = [pyscan.grid_direc_kernel(pyscan.dp_compress(traj, alpha), chord_l, alpha) for traj in net]
        m_sample = list(pyscan.trajectories_to_labels(m_sample))
        b_sample = list(pyscan.trajectories_to_labels(b_sample))
        net_set = list(pyscan.trajectories_to_labels(pt_net))
        #net_set = list(itertools.chain.from_iterable(pt_net))


        new_reg, new_mx = pyscan.max_disk_scale_labeled(net_set, m_sample, b_sample, fast_disk, curr_disk_r, disc)
        #print("Finished")
        #print("Should match {} {}".format(new_mx, pyscan.evaluate_range(new_reg, m_sample, b_sample, disc)))
        if new_mx > mx:
            reg = new_reg
            mx = new_mx
        curr_disk_r *= 2
    return reg, mx


def multiscale_disk_fixed(min_disk_r, max_disk_r, m_sample, b_sample, net_set, disc, fast_disk):


    #net_set = list(itertools.chain.from_iterable(net))

    mx = -1
    curr_disk_r = min_disk_r
    reg = None
    while curr_disk_r < max_disk_r:
        #print(net_set)
        new_reg, new_mx = pyscan.max_disk_scale_labeled(net_set, m_sample, b_sample, fast_disk, curr_disk_r, disc)
        if new_mx > mx:
            reg = new_reg
            mx = new_mx
        curr_disk_r *= 2
    return reg, mx


def testing_full_framework(
        red, blue,
        output_file,
        l_s, h_s, count,
        vparam="eps",
        eps=.01,
        alpha=.01,
        max_disk_r=.1,
        min_disk_r=.05,
        disc_name="disc",
        region_name="halfplane",
        sample_method="halfplane",
        fast_disk = True,
        two_level_sample=True,
        max_time=None):

    """
    How do I convert the trajectories over?
    1) Just sample evenly from the length.
    2) Choose points evenly
    3) Choose
    :param trajectories:
    :param l_s:
    :param h_s:
    :param count:
    :param vparam:
    :param eps:
    :param eps_r:
    :param r:
    :param q:
    :param disc_name:
    :param region_name:
    :param input_size:
    :return:
    """

    fieldnames = ["vparam", "disc", "region", "n", "s", "n_pts", "m_pts", "b_pts", "alpha", "time",
                  "m_disc", "m_disc_approx", "sample_method"]
    with open(output_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for i in np.logspace(l_s, h_s, count):

            if vparam == "eps":
                eps = i
            elif vparam == "alpha":
                alpha = i
            n = 1 / eps
            s = 1 / (2 * eps * eps)
            n = int(round(n) + .1)
            s = int(round(s) + .1)

            disc = utils.disc_to_func(disc_name)

            red_sample = pyscan.my_sample(red, s)
            blue_sample = pyscan.my_sample(blue, s)
            if two_level_sample:
                red_net = pyscan.my_sample(red, n)
                blue_net = pyscan.my_sample(blue, n)
            else:
                red_net = red_sample
                blue_net = blue_sample

            net = red_net + blue_net

            print("Running: {} {}".format(n, s))

            start_time = time.time()

            if region_name == "multiscale_disk":
                reg, mx = multiscale_disk(min_disk_r, max_disk_r, alpha, red_sample, blue_sample, net, disc, fast_disk)
                m_sample, b_sample, net_set = [], [], []
            else:
                if sample_method == "halfplane":
                    m_sample = [pyscan.halfplane_kernel([pyscan.Point(pt[0], pt[1], 1.0) for pt in traj], alpha) for traj in red_sample]
                    b_sample = [pyscan.halfplane_kernel([pyscan.Point(pt[0], pt[1], 1.0) for pt in traj], alpha) for traj in blue_sample]
                    pt_net = [pyscan.halfplane_kernel([pyscan.Point(pt[0], pt[1], 1.0) for pt in traj], alpha) for traj in net]
                elif sample_method == "dp":
                    m_sample = [pyscan.dp_compress(traj, alpha) for traj in red_sample]
                    b_sample = [pyscan.dp_compress(traj, alpha) for traj in blue_sample]
                    pt_net = [pyscan.dp_compress(traj, alpha) for traj in net]
                elif sample_method == "hull":
                    m_sample = [pyscan.convex_hull([pyscan.Point(pt[0], pt[1], 1.0) for pt in traj]) for traj in red_sample]
                    b_sample = [pyscan.convex_hull([pyscan.Point(pt[0], pt[1], 1.0) for pt in traj]) for traj in blue_sample]
                    pt_net = [pyscan.convex_hull([pyscan.Point(pt[0], pt[1], 1.0) for pt in traj]) for traj in net]
                elif sample_method is None:
                    #just takes the waypoints.
                    m_sample = [[pyscan.Point(pt[0], pt[1], 1.0) for pt in traj] for traj in red_sample]
                    b_sample = [[pyscan.Point(pt[0], pt[1], 1.0) for pt in traj] for traj in blue_sample]
                    pt_net = [[pyscan.Point(pt[0], pt[1], 1.0) for pt in traj] for traj in net]
                elif sample_method == "grid":
                    m_sample = [pyscan.grid_kernel(traj, alpha) for traj in red_sample]
                    b_sample = [pyscan.grid_kernel(traj, alpha) for traj in blue_sample]
                    pt_net = [pyscan.grid_kernel(traj, alpha) for traj in net]
                elif sample_method == "lifting":
                    m_sample = [pyscan.lifting_kernel(traj, alpha) for traj in red_sample]
                    b_sample = [pyscan.lifting_kernel(traj, alpha) for traj in blue_sample]
                    pt_net = [pyscan.lifting_kernel(traj, alpha) for traj in net]
                elif sample_method == "grid_direc":
                    chord_l = math.sqrt(4 * alpha * min_disk_r - 2 * alpha * alpha)
                    m_sample = [pyscan.grid_direc_kernel(pyscan.dp_compress(traj, alpha), chord_l, alpha) for traj in
                                red_sample]
                    b_sample = [pyscan.grid_direc_kernel(pyscan.dp_compress(traj, alpha), chord_l, alpha) for traj in
                                blue_sample]
                    pt_net = [pyscan.grid_direc_kernel(pyscan.dp_compress(traj, alpha), chord_l, alpha) for traj in net]
                elif sample_method == "even":
                    m_sample = [pyscan.even_sample_error(traj, alpha, False) for traj in red_sample]
                    b_sample = [pyscan.even_sample_error(traj, alpha, False) for traj in blue_sample]
                    pt_net = [pyscan.even_sample_error(traj, alpha, False) for traj in net]
                else:
                    return

                if region_name == "multiscale_disk_fixed":
                    m_sample = list(pyscan.trajectories_to_labels(red_sample))
                    b_sample = list(pyscan.trajectories_to_labels(blue_sample))
                    net_set = list(pyscan.trajectories_to_labels(net))
                    reg, mx = multiscale_disk_fixed(min_disk_r, max_disk_r, m_sample, b_sample, net_set, disc, fast_disk)
                else:
                    m_sample = list(pyscan.trajectories_to_labels(m_sample))
                    b_sample = list(pyscan.trajectories_to_labels(b_sample))
                    net_set = list(itertools.chain.from_iterable(pt_net))
                    if region_name == "halfplane":
                        reg, mx = pyscan.max_halfplane_labeled(net_set, m_sample, b_sample, disc)
                    elif region_name == "disk":
                        reg, mx = pyscan.max_disk_labeled(net_set, m_sample, b_sample, disc)
                    elif region_name == "rectangle":
                        reg, mx = pyscan.max_rect_labeled(n, max_disk_r, m_sample, b_sample, disc)
                    elif region_name == "rectangle_scale":
                        reg, mx = pyscan.max_rect_labeled(n, max_disk_r, alpha, net_set, m_sample, b_sample, disc)
                    else:
                        return

            end_time = time.time()
            actual_mx = pyscan.evaluate_range_trajectory(reg, red, blue, disc)
            row = {"vparam": vparam,
                   "disc": disc_name,
                   "region": region_name,
                   "n": n, "s": s,
                   "n_pts": len(net_set), "m_pts":len(m_sample), "b_pts":len(b_sample),
                   "alpha":alpha,
                   "time": end_time - start_time,
                   "m_disc_approx": mx,
                   "m_disc": actual_mx,
                   "sample_method": sample_method}
            writer.writerow(row)
            f.flush()
            print(row)
            if max_time is not None and end_time - start_time > max_time:
                return


def generate_disk_sets(fname, r, p, q, min_r, max_r):

    disc = utils.disc_to_func("disc")
    trajectories = paths.read_dong_csv("/data/Dong_sets/Trajectory_Sets/samples/{}.tsv".format(fname))

    while True:
        red, blue, planted_reg , planted_mx = pyscan.plant_full_disk(trajectories, r, p, q, disc)
        print(min_r, planted_reg.get_radius(), max_r)
        if min_r < planted_reg.get_radius() < max_r:
            break
    return red, blue

def generate_halfplane_sets(fname, r, p, q):
    disc = utils.disc_to_func("disc")
    trajectories = paths.read_dong_csv("/data/Dong_sets/Trajectory_Sets/samples/{}.tsv".format(fname))
    red, blue, planted_reg, planted_mx = pyscan.plant_full_halfplane(trajectories, r, p, q, disc)
    return red, blue


