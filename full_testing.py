import csv
import numpy as np
import utils
import paths
import matplotlib.pyplot as plt
import time
import pyscan
import itertools


def testing_partial_framework(
        trajectories,
        l_s, h_s, count,
        vparam="eps",
        eps=.01,
        eps_r=.01,
        r=.04,
        p=0.5,
        q=.2,
        alpha=.01,
        error_thresh=3,
        max_disk_size=.1,
        disc_name="disc",
        region_name="halfplane",
        sample_method="halfplane",
        input_size=10000):

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
    output_file = "{}_{}_full_discrepancy.csv".format(disc_name, region_name)

    fieldnames = ["vparam", "input_size", "disc", "region", "n", "s", "r", "p", "q", "alpha", "time",
                  "m_disc", "n_pts", "sample_size",
                  "m_disc_approx", "sample_method"]
    with open(output_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for i in np.logspace(l_s, h_s, count):

            if vparam == "eps":
                eps = i
            elif vparam == "r":
                r = i
            elif vparam == "q":
                q = i
            n = 1 / eps
            s = 1 / (2 * eps * eps)
            n = int(round(n) + .1)
            s = int(round(s) + .1)

            disc = utils.disc_to_func(disc_name)

            scan = utils.range_to_lfunc(region_name)
            st = time.time()
            if "disk" in region_name:
                red, blue, actual_mx = pyscan.plant_trajectory_disk(trajectories, r, p, q, disc)
            elif region_name == "halfplane":
                red, blue, actual_mx = pyscan.plant_trajectory_halfplane(trajectories, r, p, q, disc)
            print("got here")
            red_sample = pyscan.my_sample(red, s)
            blue_sample = pyscan.my_sample(blue, s)
            red_net = pyscan.my_sample(red, s)
            blue_net = pyscan.my_sample(blue, n)
            net = red_net + blue_net

            et = time.time()
            print("Time to plant region {}".format(et - st))

            start_time = time.time()
            if sample_method == "halfplane":
                m_sample = [pyscan.halfplane_kernel(traj, alpha) for traj in red_sample]
                b_sample = [pyscan.halfplane_kernel(traj, alpha) for traj in blue_sample]
                net = [pyscan.halfplane_kernel(traj, alpha) for traj in net]
            elif sample_method == "grid":
                m_sample = [pyscan.grid_kernel(traj, alpha) for traj in red_sample]
                b_sample = [pyscan.grid_kernel(traj, alpha) for traj in blue_sample]
                net = [pyscan.grid_kernel(traj, alpha) for traj in net]
            elif sample_method == "grid_direc":
                m_sample = [pyscan.grid_direc_kernel(traj, alpha) for traj in red_sample]
                b_sample = [pyscan.grid_direc_kernel(traj, alpha) for traj in blue_sample]
                net = [pyscan.grid_direc_kernell(traj, alpha) for traj in net]
            elif sample_method == "lifting":
                m_sample = [pyscan.lifting_kernel(traj, alpha) for traj in red_sample]
                b_sample = [pyscan.lifting_kernel(traj, alpha) for traj in blue_sample]
                net = [pyscan.lifting_kernel(traj, alpha) for traj in net]
            elif sample_method == "dp":
                m_sample = [pyscan.dp_compress(traj, alpha) for traj in red_sample]
                b_sample = [pyscan.dp_compress(traj, alpha) for traj in blue_sample]
                net = [pyscan.dp_compress(traj, alpha) for traj in net]

            m_sample = list(pyscan.trajectories_to_labels(m_sample))
            b_sample = list(pyscan.trajectories_to_labels(b_sample))
            net_set = list(itertools.chain.from_iterable(net))

            if region_name == "halfplane":
                reg, mx = pyscan.max_halfplane_labeled(net_set, m_sample, b_sample, disc)
            elif region_name == "disk":
                reg, mx = pyscan.max_disk_labeled(net_set, m_sample, b_sample, disc)
            elif region_name == "small_disk":
                reg, mx = pyscan.max_disk_scale_labeled(net_set, m_sample, b_sample, alpha, int(1 / max_disk_size + 1), disc)

            end_time = time.time()

            s_prime = int(10 ** (2 * error_thresh) + .5)
            print(s_prime)

            row = {"vparam": vparam,
                   "input_size": input_size,
                   "disc": disc_name,
                   "region": region_name,
                   "n": n, "s": s, "n_pts":len(net_set), "sample_size": len(m_sample) + len(b_sample),
                   "r": r, "q": q,"p":p,
                   "alpha":alpha,
                   "time": end_time - start_time,
                   "m_disc_approx": mx,
                   "m_disc": actual_mx,
                   "sample_method": sample_method}
            writer.writerow(row)
            print(row)


if __name__ == "__main__":

    trajectories = paths.read_geolife_files(1000)
    testing_partial_framework(trajectories, -1, -3, 10, r=.01, q=0.0, region_name="small_disk")