import csv
import numpy as np
import utils
import paths
import matplotlib.pyplot as plt
import time
import pyscan
import itertools
import math


def testing_full_framework(
        trajectories,
        l_s, h_s, count,
        vparam="eps",
        eps=.01,
        r=.04,
        p=0.5,
        q=.2,
        alpha=.01,
        planted_points=None,
        actual_mx=None,
        max_disk_r=.1,
        min_disk_r=.05,
        disc_name="disc",
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
    output_file = "{0}_multi_disk_{1:.2f}_{2:.2f}_full_discrepancy.csv".format(disc_name,  min_disk_r, max_disk_r)

    fieldnames = ["vparam", "input_size", "region", "disc", "n", "s", "r", "p", "q", "alpha", "time",
                  "m_disc",
                  "m_disc_approx"]
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


            st = time.time()
            if planted_points is None:
                red, blue, _, actual_mx = pyscan.plant_full_disk(trajectories, r, p, q, disc)
            else:
                red, blue = planted_points

            red_sample = pyscan.my_sample(red, s)
            blue_sample = pyscan.my_sample(blue, s)
            red_net = pyscan.my_sample(red, n)
            blue_net = pyscan.my_sample(blue, n)
            net = red_net + blue_net

            et = time.time()
            print("Time to plant region {}".format(et - st))

            start_time = time.time()
            reg, mx = pyscan.max_disk_traj_grid(net, [pyscan.WTrajectory(1.0, traj) for traj in red_sample],
                                                    [pyscan.WTrajectory(1.0, traj) for traj in blue_sample], min_disk_r, max_disk_r, disc)
            end_time = time.time()

            row = {"vparam": vparam,
                   "input_size": input_size,
                   "disc": disc_name,
                   "region": "multi_disk",
                   "n": n, "s": s,
                   "r": r, "q": q,"p":p,
                   "alpha":alpha,
                   "time": end_time - start_time,
                   "m_disc_approx": mx,
                   "m_disc": actual_mx}
            writer.writerow(row)
            f.flush()
            print(row)


if __name__ == "__main__":

    trajectories = paths.read_geolife_files(10000)
    # print(len(trajectories))

    r = .05
    p = .5
    q = .8
    alpha = .01
    max_r = .05
    disc = utils.disc_to_func("disc")

    red, blue, _ ,actual_mx = pyscan.plant_full_disk(trajectories, r, p, q, disc)
    #
    for min_r in np.logspace(math.log(alpha, 10), math.log(max_r, 10), num=5, endpoint=False):
        print(min_r, max_r)
        testing_full_framework(trajectories, -1, -1.8, 20, r=r, q=q, p=p,
                       planted_points=(red,blue),
                        actual_mx=actual_mx,
                       min_disk_r = min_r,
                       max_disk_r = max_r)
