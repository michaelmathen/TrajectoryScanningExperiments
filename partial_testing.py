import csv
import numpy as np
import utils
import paths
import matplotlib.pyplot as plt
import time
import pyscan


def testing_partial_framework(
        trajectories,
        l_s, h_s, count,
        vparam="eps",
        eps=.01,
        eps_r=.01,
        r=.04,
        p=0.5,
        q=.2,
        error_thresh=3,
        disc_name="disc",
        region_name="disk",
        sample_method="block",
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
    output_file = "{}_{}_partial_discrepancy.csv".format(disc_name, region_name)

    fieldnames = ["vparam", "input_size", "disc", "region", "n", "s", "r", "p", "q", "time",
                  "m_disc",
                  "m_disc_approx", "sample_method"]
    with open(output_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        disc = utils.disc_to_func(disc_name)
        scan = utils.range_to_func(region_name)

        if vparam == "eps":
            if region_name == "disk":
                red, blue, _ = pyscan.plant_trajectory_disk(trajectories, r, p, q, disc)
            elif region_name == "halfplane":
                red, blue, _ = pyscan.plant_trajectory_halfplane(trajectories, r, p, q, disc)

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

            st = time.time()
            if vparam != "eps":
                if region_name == "disk":
                    red, blue, _ = pyscan.plant_trajectory_disk(trajectories, r, p, q, disc)
                elif region_name == "halfplane":
                    red, blue, _ = pyscan.plant_trajectory_halfplane(trajectories, r, p, q, disc)

            et = time.time()
            print("Time to plant region {}".format(et - st))

            start_time = time.time()
            if sample_method == "block":
                m_sample = pyscan.block_sample(red, s, False)
                b_sample = pyscan.block_sample(blue, s, False)

                net_set1 = pyscan.block_sample(red, n, False)
                net_set2 = pyscan.block_sample(blue, n, False)
            elif sample_method == "even":
                m_sample = pyscan.even_sample(red, s, False)
                b_sample = pyscan.even_sample(blue, s, False)

                net_set1 = pyscan.even_sample(red, n, False)
                net_set2 = pyscan.even_sample(blue, n, False)
            elif sample_method == "uniform":
                m_sample = pyscan.uniform_sample(red, s, False)
                b_sample = pyscan.uniform_sample(blue, s, False)

                net_set1 = pyscan.uniform_sample(red, n, False)
                net_set2 = pyscan.uniform_sample(blue, n, False)

            net_set = net_set1 + net_set2

            # def plot_tuples(ax, pts, c):
            #     xs = []
            #     ys = []
            #     for pt in pts:
            #         xs.append(pt[0])
            #         ys.append(pt[1])
            #     ax.plot(xs, ys, color=c)
            # ax = plt.subplot()
            # for traj in trajectories:
            #     plot_tuples(ax, traj, "g")
            # ax.scatter([x[0] for x in m_sample], [x[1] for x in m_sample], c="r")
            # ax.scatter([x[0] for x in b_sample], [x[1] for x in b_sample], c="b")
            # plt.show()
            m_sample = pyscan.to_weighted(m_sample)
            b_sample = pyscan.to_weighted(b_sample)
            reg, mx = scan(net_set, m_sample, b_sample, disc)
            end_time = time.time()

            s_prime = int(10 ** (2 * error_thresh) + .5)
            print(s_prime)
            m_sample_prime = pyscan.uniform_sample(red, s_prime, False)
            b_sample_prime = pyscan.uniform_sample(blue, s_prime, False)
            #print(m_sample_prime)
            actual_mx = pyscan.evaluate_range(reg, pyscan.to_weighted(m_sample_prime), pyscan.to_weighted(b_sample_prime), disc)

            row = {"vparam": vparam,
                   "input_size": input_size,
                   "disc": disc_name,
                   "region": region_name,
                   "n": n, "s": s, "r": r, "q": q,"p":p,
                   "time": end_time - start_time,
                   "m_disc_approx": mx,
                   "m_disc": actual_mx,
                   "sample_method": sample_method}
            writer.writerow(row)
            print(row)


if __name__ == "__main__":

    trajectories = paths.read_geolife_files(1000)
    for range in ["halfspace", "disk"]:
        for sample in ["block", "even", "uniform"]:
            testing_partial_framework(trajectories, -1, -2, 10, r=.01, q=0.0, region_name="disk", sample_method=sample)