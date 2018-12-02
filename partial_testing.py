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
        planted_points=None,
        input_size=10000,
        fast_halfplane = False):

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
    if fast_halfplane:
        output_file = "Fast_{}_{}_{}_partial_discrepancy.csv".format(disc_name, region_name, sample_method)
    else:
        output_file = "{}_{}_{}_partial_discrepancy.csv".format(disc_name, region_name, sample_method)

    fieldnames = ["vparam", "input_size", "disc", "region", "n", "s", "r", "p", "q", "time",
                  "m_disc",
                  "m_disc_approx", "sample_method"]
    with open(output_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        disc = utils.disc_to_func(disc_name)
        scan = utils.range_to_func(region_name)


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
            #if vparam != "eps":
            if planted_points is None:
                if region_name == "disk":
                    red, blue, _ = pyscan.plant_partial_disk(trajectories, r, p, q, eps_r, disc)
                elif region_name == "halfplane":
                    red, blue, _ = pyscan.plant_partial_halfplane(trajectories, r, p, q, eps_r, disc)
            else:
                red, blue = planted_points

            et = time.time()
            print("Time to plant region {}".format(et - st))


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

            start_time = time.time()
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
            f.flush()

if __name__ == "__main__":

    #trajectories = paths.read_geolife_files(100)
    trajectories = paths.read_dong_csv("/data/Dong_sets/Trajectory_Sets/samples/bjtaxi_samples_1m.tsv")
    print(len(trajectories))
    r = .01
    q = .2
    p = .5
    eps_r = .001
    disc = utils.disc_to_func("disc")
    for region_name in ["disk"]:
        if region_name == "disk":
            red, blue, _, _ = pyscan.plant_partial_disk(trajectories, r, p, q, eps_r, disc)
        elif region_name == "halfplane":
            red, blue, _, _ = pyscan.plant_partial_halfplane(trajectories, r, p, q, eps_r, disc)

        for sample in ["block", "even", "uniform"]:
            disc = utils.disc_to_func("disc")
            testing_partial_framework(trajectories, -1, -2, 40, r=r, q=q, p=p,
                                    region_name=region_name,
                                    sample_method=sample,
                                    planted_points=(red, blue), fast_halfplane=False)
