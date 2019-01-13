import csv
import numpy as np
import utils
import paths
import matplotlib.pyplot as plt
import time
import pyscan
import math

def testing_partial_framework(
        red, blue,
        output_file,
        l_s, h_s, count,
        r=.04,
        p=0.5,
        q=.2,
        error_thresh=3,
        two_level_sample=True,
        ham_sample=True,
        disc_name="disc",
        region_name="disk",
        sample_method="block",
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
    fieldnames = ["disc", "region", "n", "s", "r", "p", "q", "time",
                  "m_disc",
                  "m_disc_approx", "sample_method"]
    with open(output_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        disc = utils.disc_to_func(disc_name)
        s_prime = int(10 ** (2 * error_thresh) + .5)
        m_sample_prime = pyscan.uniform_sample(red, s_prime, False)
        b_sample_prime = pyscan.uniform_sample(blue, s_prime, False)

        for i in np.logspace(l_s, h_s, count):
            eps = i
            n = 1 / eps
            s = 1 / (2 * eps * eps)
            n = int(round(n) + .1)
            s = int(round(s) + .1)

            start_time = time.time()
            if sample_method == "block":
                f_sample = pyscan.block_sample
            elif sample_method == "even":
                f_sample = pyscan.even_sample
            elif sample_method == "uniform":
                f_sample = pyscan.uniform_sample

            m_sample = f_sample(red, s, False)
            b_sample = f_sample(blue, s, False)
            m_sample = pyscan.to_weighted(m_sample)
            b_sample = pyscan.to_weighted(b_sample)

            if two_level_sample:
                net_set1 = pyscan.my_sample(m_sample, n)
                net_set2 = pyscan.my_sample(b_sample, n)
                if ham_sample:
                    s = int(1 / (2 * eps ** (4.0 / 3)) * math.log(1/eps)**(2/3.0))

                    m_sample = pyscan.ham_tree_sample(m_sample, s)
                    b_sample = pyscan.ham_tree_sample(b_sample, s)
            else:
                net_set1 = [pyscan.Point(p[0], p[1], p[2]) for p in m_sample]
                net_set2 = [pyscan.Point(p[0], p[1], p[2]) for p in b_sample]
                n = s

            net_set1 = [pyscan.Point(p[0], p[1], p[2]) for p in net_set1]
            net_set2 = [pyscan.Point(p[0], p[1], p[2]) for p in net_set2]
            net_set = net_set1 + net_set2


            if region_name == "halfplane":
                reg, mx = pyscan.max_halfplane(net_set, m_sample, b_sample, disc)
            elif region_name == "disk":
                reg, mx = pyscan.max_disk(net_set, m_sample, b_sample, disc)
            elif region_name == "rectangle":
                grid = pyscan.Grid(n, m_sample, b_sample)
                s1 = pyscan.max_subgrid_linear(grid, -1.0, 1.0)
                s2 = pyscan.max_subgrid_linear(grid, 1.0, -1.0)
                if s1.fValue() > s2.fValue():
                    reg = grid.toRectangle(s1)
                    mx = s1.fValue()
                else:
                    reg = grid.toRectangle(s2)
                    mx = s2.fValue()

            end_time = time.time()


            actual_mx = pyscan.evaluate_range(reg, pyscan.to_weighted(m_sample_prime), pyscan.to_weighted(b_sample_prime), disc)

            row = { "disc": disc_name,
                   "region": region_name,
                   "n": n, "s": s, "r": r, "q": q,"p":p,
                   "time": end_time - start_time,
                   "m_disc_approx": mx,
                   "m_disc": actual_mx,
                   "sample_method": sample_method}
            writer.writerow(row)
            print(row)
            f.flush()
            if max_time is not None and end_time - start_time > max_time:
                return

if __name__ == "__main__":

    #trajectories = paths.read_geolife_files(100)
    trajectories = paths.read_dong_csv("/data/Dong_sets/Trajectory_Sets/samples/bjtaxi_samples_100k.tsv")
    print(len(trajectories))
    r = .0025
    q = .2
    p = .5
    eps_r = .001
    region_name = "halfplane"
    disc = utils.disc_to_func("disc")
    if region_name == "disk":
        red, blue, _, _ = pyscan.plant_partial_disk(trajectories, r, p, q, eps_r, disc)
    elif region_name == "halfplane":
        red, blue, _, _ = pyscan.plant_partial_halfplane(trajectories, r, p, q, eps_r, disc)
    elif region_name == "rectangle":
        red, blue, _, _ = pyscan.plant_partial_rectangle(trajectories, r, p, q, eps_r, disc)


    for two_level_sample, ham_sand in [(True, True), (False, False), (True, False)]:

        output_file = "partial_alg_progression_{}_{}_{}.csv".format(region_name, "2" if two_level_sample else "1", "ham" if ham_sand else "rand")
        testing_partial_framework(output_file, trajectories, -1, -4, 80, r=r, q=q, p=p,
                                region_name=region_name,
                                  two_level_sample=two_level_sample,
                                  ham_sample=ham_sand,
                                sample_method="uniform",
                                max_time=20,
                                planted_points=(red, blue))
