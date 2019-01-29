import partial_testing
import paths
import utils
import pyscan

import itertools
import matplotlib.pyplot as plt
import random
import numpy as np


def traj_length_distribution(trajectories):

    traj_length = [pyscan.Trajectory(traj).get_length() for traj in trajectories]
    plt.step(sorted(traj_length), np.arange(len(traj_length)))
    #plt.hist(traj_length, bins=50, weights=np.array(traj_length) / sum(traj_length), density=True)
    #plt.semilogy()
    plt.show()

if __name__ == "__main__":

    #trajectories = paths.read_geolife_files(100)
    trajectories = paths.read_dong_csv("/data/Dong_sets/Trajectory_Sets/samples/osm_eu_sample_100k_nw.tsv", filter_long=True)

    r = .005
    q = .2
    p = 0.5
    eps_r = .001
    disc = utils.disc_to_func("disc")

    for region_name in ["halfplane", "disk", "rectangle"]:
        if region_name == "disk":
            red, blue, reg, d = pyscan.plant_partial_disk(trajectories, r, p, q, eps_r, disc)
        elif region_name == "rectangle":
            red, blue, reg, d = pyscan.plant_partial_rectangle(trajectories, r, p, q, eps_r, disc)
        else:
            red, blue, reg, d = pyscan.plant_partial_halfplane(trajectories, r, p, q, eps_r, disc)
        print(d)

        # m_sample_prime = pyscan.uniform_sample(red, 100000, False)
        # print("here")
        # i = 0
        # while True:
        #     i += 1
        #     wpts = pyscan.to_weighted(m_sample_prime)
        #     #n_wpts = pyscan.VectorWP()
        #     #n_wpts.extend(p for p in wpts)
        #     actual_mx = pyscan.evaluate_range(reg, wpts, wpts, disc)

        #     if i > 5:
        #         i = 0
        #         print(i)
        # def plot_tuples(ax, pts, c):
        #     xs = []
        #     ys = []
        #     for pt in pts:
        #         xs.append(pt[0])
        #         ys.append(pt[1])
        #     ax.plot(xs, ys, color=c, marker='.')
        #
        #
        # ax = plt.subplot()
        # rtraj = random.sample(red, 1000)
        # for traj in [traj.get_pts() for traj in rtraj]:
        #     plot_tuples(ax, traj, "r")
        #
        # btraj = random.sample(blue, 1000)
        # for traj in [traj.get_pts() for traj in btraj]:
        #     plot_tuples(ax, traj, "b")
        #
        # ax.set_xlim(0, 1.0)
        # ax.set_ylim(0, 1.0)
        # plt.tight_layout()
        # plt.axis('off')
        # plt.show()

        for approx in ["even", "uniform", "block"]:
            if region_name == "halfplane":
                for ham_sand in [True, False]:
                    output_file = "partial_error_{}_{}_{}.csv".format(region_name, "ham" if ham_sand else "rand", approx)
                    partial_testing.testing_partial_framework(red, blue, output_file, -1, -4.5, 30, r=r, q=q, p=p,
                                              region_name=region_name,
                                              two_level_sample=True,
                                              ham_sample=ham_sand,
                                              sample_method=approx,
                                              max_time=5)
            else:
                ham_sand = False
                output_file = "partial_error_{}_{}_{}.csv".format(region_name, "ham" if ham_sand else "rand", approx)
                partial_testing.testing_partial_framework(red, blue, output_file, -1, -4.5, 30, r=r, q=q, p=p,
                                                          region_name=region_name,
                                                          two_level_sample=True,
                                                          ham_sample=ham_sand,
                                                          sample_method=approx,
                                                          max_time=5)