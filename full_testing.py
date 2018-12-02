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
        eps_r=.01,
        r=.04,
        p=0.5,
        q=.2,
        alpha=.01,
        error_thresh=3,
        planted_points=None,
        max_disk_r=.1,
        min_disk_r=.05,
        max_w = .1,
        disc_name="disc",
        region_name="halfplane",
        sample_method="halfplane",
        input_size=10000,
        prefix=""):

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
    output_file = "{5}_{0}_{1}_{2}_{3:.2f}_{4:.2f}_full_discrepancy.csv".format(disc_name, region_name, sample_method, min_disk_r,max_disk_r, prefix)

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

            if planted_points is None:
                if "disk" in region_name or "rectangle" in region_name:
                    red, blue, _, actual_mx = pyscan.plant_full_disk(trajectories, r, p, q, disc)
                elif region_name == "halfplane":
                    red, blue, _, actual_mx = pyscan.plant_full_halfplane(trajectories, r, p, q, disc)
            else:
                red, blue = planted_points


            red_sample = pyscan.my_sample(red, s)
            blue_sample = pyscan.my_sample(blue, s)
            red_net = pyscan.my_sample(red, n)
            blue_net = pyscan.my_sample(blue, n)
            net = red_net + blue_net

            print("Running: {} {}".format(n, s))
            start_time = time.time()

            if sample_method == "halfplane":
                m_sample = [pyscan.halfplane_kernel([pyscan.Point(pt[0], pt[1], 1.0) for pt in traj], alpha) for traj in red_sample]
                b_sample = [pyscan.halfplane_kernel([pyscan.Point(pt[0], pt[1], 1.0) for pt in traj], alpha) for traj in blue_sample]
                pt_net = [pyscan.halfplane_kernel([pyscan.Point(pt[0], pt[1], 1.0) for pt in traj], alpha) for traj in net]
            elif sample_method == "grid":
                m_sample = [pyscan.grid_kernel(traj, alpha) for traj in red_sample]
                b_sample  = [pyscan.grid_kernel(traj, alpha) for traj in blue_sample]
                pt_net = [pyscan.grid_kernel(traj, alpha) for traj in net]
            elif sample_method == "grid_direc":
                chord_l = math.sqrt(4 * alpha * min_disk_r - 2 * alpha * alpha)
                m_sample = [pyscan.grid_direc_kernel(pyscan.dp_compress(traj, alpha), chord_l, alpha) for traj in red_sample]
                b_sample = [pyscan.grid_direc_kernel(pyscan.dp_compress(traj, alpha), chord_l, alpha) for traj in blue_sample]
                pt_net = [pyscan.grid_direc_kernel(pyscan.dp_compress(traj, alpha), chord_l, alpha) for traj in net]
            elif sample_method == "lifting":
                m_sample = [pyscan.lifting_kernel(traj, alpha) for traj in red_sample]
                b_sample = [pyscan.lifting_kernel(traj, alpha) for traj in blue_sample]
                pt_net = [pyscan.lifting_kernel(traj, alpha) for traj in net]
            elif sample_method == "dp":
                m_sample = [pyscan.dp_compress(traj, alpha) for traj in red_sample]
                b_sample = [pyscan.dp_compress(traj, alpha) for traj in blue_sample]
                pt_net = [pyscan.dp_compress(traj, alpha) for traj in net]
            else:
                print("Undefined sample method")
                return

            m_sample = list(pyscan.trajectories_to_labels(m_sample))
            b_sample = list(pyscan.trajectories_to_labels(b_sample))
            net_set = list(itertools.chain.from_iterable(pt_net))
            print("Points per trajectory {} {} {}".format(len(m_sample) / s, len(b_sample) / s, len(net_set) / n))
            print("Total Points {} {} {}".format(len(m_sample), len(b_sample), len(net_set)))

            if region_name == "halfplane":
                reg, mx = pyscan.max_halfplane_labeled(net_set, m_sample, b_sample, disc)
            elif region_name == "disk":
                reg, mx = pyscan.max_disk_labeled(net_set, m_sample, b_sample, disc)
            elif region_name == "small_disk":
                #net_set = list(pyscan.trajectories_to_labels(pt_net))
                mx = -1
                curr_disk_r = min_disk_r
                reg = None
                while curr_disk_r < max_disk_r:
                    new_reg, new_mx = pyscan.max_disk_scale_labeled_alt(net_set, m_sample, b_sample, -1, curr_disk_r, disc)
                    if new_mx > mx:
                        reg = new_reg
                        mx = new_mx
                    curr_disk_r *= 2

            elif region_name == "rectangle":
                reg, mx = pyscan.max_rect_labeled(2 * n, max_w, m_sample, b_sample, disc)
            else:
                print("Undefined region_name")
                return

            end_time = time.time()
            if "disk" in region_name:
                red_w = 0
                red_total = 0
                for traj in red:
                    if pyscan.Trajectory(traj).intersects_disk(reg):
                        red_w += 1
                    red_total += 1
                blue_w = 0
                blue_total = 0
                for traj in blue:
                    if pyscan.Trajectory(traj).intersects_disk(reg):
                        blue_w += 1
                    blue_total += 1
            elif "halfplane" in region_name:
                red_w = 0
                red_total = 0
                for traj in red:
                    if pyscan.Trajectory(traj).intersects_halfplane(reg):
                        red_w += 1
                    red_total += 1
                blue_w = 0
                blue_total = 0
                for traj in blue:
                    if pyscan.Trajectory(traj).intersects_halfplane(reg):
                        blue_w += 1
                    blue_total += 1
            else:
                return

            actual_mx = pyscan.evaluate(disc, red_w, red_total, blue_w, blue_total)
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
            f.flush()
            print(row)


if __name__ == "__main__":

    for fname in ["osm_eu_sample_10k_nw", "osm_eu_sample_100k_nw", "osm_eu_sample_1m_nw"]:

        trajectories = paths.read_dong_csv("/data/Dong_sets/Trajectory_Sets/samples/{}.tsv".format(fname))
        print(len(trajectories))

        r = .05
        p = .5
        q = .2
        #bjtax is 1/500 and 1/50
        #osm is 1/6000 and 1/300
        alpha = 1/6000
        max_r = 1/300
        min_r = alpha
        disc = utils.disc_to_func("disc")
        for region_name, sample_method in [("small_disk", "grid_direc")]:
            if "disk" in region_name or "rectangle" in region_name:
                count = 0
                while count < 10:
                    red, blue, reg ,actual_mx = pyscan.plant_full_disk(trajectories, r, p, q, disc)
                    print(min_r, reg.get_radius(), max_r)
                    if min_r < reg.get_radius() < max_r:
                        break
                    else:
                        count += 1
                if count == 10:
                    import sys
                    print("region is too small or big")
                    sys.exit()
            elif region_name == "halfplane":
                red, blue, _, actual_mx = pyscan.plant_full_halfplane(trajectories, r, p, q, disc)
            print(actual_mx)
            #for min_r in np.linspace(alpha, max_r, 1):
            testing_full_framework(trajectories, -1, -3, 10, r=r, q=q, p=p, region_name=region_name,
                           sample_method=sample_method,
                            alpha=alpha,
                           planted_points=(red,blue),
                           prefix=fname,
                           min_disk_r = min_r,
                           max_disk_r = max_r)
            print(actual_mx)
