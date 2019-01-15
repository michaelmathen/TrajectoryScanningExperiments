import numpy as np
import time
import csv
import pyscan
import paths
import utils
import math


def testing_flux_framework(
        output_file,
        red, blue,
        l_s, h_s, count,
        r=.04,
        q=.2,
        region_name="disk",
        two_level_sample=True,
        ham_sample=False,
        max_time=None):



    fieldnames = ["disc", "region", "n", "s", "r", "q", "time",
                  "m_disc",
                  "m_disc_approx"]
    with open(output_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for i in np.logspace(l_s, h_s, count):

            eps = i
            n = 1 / eps
            s = 1 / (2 * eps * eps)
            n = int(round(n) + .1)
            s = int(round(s) + .1)
      
            disc = utils.disc_to_func("disc")
            start_time = time.time()
            
            m_sample = [pyscan.WPoint(1.0, p[0], p[1], 1.0) for p in pyscan.my_sample(red, s)]
            b_sample = [pyscan.WPoint(1.0, p[0], p[1], 1.0) for p in pyscan.my_sample(blue, s)]

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
            else:
                return
            end_time = time.time()

            st = time.time()
            actual_mx = pyscan.evaluate_range(reg, red, blue, disc)
            et = time.time()
            print("Time to evaluate region {}".format(et - st))

            row = { "disc": "disc",
                   "region":region_name,
                   "n":n, "s":s, "r": r, "q":q,
                   "time":end_time - start_time,
                   "m_disc_approx": mx,
                   "m_disc": actual_mx}
            writer.writerow(row)
            f.flush()
            print(row)
            if max_time is not None and end_time - start_time > max_time:
                return


if __name__ == "__main__":

    #trajectories = paths.read_geolife_files(10000)
    trajectories = paths.read_dong_csv("/data/Dong_sets/Trajectory_Sets/samples/osm_eu_sample_100k_nw.tsv")
    st_pts, end_pts = pyscan.trajectories_to_flux(trajectories)
    st_pts = [pyscan.WPoint(1.0, float(p[0]), float(p[1]), 1.0) for p in st_pts]
    end_pts = [pyscan.WPoint(1.0, float(p[0]), float(p[1]), 1.0) for p in end_pts]

    r=.01
    q=.8
    eps_r=.01


    for region_name in ["halfplane", "disk"]:
        scan_f = utils.range_to_func(region_name)
        print("Started planting")
        red, blue, region = pyscan.paired_plant_region(st_pts, end_pts, r, q, eps_r, scan_f)

        print("got here")
        testing_flux_framework(st_pts, end_pts, -1, -3, 40, r=r, q=q, region_name=region_name, planted_points=(red, blue))
