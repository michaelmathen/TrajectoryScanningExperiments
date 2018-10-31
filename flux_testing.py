import numpy as np
import random
import time
import csv
import pyscan
import paths
import utils
import matplotlib.pyplot as plt



def testing_flux_framework(
        st_pts, end_pts,
        l_s, h_s, count,
        vparam="eps",
        eps=.01,
        eps_r=.01,
        r=.04,
        q=.2,
        disc_name="disc",
        region_name="disk",
        input_size=10000):


    output_file = "{}_{}_flux_discrepancy.csv".format(disc_name, region_name)


    
    fieldnames = ["vparam", "input_size", "disc", "region", "n", "s", "r", "q", "time",
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
            print(eps)
            n = 1 / eps
            s = 1 / (2 * eps * eps)
            n = int(round(n) + .1)
            s = int(round(s) + .1)
      
            disc = utils.disc_to_func(disc_name)
            scan = utils.range_to_func(region_name)
            st = time.time()
            red, blue = pyscan.paired_plant_region(st_pts, end_pts, r, q, eps_r, scan)

            plt.scatter([x[0] for x in red], [x[1] for x in red], c="r")
            plt.scatter([x[0] for x in blue], [x[1] for x in blue], c="b")
            plt.show()

            et = time.time()
            print("Time to plant region {}".format(et - st))
            print(len(red), len(blue))

            start_time = time.time()
            
            net_set = pyscan.my_sample(red, n) + pyscan.my_sample(blue, n)
            m_sample = [pyscan.WPoint(1.0, p[0], p[1], 1.0) for p in pyscan.my_sample(red, s)]
            b_sample = [pyscan.WPoint(1.0, p[0], p[1], 1.0) for p in pyscan.my_sample(blue, s)]
            print(len(net_set), len(m_sample), len(b_sample))

            reg, mx = scan(net_set, m_sample, b_sample, disc)
            end_time = time.time()

            st = time.time()
            actual_mx = pyscan.evaluate_range(reg, red, blue, disc)
            et = time.time()
            print("Time to evaluate region {}".format(et - st))

            row = {"vparam": vparam,
                   "input_size":input_size,
                   "disc": disc_name,
                   "region":region_name,
                   "n":n, "s":s, "r": r, "q":q,
                   "time":end_time - start_time,
                   "m_disc_approx": mx,
                   "m_disc": actual_mx}
            writer.writerow(row)
            print(row)

if __name__ == "__main__":

    trajectories = paths.read_geolife_files(1000)
    st_pts, end_pts = pyscan.trajectories_to_flux(trajectories)
    print (st_pts, end_pts)
    st_pts = [pyscan.WPoint(1.0, float(p[0]), float(p[1]), 1.0) for p in st_pts]
    end_pts = [pyscan.WPoint(1.0, float(p[0]), float(p[1]), 1.0) for p in end_pts]
    print(st_pts)
    for region in ["disk"]:
        testing_flux_framework(st_pts, end_pts, -1, -3, 10, r=.1, q=0.0, region_name=region)
