import numpy as np
import random
import time
import csv
import paths
import pyscan
import line_discrepancy as ld
import math


def testing_flux_framework(l_s, h_s, count,
                            vparam="eps",
                            eps=.01,
                            c=1,
                            n=100,
                            s=5000,
                            r=.04,
                            q=.2,
                            scan_function="disc",
                            region="disk",
                            input_size=1000):


    output_file = "{}_{}_discrepancy.csv".format(scan_function, region)


    fieldnames = ["vparam", "input_size", "scan_func", "region", "n", "s", "r", "q", "time",  "m_disc"]
    with open(output_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for i in np.logspace(l_s, h_s, count):

            if vparam == "eps":
                n = 1 / i
                s = 1 / (i * i)
            elif vparam == "n":
                n = i
            elif vparam == "s":
                s = i
            elif vparam == "alpha":
                alpha = i
            elif vparam == "r":
                r = i
            elif vparam == "q":
                q = i

            n = int(round(n) + .1)
            s = int(round(s) + .1)

            trajectories = paths.read_geolife_files(alpha, input_size, include_end_points)


            print("red {}, blue {}, rl {}, bl {}".format(len(rp), len(bp), len(set(bl)), len(set(rl))))
            rlabel_unique = set(rl)
            blabel_unique = set(bl)
            r_label_sample = set(random.sample(rlabel_unique, s))
            b_label_sample = set(random.sample(blabel_unique, s))
            rp, rl = zip(*[(p, l) for p, l in zip(rp, rl) if l in r_label_sample])
            bp, bl = zip(*[(p, l) for p, l in zip(bp, bl) if l in b_label_sample])
            nr_label_sample = set(random.sample(rlabel_unique, n))
            nb_label_sample = set(random.sample(blabel_unique, n))
            nrp = [p for p, l in zip(rp, rl) if l in nr_label_sample]
            nbp = [p for p, l in zip(bp, bl) if l in nb_label_sample]
            net = nrp + nbp

            start_time = time.time()

            if scan_function == "disc":
                scan = pyscan.DISC
            elif scan_function == "kul":
                scan = pyscan.KULLDORF

            if region == "line":
                (reg, mx) = pyscan.max_halfplane_labels(net, list(rp), [1.0 for _ in rl],
                                                        list(rl), list(bp),
                                                        [1.0 for _ in bl] ,list(bl), scan)
            elif region == "disk":
                (reg, mx) = pyscan.max_disk_labels(net, list(rp), [1.0 for _ in rl], list(rl),
                                                   list(bp), [1.0 for _ in bl], list(bl), scan)
            elif region == "rect":
                pass

            #TODO measure the exact discrepancy on the rp, rl, bp, bl sets
            m, b = ld.evaluate_region(reg, rp, rl, bp, bl, region_t=region)
            end_time = time.time()
            actual_mx = pyscan.evaluate(scan, m, b)
            row = {"vparam": vparam,
                   "input_size":input_size,
                   "scan_func": scan_function,
                   "region":region,
                   "n":n, "s":s, "r": r, "q":q,
                   "time":end_time - start_time,
                   "m_disc_approx": mx,
                   "m_disc": actual_mx}
            writer.writerow(row)
            print(row)

def testing_labels_framework(l_s, h_s, count,
                            vparam="eps",
                            n=100,
                            s=5000,
                            r=.04,
                            q=.2,
                            alpha=.01,
                            include_end_points=False,
                            scan_function="disc",
                            alpha_simpl = True,
                            region="disk",
                            input_size=1000):


    output_file = "{}_{}_discrepancy.csv".format(scan_function, region)


    fieldnames = ["vparam", "input_size", "scan_func", "region", "n", "s", "r", "q", "time",  "m_disc", "m_disc_approx"]
    with open(output_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        st = time.time()
        trajectories = paths.read_geolife_files(input_size)
        print("time to sample trajectories {}".format(time.time() - st))

        for i in np.linspace(l_s, h_s, count):

            if vparam == "eps":
                n = 1.0 / i
                s = 1.0 / (i * i)
            elif vparam == "n":
                n = i
            elif vparam == "s":
                s = i
            elif vparam == "alpha":
                alpha = i
            elif vparam == "r":
                r = i
            elif vparam == "q":
                q = i
            print(i)
            n = min(int(round(n) + .1), len(trajectories))
            s = min(int(round(s) + .1), len(trajectories))
            print("n = {}, s= {}".format(n, s))



            net_sample = random.sample(trajectories, n)
            red_sample = random.sample(trajectories, s)
            blue_sample = random.sample(trajectories, s)

            # Convert to the c++ objects. Need to implement custom converters.
            net_sample = [[pyscan.Point(pt[0], pt[1], 1.0) for pt in traj] for traj in net_sample]
            red_sample = [[pyscan.Point(pt[0], pt[1], 1.0) for pt in traj] for traj in red_sample]
            blue_sample = [[pyscan.Point(pt[0], pt[1], 1.0) for pt in traj] for traj in blue_sample]

            # st = time.time()
            # rlabel_unique = set(rl)
            # blabel_unique = set(bl)
            # print(s)
            # r_label_sample = set(random.sample(rlabel_unique, min(s, len(rlabel_unique))))
            # b_label_sample = set(random.sample(blabel_unique, min(s, len(blabel_unique))))
            #
            #
            # rp, rl = zip(*[(p, l) for p, l in zip(rp, rl) if l in r_label_sample])
            # bp, bl = zip(*[(p, l) for p, l in zip(bp, bl) if l in b_label_sample])
            # nr_label_sample = set(random.sample(rlabel_unique, n))
            # nb_label_sample = set(random.sample(blabel_unique, n))
            # nrp = [p for p, l in zip(rp, rl) if l in nr_label_sample]
            # nbp = [p for p, l in zip(bp, bl) if l in nb_label_sample]
            # net = nrp + nbp
            # print("time to generate net and samples {}".format(time.time() - st))

            start_time = time.time()
            if scan_function == "disc":
                scan = pyscan.DISC
            elif scan_function == "kul":
                scan = pyscan.KULLDORF

            if region == "disk":
                scale = alpha
                while scale < 1 / 9.0:

                    (reg, mx) = pyscan.max_traj_disk(net_sample,
                                                            red_sample,
                                                            [1.0  for p in red_sample],
                                                            blue_sample,
                                                             [1.0 for p in blue_sample],
                                                             alpha,
                                                            scale,
                                                            scan)
                    scale *= 2
            elif region == "rect":
                pass
            end_time = time.time()

            #TODO measure the exact discrepancy on the rp, rl, bp, bl sets
            #m, b = ld.evaluate_region(reg, rp, rl, bp, bl, region_t=region)
            #actual_mx = pyscan.evaluate(scan, m, b)
            row = {"vparam": vparam,
                   "input_size":input_size,
                   "scan_func": scan_function,
                   "region":region,
                   "n":n, "s":s, "r": r, "q":q,
                   "time":end_time - start_time,
                   "m_disc_approx": mx,
                   "m_disc": 0}
            writer.writerow(row)
            print(row)



if __name__ == "__main__":
    testing_labels_framework(.3, .01, 30, alpha=.001)