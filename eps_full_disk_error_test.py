import full_testing

if __name__ == "__main__":

    for fname in ["osm_eu_sample_100k_nw"]:

        r = .005
        p = .5
        q = .2
        #bjtax is 1/500 and 1/50
        #osm is 1/6000 and 1/300
        #alpha = 1/6000
        #max_r = 1/50
        alpha = 1 / 6000
        max_r = 1 / 300
        min_r = alpha
        c = 0

        red, blue = full_testing.generate_disk_sets(fname, r, p, q, min_r, max_r)
        for region, two_l_samp, sample_method, fast_disk in [("disk", True, "even", False), ("multiscale_disk_fixed", True, "even", False), ("multiscale_disk", True, "grid_direc", False), ("multiscale_disk", True, "grid_direc", True)]:
            full_testing.testing_full_framework(red, blue,
                                                "full_eps_disk_error_{}_{}_{}.csv".format(c, fname, region), -1, -4, 80,
                                                vparam="eps",
                                                region_name=region,
                                                sample_method=sample_method,
                                                alpha=alpha,
                                                two_level_sample=two_l_samp,
                                                fast_disk=fast_disk,
                                                min_disk_r=min_r,
                                                max_disk_r=max_r,
                                                max_time = 1000)
            c += 1


