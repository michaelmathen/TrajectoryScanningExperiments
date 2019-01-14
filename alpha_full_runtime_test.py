import full_testing

if __name__ == "__main__":

    for fname in ["osm_eu_sample_100k_nw"]:

        r = .05
        p = .5
        q = .2
        #bjtax is 1/500 and 1/50
        #osm is 1/6000 and 1/300
        #alpha = 1/6000
        #max_r = 1/50
        alpha = 1 / 500
        max_r = float("inf")
        min_r = alpha
        c = 0

        red, blue = full_testing.generate_disk_sets(fname, r, p, q, min_r, max_r)
        for region, two_l_samp, sample_method, fast_disk in [("disk", True, "even", False), ("halfspace", True, "even", False), ("rectangle", True, "even", False)]:
            full_testing.testing_full_framework(red, blue,
                                                "full_runtime_{}_{}_{}.csv".format(c, fname, region), -.1, -5, 80,
                                                eps=.025,
                                                vparam="alpha",
                                                region_name=region,
                                                sample_method=sample_method,
                                                alpha=alpha,
                                                two_level_sample=two_l_samp,
                                                fast_disk=fast_disk,
                                                min_disk_r=min_r,
                                                max_disk_r=max_r,
                                                max_time = 1000)
            c += 1


