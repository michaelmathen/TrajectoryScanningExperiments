import full_testing
import pyscan

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
        max_r = 1 / 100
        min_r = alpha
        c = 0


        for region, two_l_samp, sample_method, fast_disk in [("multiscale_disk", True, "grid_direc", True)]:
            if region == "multiscale_disk":
                red, blue = full_testing.generate_disk_sets(fname, r, p, q, min_r, max_r)
            full_testing.testing_full_framework(red, blue,
                                                "full_alpha_bc_error_{}_{}_{}.csv".format(c, fname, region), -1, -5, 40,
                                                vparam="alpha",
                                                eps=.005,
                                                region_name=region,
                                                sample_method=sample_method,
                                                alpha=alpha,
                                                two_level_sample=two_l_samp,
                                                fast_disk=fast_disk,
                                                min_disk_r=min_r,
                                                max_disk_r=max_r,
                                                max_time = 10000)
            full_testing.testing_full_framework(red, blue,
                                                "full_alpha_bc_error_{}_{}_{}.csv".format(c, fname, region), -1, -4, 40,
                                                vparam="eps",
                                                eps=.005,
                                                region_name=region,
                                                sample_method=sample_method,
                                                alpha=alpha,
                                                two_level_sample=two_l_samp,
                                                fast_disk=fast_disk,
                                                min_disk_r=min_r,
                                                max_disk_r=max_r,
                                                max_time = 10000)
            c += 1


