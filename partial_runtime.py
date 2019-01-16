import partial_testing
import paths
import utils
import pyscan

if __name__ == "__main__":

    #trajectories = paths.read_geolife_files(100)
    trajectories = paths.read_dong_csv("/data/Dong_sets/Trajectory_Sets/samples/osm_eu_sample_100k.tsv")
    r = .0025
    q = .2
    p = .5
    eps_r = .001
    disc = utils.disc_to_func("disc")
    red, blue, _, _ = pyscan.plant_partial_disk(trajectories, r, p, q, eps_r, disc)

    for region_name, two_level_sample, ham_sand in [("disk", True, False), ("halfplane", True, False), ("halfplane", True, True), ("rectangle", True, False)]:

        output_file = "partial_runtime_{}_{}_{}.csv".format(region_name, "2" if two_level_sample else "1", "ham" if ham_sand else "rand")
        partial_testing.testing_partial_framework(red, blue, output_file, -1, -4, 20, r=r, q=q, p=p,
                                  region_name=region_name,
                                  two_level_sample=two_level_sample,
                                  ham_sample=ham_sand,
                                  sample_method="even",
                                  max_time=None)
