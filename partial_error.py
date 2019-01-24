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

    for region_name in ["disk", "halfplane", "rectangle"]:
        if region_name == "disk":
            red, blue, _, _ = pyscan.plant_partial_disk(trajectories, r, p, q, eps_r, disc)
        elif region_name == "rectangle":
            red, blue, _, _ = pyscan.plant_partial_rectangle(trajectories, r, p, q, eps_r, disc)
        else:
            red, blue, _, _ = pyscan.plant_partial_halfplane(trajectories, r, p, q, eps_r, disc)
        for approx in ["even", "uniform", "block"]:
            if region_name == "halfplane":
                for ham_sand in [True, False]:
                    output_file = "partial_error_{}_{}_approx.csv".format(region_name, "2", "ham" if ham_sand else "rand", approx)
                    partial_testing.testing_partial_framework(red, blue, output_file, -1, -4, 20, r=r, q=q, p=p,
                                              region_name=region_name,
                                              two_level_sample=True,
                                              ham_sample=ham_sand,
                                              sample_method=approx,
                                              max_time=1000)
            else:
                ham_sand = False
                output_file = "partial_error_{}_{}_approx.csv".format(region_name, "2", "ham" if ham_sand else "rand",
                                                                      approx)
                partial_testing.testing_partial_framework(red, blue, output_file, -1, -4, 20, r=r, q=q, p=p,
                                                          region_name=region_name,
                                                          two_level_sample=True,
                                                          ham_sample=ham_sand,
                                                          sample_method=approx,
                                                          max_time=1000)