import flux_testing
import paths
import utils
import pyscan
import file

if __name__ == "__main__":

    #trajectories = paths.read_geolife_files(100)
    trajectories = paths.read_dong_csv(file.PATH + "bjtaxi_samples_100k_nw.tsv")

    st_pts, end_pts = pyscan.trajectories_to_flux(trajectories)
    st_pts = [pyscan.WPoint(1.0, float(p[0]), float(p[1]), 1.0) for p in st_pts]
    end_pts = [pyscan.WPoint(1.0, float(p[0]), float(p[1]), 1.0) for p in end_pts]

    r = .05
    q = .2
    p = .5
    eps_r = .001
    disc = utils.disc_to_func("disc")
    #red, blue, _, _ = pyscan.plant_partial_disk(trajectories, r, p, q, eps_r, disc)

    for region_name in ["disk"]:#, "rectangle", "disk"]:
        if region_name == "disk":
            method = pyscan.plant_disk
        elif region_name == "rectangle":
            method = pyscan.plant_rectangle
        else:
            method = pyscan.plant_halfplane

        red, blue, reg = pyscan.paired_plant_region(st_pts, end_pts, r, q, method)
        if region_name == "halfplane":
            for ham_sand in [False, True]:
                output_file = "flux_error_{}_{}.csv".format(region_name, "ham" if ham_sand else "rand")
                flux_testing.testing_flux_framework(output_file, red, blue, -1, -5, 20,
                                                    region_name=region_name,
                                                    two_level_sample=True,
                                                    ham_sample=ham_sand,
                                                    max_time=1000)
        else:
            output_file = "flux_error_{}_{}.csv".format(region_name, "rand")
            flux_testing.testing_flux_framework(output_file, red, blue, -1, -5, 20,
                                              region_name=region_name,
                                              two_level_sample=True,
                                              ham_sample=False,
                                              max_time=1000)
