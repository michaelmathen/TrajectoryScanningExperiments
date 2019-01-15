import flux_testing
import paths
import utils
import pyscan

if __name__ == "__main__":

    #trajectories = paths.read_geolife_files(100)
    trajectories = paths.read_dong_csv("/data/Dong_sets/Trajectory_Sets/samples/bjtaxi_samples_100k.tsv")

    st_pts, end_pts = pyscan.trajectories_to_flux(trajectories)
    st_pts = [pyscan.WPoint(1.0, float(p[0]), float(p[1]), 1.0) for p in st_pts]
    end_pts = [pyscan.WPoint(1.0, float(p[0]), float(p[1]), 1.0) for p in end_pts]

    r = .0025
    q = .2
    p = .5
    eps_r = .001
    disc = utils.disc_to_func("disc")
    #red, blue, _, _ = pyscan.plant_partial_disk(trajectories, r, p, q, eps_r, disc)

    for region_name, two_level_sample, ham_sand in [ ("halfplane", True, False), ("disk", True, False),("halfplane", True, True), ("rectangle", True, False)]:


        # if region_name == "halfplane":
        #     scan = pyscan.max_halfplane
        # elif region_name == "disk":
        #     scan = pyscan.max_disk
        # elif region_name == "rectangle":
        #     def scan(n_s, m_s, b_s, disc):
        #         grid = pyscan.Grid(len(n_s), m_s, b_s)
        #         s1 = pyscan.max_subgrid_linear(grid, -1.0, 1.0)
        #         s2 = pyscan.max_subgrid_linear(grid, 1.0, -1.0)
        #         if s1.fValue() > s2.fValue():
        #             reg = grid.toRectangle(s1)
        #             mx = s1.fValue()
        #         else:
        #             reg = grid.toRectangle(s2)
        #             mx = s2.fValue()
        #         return reg, mx

        output_file = "flux_runtime_{}_{}_{}.csv".format(region_name, "2" if two_level_sample else "1", "ham" if ham_sand else "rand")
        flux_testing.testing_flux_framework(output_file, st_pts, end_pts, -1, -4, 80, r=r, q=q,
                                  region_name=region_name,
                                  two_level_sample=two_level_sample,
                                  ham_sample=ham_sand,
                                  max_time=100)
