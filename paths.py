import os
import uuid
import gpxpy
import gpxpy.gpx
import csv
import random
import pyscan
import math
import itertools

import matplotlib.pyplot as plt


gpx_files = "/data/gpx-planet-2013-04-09/identifiable"

geolife = "/data/Trajectory_Sets/Geolife Trajectories 1.3/Data"
output = "/data/Trajectory_Sets/gpx_traces_near_ny"

def all_files(root):
    """
    Recursively gets all the subfiles.
    :param root:
    :return:
    """
    for root, directories, filenames in os.walk(root):
        for filename in filenames:
            yield os.path.join(root,filename)


def only_gpx(root):
    """
    Recursively gets all the subfiles, but only if they end with gpx.
    :param root:
    :return:
    """
    for f in all_files(root):
        if f.endswith(".gpx"):
            yield f

def only_plt(root):
    """
    Recursively gets all the subfiles, but only if they end with gpx.
    :param root:
    :return:
    """
    for f in all_files(root):
        if f.endswith(".plt"):
            yield f


def all_traces(root):
    """
    Keeps retrieving gps traces until all have been exhausted.
    :param root:
    :return:
    """
    for file_name in only_gpx(root):
        with open(file_name, 'r') as gpx_file:
            gpx = gpxpy.parse( gpx_file )
            for track in gpx.tracks:
                traj = []
                for segment in track.segments:
                    traj.extend([(point.latitude, point.longitude, point.time) for point in segment.points])
                yield traj




def inside_region(lat_lower, lat_upper, lon_lower, lon_upper, trace):
    for pt in trace:
        if lat_lower <= pt[0] <= lat_upper and lon_lower <= pt[1] <= lon_upper:
            continue
        else:
            return False
    return True


def traces_in_region(lat_lower, lat_upper, lon_lower, lon_upper, root):
    """
    Keeps all traces that are contained completely inside of this region.
    :param lat_lower:
    :param lat_upper:
    :param lon_lower:
    :param lon_upper:
    :param root:
    :return:
    """
    ix = 0
    for trace in all_traces(root):
        if inside_region(lat_lower, lat_upper, lon_lower, lon_upper, trace):
            print("{0} number current".format(ix))
            ix += 1
            yield trace


def write_trace_file(lat_lower, lat_upper, lon_lower, lon_upper, output_dir):

    for trace in traces_in_region(lat_lower, lat_upper, lon_lower, lon_upper, gpx_files):
        uid = uuid.uuid4()
        with open(output_dir + "/" + str(uid) + ".csv", 'w') as f:
            writer = csv.DictWriter(f, fieldnames=["latitude", "longitude", "time"])
            for pt in trace:
                writer.writerow({"latitude":pt[0], "longitude":pt[1], "time":str(pt[2])})


def center(lat_lower, lat_upper):
    return (lat_lower + lat_upper) / 2


LAT_L = 33.0
LAT_U = 42.0
LON_U = -71.0
LON_L = -90.0
RADIUS = 6.37e6


def equirectangular_projection(latitude, longitude):
    lamb = math.radians(latitude - LAT_L)
    gamma = math.radians(longitude - LON_L)
    x = RADIUS * lamb * math.cos(math.radians(center(LAT_L, LAT_U)))
    y = RADIUS * gamma
    return x, y


def equirectangular_projection_box(latitude, longitude, lat_l, lat_u, lon_l):
    lamb = math.radians(latitude - LAT_L)
    gamma = math.radians(longitude - LON_L)
    x = RADIUS * lamb * math.cos(math.radians(center(LAT_L, LAT_U)))
    y = RADIUS * gamma
    return x, y


def read_trace_file(trace_file_name):
    trace = []
    with open(trace_file_name, 'r') as f:
        reader = csv.DictReader(f, fieldnames=["latitude", "longitude", "time"])
        ignore = False
        for seg_pt in reader:
            if not ignore:
                x, y = equirectangular_projection(float(seg_pt["latitude"]), float(seg_pt["longitude"]))
                trace.append((x, y))
            ignore = not ignore
    return trace


def sample_trace_files(num_trace_files):
    fname_sample = random.sample([fname for fname in all_files(output)], num_trace_files)
    return [read_trace_file(trace_name) for trace_name in fname_sample]


def dist(pt1, pt2):
    difx = pt1[0] - pt2[0]
    dify = pt1[1] - pt2[1]
    return math.sqrt(difx ** 2 + dify ** 2)


def normalize(pt, mxx, mnx, mxy, mny):
    return (pt[0] - mnx[0]) / (mxx[0] - mnx[0]), (pt[1] - mny[1]) / (mxy[1] - mny[1])


def normalize_all_projection(traces):
    mxx = max(list(itertools.chain(*[trace for trace in traces])), key=lambda x: x[0])
    mnx = min(list(itertools.chain(*[trace for trace in traces])), key=lambda x: x[0])
    mxy = max(list(itertools.chain(*[trace for trace in traces])), key=lambda x: x[1])
    mny = min(list(itertools.chain(*[trace for trace in traces])), key=lambda x: x[1])

    proj_traces = []
    for trace in traces:
        proj_trace = []
        for pt in trace:
            x, y = equirectangular_projection_box(pt[0], pt[1], mnx, mxx, mny)
            proj_trace.append(normalize((x, y), mxx, mnx, mxy, mny))
        proj_traces.append(proj_trace)
    mxx = max(list(itertools.chain(*[trace for trace in proj_traces])), key=lambda x: x[0])
    mnx = min(list(itertools.chain(*[trace for trace in proj_traces])), key=lambda x: x[0])
    mxy = max(list(itertools.chain(*[trace for trace in proj_traces])), key=lambda x: x[1])
    mny = min(list(itertools.chain(*[trace for trace in proj_traces])), key=lambda x: x[1])

    norm_traces = []
    for trace in proj_traces:
        norm_trace = []
        for pt in trace:
            x, y = pt
            norm_trace.append(normalize(pyscan.Point(x, y, 1.0), mxx, mnx, mxy, mny))
        norm_traces.append(norm_trace)
    return norm_traces


def normalize_all(traces):
    mxx = max(list(itertools.chain(*[trace for trace in traces])), key=lambda x: x[0])
    mnx = min(list(itertools.chain(*[trace for trace in traces])), key=lambda x: x[0])
    mxy = max(list(itertools.chain(*[trace for trace in traces])), key=lambda x: x[1])
    mny = min(list(itertools.chain(*[trace for trace in traces])), key=lambda x: x[1])

    norm_traces = []
    for trace in traces:
        norm_trace = []
        for pt in trace:
            norm_trace.append(normalize(pt, mxx, mnx, mxy, mny))
        norm_traces.append(norm_trace)
    return norm_traces


def alpha_simplification(alpha, trace, include_end_points = False):
    new_trace = []
    curr_alpha = 0
    for lpt, rpt in zip(trace, trace[1:]):
        curr_pt_seg = dist(lpt, rpt)
        while curr_pt_seg + curr_alpha > alpha:
            curr_pt_seg -= alpha - curr_alpha
            a = curr_pt_seg / dist(lpt, rpt)
            new_trace.append(pyscan.Point(lpt[0] * a + rpt[0] * (1 - a), lpt[1] * a + rpt[1] * (1 - a), 1.0))
            curr_alpha = 0

        curr_alpha += curr_pt_seg
        if include_end_points:
            new_trace.append(pyscan.Point(lpt[0], lpt[1], 1.0))
    if include_end_points:
        new_trace.append(pyscan.Point(trace[-1][0], trace[-1][1], 1.0))

    return new_trace


# TODO
# def alpha_simplification_kernel2d(alpha, trace):
#     new_trace = []
#     curr_alpha = 0
#     pyscan.approximate_hull()
# 
#     return new_trace


def plot_trace(ax, trace):
    if trace:
        x, y = zip(*trace)
        ax.plot(x, y, '.')



def read_geolife_files(count):
    traj_set = list(only_plt('/data/Trajectory_Sets/Geolife Trajectories 1.3'))
    #print(len(traj_set))
    trajectory_files = random.sample(traj_set, count)

    all_traces = []
    for fname in trajectory_files:
        with open(fname, 'r') as f:
            ix = 0
            for line in f:
                ix += 1
                if ix >= 6:
                    break
            reader = csv.reader(f)
            trace = []
            for row in reader:
                trace.append(pyscan.Point(float(row[0]), float(row[1]), 1))
            all_traces.append(trace)
    normed_traces = normalize_all_projection(all_traces)
    return normed_traces


# for waypoint in gpx.waypoints:
#     print 'waypoint {0} -> ({1},{2})'.format( waypoint.name, waypoint.latitude, waypoint.longitude )
#
# for route in gpx.routes:
#     print 'Route:'
#     for point in route:
#         print 'Point at ({0},{1}) -> {2}'.format( point.latitude, point.longitude, point.elevation )
#
# # There are more utility methods and functions...
#
# # You can manipulate/add/remove tracks, segments, points, waypoints and routes and
# # get the GPX XML file from the resulting object:
#
# print 'GPX:', gpx.to_xml()


def norm_sample(count):
    return list(normalize_all(sample_trace_files(count)))


def traj_to_labels(traj):
    ix = 0
    pts = []
    labels = []
    for trace in traj:
        pts.extend([pyscan.Point(pt[0], pt[1], 1) for pt in trace])
        labels.extend([ix for _ in trace])
        ix += 1
    return pts, labels

print(sum(1 for f in all_files("/data/Trajectory_Sets/cabspottingdata")))


if __name__ == "__main__":

    f, ax = plt.subplots()#
    for trace in normalize_all(sample_trace_files(100)):
        plot_trace(ax, trace)
    plt.show()

