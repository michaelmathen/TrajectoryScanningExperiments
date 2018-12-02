from libpyscan import *
import random
import itertools

def to_weighted(points):
    return [WPoint(1.0, pt[0], pt[1], 1.0) for pt in points]

def trajectories_to_flux(trajectories):
    return [t[0] for t in trajectories], [t[-1] for t in trajectories]

def trajectories_to_labels(trajectories):
    return itertools.chain.from_iterable([[LPoint(label, 1.0, pt[0], pt[1], 1.0) for pt in traj]
            for label, traj in zip(range(len(trajectories)), trajectories)])

def my_sample(samp, count):
    return random.sample(samp, min(len(samp), int(count + .5)))

def evaluate_range(range, mp, bp, disc_f):
    """
    Evaluates this range to compute the total discrepancy.
    :param range:
    :param mp:
    :param bp:
    :param disc_f:
    :return:
    """
    total_mw = sum(p.get_weight() for p in mp)
    range_mw = sum(p.get_weight() for p in mp if range.contains(p))
    total_bw = sum(p.get_weight() for p in bp)
    range_bw = sum(p.get_weight() for p in bp if range.contains(p))
    return evaluate(disc_f, range_mw, total_mw, range_bw, total_bw)



def split_set(pts, rate):
    """
    Divides the point set into two random sets where one contains approximately rate * len(pts) number of
    points and the other has (1 - rate) * len(pts)
    :param pts:
    :param rate:
    :return:
    """
    red_set = my_sample(pts, len(pts) * rate)
    red_set_set = set(red_set)
    return red_set, [item for item in pts if item not in red_set_set]


"""
This plants a region where every trajectory:
Completely outside or inside of the region has an endpoint chosen at random.
Every trajectory with one endpoint inside the region has an endpoint chosen inside
with probability q (exactly q fraction have one endpoint in the region)

r controls how many points the region contains.

"""
def paired_plant_region(traj_start, traj_end, r, q, eps, scan_f):

    while True:
        net_size = int( 1 / min(r, eps) + 1)
        s_size = int(1 / min(eps, r) ** 2 + 1)

        net = my_sample(traj_start + traj_end, net_size)
        sample = my_sample(traj_start + traj_end, s_size)
        sample = [WPoint(1.0, p[0], p[1], 1.0) for p in sample]
        disc = size_region(r)
        reg, _ = scan_f(net, sample, [], disc)

        val = evaluate_range(reg, traj_start + traj_end, [], disc)
        if 1 - val < eps:
            break

    flux_region = []
    out_region = []
    for st_pt, end_pt in zip(traj_start, traj_end):
        if reg.contains(st_pt) and not reg.contains(end_pt):
            flux_region.append((st_pt, end_pt))
        elif reg.contains(end_pt) and not reg.contains(st_pt):
            flux_region.append((end_pt, st_pt))
        else:
            out_region.append((st_pt, end_pt))

    q_fraction, remainder = split_set(flux_region, q)
    remainder = [(ep, sp) for (sp, ep) in remainder]

    q_fraction_o, remainder_o = split_set(out_region, .5)
    remainder = [(ep, sp) for (sp, ep) in remainder]
    red, blue = zip(*(q_fraction + q_fraction_o + remainder + remainder_o))
    return red, blue, reg



"""
This takes a scanning function and two point sets and then computes a planted region that contains some fraction r
of the points with some tolerance.
This then computes a red and blue set of points based on the planted region.

Inside the region q fraction of points are red.
Outside the region p fraction of points are red
"""
def plant_region(points, r, p, q, eps, scan_f):

    while True:

        net_size = int(1 / min(r, eps) + 1)
        s_size = int(1 / min(r, eps)** 2 + 1)

        net = my_sample(points, net_size)
        sample = my_sample(points, s_size)
        disc = size_region(r)
        reg, _ = scan_f(net, sample, [], disc)

        val = evaluate_range(reg, points, [], disc)
        if 1 - val < eps:
            break

    in_region = []
    out_region = []
    for pt in pts:
        if reg.contains(pt):
            in_region.append(pt)
        else:
            out_region.append(pt)

    red_in, blue_in = split_set(in_region, q)

    red_out, blue_out = split_set(out_region, p)
    return red_in + red_out, blue_in + blue_out, reg


def plant_trajectory_rectangles(trajectories, r, p, q):
    """
    TODO
    Choose a point at random from a trajectory and then expand outward from there.
    :param trajectories this consists of lists of lists of points where the points are type Pyscan.Point
    :param r:
    :param p:
    :param q:
    :return:
    """
    pass




def plant_full_halfplane(trajectories, r, p, q, disc):
    """
    Choose a point at random from a trajectory and then expand outward from there.
    :param trajectories this consists of lists of lists of points where the points are type Pyscan.Point
    :param r:
    :param p:
    :param q:
    :return:
    """
    def min_distance(pts, direc):
        return min([direc[0] * pt[0] + direc[1] * pt[1] for pt in pts])
    rand_direc = (random.gauss(0,1), random.gauss(0, 1))
    trajectories = sorted(trajectories, key=lambda x: min_distance(x, rand_direc))

    pt = min(trajectories[int(r * len(trajectories))], key=lambda pt: rand_direc[0] * pt[0] + rand_direc[1] * pt[1])
    lc = pt[0] * rand_direc[0] + pt[1] * rand_direc[1]
    plant_region = Halfplane(Point(rand_direc[0], rand_direc[1], -lc))

    inside_plane = trajectories[:int(r * len(trajectories))]
    outside_plane = trajectories[int(r * len(trajectories)):]
    red_in, blue_in = split_set([tuple(traj) for traj in inside_plane], q)
    red_out, blue_out = split_set([tuple(traj) for traj in outside_plane], p)

    diff = evaluate(disc, len(red_in), len(red_in) + len(red_out), len(blue_in), len(blue_in) + len(blue_out))
    return red_in + red_out, blue_in + blue_out, plant_region, diff


def measure_halfplane_full(region, traj_red, traj_blue, disc):
    red_w = 0
    blue_w = 0
    for traj in traj_red:
        if Trajectory(traj).intersects_halfplane(region):
            red_w += 1
    for traj in traj_blue:
        if Trajectory(traj).intersects_halfplane(region):
            blue_w += 1
    return evaluate(disc, red_w, len(traj_red), blue_w, len(traj_blue))


def measure_disk_full(region, traj_red, traj_blue, disc):
    red_w = 0
    blue_w = 0
    for traj in traj_red:
        if Trajectory(traj).intersects_disk(region):
            red_w += 1
    for traj in traj_blue:
        if Trajectory(traj).intersects_disk(region):
            blue_w += 1
    return evaluate(disc, red_w, len(traj_red), blue_w, len(traj_blue))

def plant_partial_halfplane(trajectories, r, p, q, eps, disc):
    """
    Choose a point at random from a trajectory and then expand outward from there.
    :param trajectories this consists of lists of lists of points where the points are type Pyscan.Point
    :param r:
    :param p:
    :param q:
    :return:
    """
    def min_distance(pt, direc):
        return direc[0] * pt[0] + direc[1] * pt[1]
    rand_direc = (random.gauss(0,1), random.gauss(0, 1))

    trajectory_obj = [Trajectory(pts) for pts in trajectories]
    all_pts = uniform_sample(trajectory_obj, int(1 / eps ** 2 + 1), False)
    all_pts = sorted(all_pts, key=lambda x: min_distance(x, rand_direc))

    pt = all_pts[int(r * len(trajectories))]
    lc = pt.pdot(rand_direc)
    plant_region = Halfplane(Point(rand_direc[0], rand_direc[1], -lc))

    inside_disk = [traj for traj in trajectory_obj if traj.intersects_halfplane(plant_region)]
    outside_disk = [traj for traj in trajectory_obj if not traj.intersects_halfplane(plant_region)]

    red_in, blue_in = split_set(inside_disk, q)
    red_out, blue_out = split_set(outside_disk, p)
    diff = evaluate(disc, len(red_in), len(red_in) + len(red_out), len(blue_in), len(blue_in) + len(blue_out))
    return red_in + red_out, blue_in + blue_out, plant_region, diff


def plant_full_disk(trajectories, r, p, q, disc):
    """
    Choose a point at random from a trajectory and then expand outward from there.
    :param trajectories this consists of lists of lists of points where the points are type Pyscan.Point
    :param r:
    :param p:
    :param q:
    :return:
    """
    origin = random.choice(list(itertools.chain.from_iterable(trajectories)))
    trajectory_obj = [Trajectory(pts) for pts in trajectories]

    trajectories = sorted(trajectory_obj, key=lambda el: el.point_dist(origin))

    inside_disk = trajectories[:int(r * len(trajectories))]
    outside_disk = trajectories[int(r * len(trajectories)):]

    max_disk = [Disk(origin[0], origin[1], traj.point_dist(origin)) for traj in inside_disk]
    print([dsk.get_radius() for dsk in max_disk])
    max_disk = max_disk[-1]

    red_in, blue_in = split_set([tuple(traj.get_pts()) for traj in inside_disk], q)
    red_out, blue_out = split_set([tuple(traj.get_pts()) for traj in outside_disk], p)
    diff = evaluate(disc, len(red_in), len(red_in) + len(red_out), len(blue_in), len(blue_in) + len(blue_out))
    return red_in + red_out, blue_in + blue_out, max_disk, diff


def plant_partial_disk(trajectories, r, p, q, eps, disc):
    """
    Choose a point at random from a trajectory and then expand outward from there.
    :param trajectories this consists of lists of lists of points where the points are type Pyscan.Point
    :param r:
    :param p:
    :param q:
    :return:
    """
    #Compute fraction of mass inside disk
    #For each segment inside disk add mass.
    # For each segment partly inside disk compute overlap and add to mass.
    #We can do bisection on this amount since it is a monotonic function to compute within eps fraction.

    origin_tuple = random.choice(list(itertools.chain.from_iterable(trajectories)))
    origin = Point(origin_tuple[0], origin_tuple[1], 1.0)
    trajectory_obj = [Trajectory(pts) for pts in trajectories]
    all_pts = uniform_sample(trajectory_obj, int(1 / eps ** 2 + 1), False)
    trajectories = sorted(all_pts, key=lambda x: x.dist(origin))

    disk_boundary = trajectories[int(r * len(trajectories))]
    max_disk = Disk(origin[0], origin[1], origin.dist(disk_boundary))
    inside_disk = [traj for traj in trajectory_obj if traj.intersects_disk(max_disk)]
    outside_disk = [traj for traj in trajectory_obj if not traj.intersects_disk(max_disk)]

    red_in, blue_in = split_set(inside_disk, q)
    red_out, blue_out = split_set(outside_disk, p)
    diff = evaluate(disc, len(red_in), len(red_in) + len(red_out), len(blue_in), len(blue_in) + len(blue_out))
    return red_in + red_out, blue_in + blue_out, max_disk, diff

"""
Create a distribution of the null distribution to measure the significance of the region.
"""
def distribution(points, p, scan_f, n, s, disc=DISC):
    while True:
        red, blue = split_set(points, p)

        net_set = my_sample(red, min(len(red), n)) + my_sample(blue, min(len(blue), n))
        m_sample = my_sample(red, min(len(red), s), red)
        b_sample = my_sample(blue, min(len(blue), s), blue)
        reg, val = scan_f(net_set, m_sample, b_sample, disc)

        yield val


def null_cdf(observations):
    """
    Generates a cdf object using a certain number of observations.
    :param observations
    :return:
    """
    values = sorted(observations)
    prob = [x / float(len(observations)) for x in range(len(observations))]
    return values, prob



