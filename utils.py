import random
import pyscan



def range_to_lfunc(name):
    if name == "halfplane":
        return pyscan.max_halfplane_labeled
    elif name == "disk":
        return pyscan.max_disk_labeled
    elif name == "rectangle":
        pass


def disc_to_func(name):
    if name == "disc":
        return pyscan.DISC
    elif name == "kulldorff":
        return pyscan.KULLDORF
    else:
        return pyscan.RKULLDORF


