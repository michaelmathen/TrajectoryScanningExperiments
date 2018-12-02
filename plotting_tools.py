import sys
import numpy
import csv
from matplotlib import pyplot
from scipy.ndimage.filters import gaussian_filter1d
from scipy.interpolate import interp1d
import csv

def error_plot(x, y, n, k, x_name, y_name, alpha=.1, color="red", marker="-", ax=None, plot_points=False):
    mx = max(*x)
    mn = min(*x)
    my = max(*y)
    mny = min(*y)
    poly = numpy.poly1d(numpy.polyfit(x, y, n))
    xp = numpy.linspace(mn, mx, 25)
    if ax is None:
        fig, ax = pyplot.subplots(1, 1, sharex=True)
    
    error = numpy.poly1d(numpy.polyfit(x, numpy.abs(y - poly(x)), k))
    pyplot.ylim([mny, my])
    pyplot.xlim([mn, mx])
    pyplot.xlabel(x_name)
    pyplot.ylabel(y_name)
    ax.fill_between(x, poly(x), error(x) + poly(x), alpha=alpha, color=color)
    ax.fill_between(q, poly(x) - error(x), poly(x), alpha=alpha, color=color)
    ax.plot(x, poly(x), '-')
    ax.plot(x, error(x) + poly(x), m)
    ax.plot(x, poly(x) - error(x), m)
    if plot_points:
        ax.plot(x, y, '.')
    


def K(u):
    return 1 / numpy.sqrt(2 * numpy.pi) * numpy.exp(-u ** 2 / 2)


def Kh(x, y, h):
    return K((x - y) / h)


def makeInterp(xs, ys, h):

    def f(x):
        num = 0.0
        den = 0.0
        for xi, yi in zip(xs, ys):
            kh = Kh(x, xi, h)
            num += kh * yi
            den += kh
        return num / den
    return f


def kernel_error_plot(ax, x, y, x_name, y_name, sig1=.005, color="r", marker='-', error_bars=False, alpha=.05, label=None):
    x = numpy.asarray(x)
    y = numpy.asarray(y)
    mx = max(*x)
    mn = min(*x)
    xp = numpy.linspace(mn, mx, 25)
    poly = makeInterp(x, y, sig1)
    error = makeInterp(x, numpy.abs(poly(x) - y), sig1)
    pyplot.xlim([mn, mx])
    pyplot.xlabel(x_name)
    pyplot.ylabel(y_name)
    if error_bars:
        ax.fill_between(xp, poly(xp) - error(xp), error(xp) + poly(xp), alpha=alpha, color=color)
    #pyplot.locator_params(nbins=4)
    ax.plot(xp, poly(xp), marker=marker, color=color, label=label)


def plot_interp(ax, x, y, sig=60, logsc=False):
    if not x:
        return None
    i_kind = "slinear"
    mx = max(*x)
    mn = min(*x)
    y_f = gaussian_filter1d(y, sig)
    if logsc:
        poly = interp1d(numpy.log(x), y_f, kind=i_kind)
        xp = numpy.linspace(numpy.log(mn), numpy.log(mx), 25)
        return ax.plot(numpy.exp(xp), poly(xp))
    else:
        poly = interp1d(x, y_f, kind=i_kind)
        xp = numpy.linspace(mn, mx, 25)
        return ax.plot(xp, poly(xp), '-')


def read_csv(fname, *args):
    var_outs = [[] for i in args]
    for name in args:
       with open(fname) as file:
           reader_obj = csv.DictReader(file)
           for row in reader_obj:
               for i, arg in enumerate(args):
                  var_outs[i].append(float(row[arg]))
    return var_outs
