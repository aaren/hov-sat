from __future__ import division

import sys
from glob import glob
import argparse

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
import Image

from mpl_toolkits.basemap import Basemap


def file_list(path='./', ext=''):
    """Create a list of files from a given path and extension.
    Default is to use all files in current directory.
    """
    files = sorted(glob(path + '*' + ext))
    return files


def pixels_from_latlon(X, Y, image, basemap, inverse=False):
    """Calculate equivalent pixels from given latitude and
    longitude for a given image representing a given basemap
    instance.

    Arguments:
        X - Longitude
        Y - Latitude
        image - A PIL image instance OR an image filename
        basemap - a Basemap instance representing the same region
                  as the image
        inverse - setting to True will compute the equivalent
                  lon / lat from given pixel coordinates.

    Returns:
        (px, py) - the image coordinates (i.e. measured from the upper-
                   left) that match the given latitude and longitude.
    """
    if type(image) is str:
        im = Image.open(image)
    elif hasattr(image, 'size'):
        im = image
    else:
        exit("The image is not a string or an image instance. Try again")

    m = basemap
    if inverse:
        # calculate lon / lat from given pixel coords
        x = X / im.size[0] * (m.xmax - m.xmin) + m.xmin
        y = (im.size[1] - Y) / im.size[1] * (m.ymax - m.ymin) + m.ymin
        lon, lat = m(x, y, inverse=True)
        return round(lon, 4), round(lat, 4)

    elif not inverse:
        x, y = m(X, Y)
        px = int(x / (m.xmax - m.xmin) * im.size[0])
        # y has to be changed to get image coords
        py = int(im.size[1] - y / (m.ymax - m.xmin) * im.size[1])
        return px, py


def extract_transect(im, coords, res=None):
    """
    From a given image array, return the values of the pixels on
    a given transect.

    See Stackoverflow_ for this exact question

    .. _Stackoverflow: http://stackoverflow.com/questions/7878398/how-to-extract-an-arbitrary-line-of-values-from-a-numpy-array

    Inputs
    im - an image array
            coords - (x0, y0, x1, y1). start and end coordinates
                      of the transect (pixels).
            res - if specified, transect will be created
                  using cubic interpolation. Otherwise
                  uses nearest neighbour (default).

    Returns: a tuple (x, y, z) of arrays
             x - x coords
             y - y coords
             z - values extracted from image
    """
    x0, y0, x1, y1 = coords
    if not res:
        # nearest neighbour
        length = int(np.hypot(x1 - x0, y1 - y0))
        x, y = np.linspace(x0, x1, length), np.linspace(y0, y1, length)
        # array indexing is transpose of image indexing!
        # TODO: check this
        zi = im[y.astype(int), x.astype(int)]

    elif res:
        x, y = np.linspace(x0, x1, res), np.linspace(y0, y1, res)
        # array indexing is transpose of image indexing!
        # TODO: check this
        coords = np.vstack((y, x))
        # cubic interpolation
        if im.ndim == 2:
            zi = ndimage.map_coordinates(im, coords)
        # cope with rgb images
        elif im.ndim > 2:
            zr = ndimage.map_coordinates(im[:, :, 0], coords)
            zg = ndimage.map_coordinates(im[:, :, 1], coords)
            zb = ndimage.map_coordinates(im[:, :, 2], coords)
            zi = np.dstack((zr, zg, zb)).squeeze()
        else:
            exit("unknown image dimensions")
    transect = (x, y, zi)
    return transect


def transect(image, transect_coords, image_coords,
             units='latlon', res=None):
    # Read in image and fit to lat lon grid
    # Return 2d array for image pixel intensities
    # and 2d array of coords
    print "Extracting %s \r" % image,
    sys.stdout.flush()
    im = Image.open(image)
    ima = np.array(im)

    lower_left, upper_right = image_coords
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    m = Basemap(projection='geos', lon_0=0, resolution=None,
                llcrnrlon=lower_left[0], llcrnrlat=lower_left[1],
                urcrnrlon=upper_right[0], urcrnrlat=upper_right[1])

    x0, y0, x1, y1 = transect_coords
    if units == 'latlon':
        X0, Y0 = pixels_from_latlon(x0, y0, im, m)
        X1, Y1 = pixels_from_latlon(x1, y1, im, m)
        coords = (X0, Y0, X1, Y1)
        X, Y, Z = extract_transect(ima, coords, res)
        # covert X, Y into lat lon
        XY = zip(X, Y)
        XYll = [pixels_from_latlon(x, y, im, m, inverse=True) for x, y in XY]
        XY_ll = zip(*XYll)
        X_ll = np.array(XY_ll[0])
        Y_ll = np.array(XY_ll[1])
        transect = X_ll, Y_ll, Z

        # imshow seems to expect what imread would give you
        # m.imshow(ima, origin='upper')
        # m.drawcoastlines()
        # m.drawcountries()
        # m.drawparallels(range(-40, 40, 10))
        # m.drawmeridians(range(-40, 40, 10))
        # X_m, Y_m = m(X_ll, Y_ll)
        # m.plot(X_m, Y_m, 'g')
        # plt.savefig('out.png')

    elif units == 'pixels':
        coords = transect_coords
        pixels_from_latlon(x0, y0, im, m, inverse=True)
        pixels_from_latlon(x1, y1, im, m, inverse=True)
        X, Y, Z = extract_transect(ima, coords, res)
    else:
        exit("invalid units")
    # fig.clf()

    # With the array, select values along a particular transect.
    # may require interpolation to a higher resolution grid.
    # return a 1d array of intensities and a 1d array of coords

    return transect


def hovmoller(transects):
    # make a hovmoller
    rgb = [t[2].reshape(1, *t[2].shape) for t in transects]
    imagearr = np.vstack(rgb)
    fig = plt.figure()
    fig.set_dpi(600)
    ax = fig.add_subplot(111)
    xmin = min(transects[0][1])
    xmax = max(transects[0][1])
    tmin = 0
    tmax = len(transects) * 15 / 60
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(tmin, tmax)
    # TODO: set interpolation for imshow
    ax.imshow(imagearr, extent=(xmin, xmax, tmin, tmax), origin='lower')
    ax.set_aspect('auto')
    ax.set_xlabel('Latitude, degrees')
    ax.set_ylabel('time, UTC')

    ax.set_title('Hovmoller')
    plt.savefig('hovmoller.png')

    return fig


def main(transect_coords=(0, 5, 10, 10), path='./images/', view='radagast'):
    if view == 'radagast':
        lower_left, upper_right = corner_coords()
    elif view == 'fennec':
        lower_left, upper_right = (-30.64, 2.60), (44.66, 46.31)
    else:
        exit("I don't know what this image is. Need some corner coords")
    image_coords = (lower_left, upper_right)

    files = file_list(path)
    args = (transect_coords, image_coords)
    kwargs = {'units': 'latlon', 'res': 500}
    transects = [transect(image, *args, **kwargs) for image in files]
    return transects


def map_test():
    # image edges for west africa seviri
    # (measured from the fennec browser hr panner)
    ll_lat = 2.60
    ll_lon = -30.64
    ur_lat = 46.31
    ur_lon = 44.66

    # ul_lat = 47.28
    # ul_lon = -54.00
    # lr_lat = 2.58
    # lr_lon = 27.08

    m = Basemap(projection='geos', llcrnrlon=ll_lon, llcrnrlat=ll_lat,
                urcrnrlon=ur_lon, urcrnrlat=ur_lat, resolution='h')
    image = Image.open('TEST.png')
    m.imshow(image)
    m.drawcountries()
    return m


def corner_coords(upper=(0, 30), right=(25, 0),
                  lower=(0, -5), left=(-25, 0)):
    """ RADAGAST
    More tricky to calculate as the region is bounded by points
    on the equator and the prime meridian, which are lat (-5,30),
    lon(-25,25). This method will work for any image for which we
    only have the coordinates of the edges rather than the corners.

    In order to compute the basemap projection that fits on top
    of the image we need to compute the corner coordinates. We do
    this by setting up a big basemap projection that is sure to
    contain all of the target area and then working out the
    corners within this.
    """
    m1 = Basemap(projection='geos', lon_0=0, resolution=None)
    # extract corner coords in lon / lat
    upper_right = m1(m1(*right)[0], m1(*upper)[1], inverse=True)
    lower_left = m1(m1(*left)[0], m1(*lower)[1], inverse=True)
    return lower_left, upper_right


def radagast_test():
    lower_left, upper_right = corner_coords()
    m = Basemap(projection='geos', lon_0=0, resolution='i',
                llcrnrlon=lower_left[0], llcrnrlat=lower_left[1],
                urcrnrlon=upper_right[0], urcrnrlat=upper_right[1])
    m.drawcoastlines()
    im = Image.open('TEST-radagast.jpg')
    m.imshow(im)
    plt.show()


# map radagast allows us to fit an image to a map projection, with
# which we can interact using either lat / lon or map projection
# coordinates.

# however, we still don't know how this translates into pixels in
# the image.

if __name__ == '__main__':
    F20_coords = [2.16, 13.5, -1.7, 19.6]

    parser = argparse.ArgumentParser(description="Produce a hovmoller from \
                                 a transect in some satellite images")
    parser.add_argument('-c', '--coordinates',
                        help='Coordinates x0 y0 x1 y1 of the \
                              transect start and end',
                        default=F20_coords, required=False, type=float, nargs=4)
    args = parser.parse_args()
    transect_coords = args.coordinates
    T = main(F20_coords, path='./images/')
    fig = hovmoller(T)

# bonus: intelligently interpolate the hovmoller such that
# coherent wave structures are picked out. requires selection of
# appropriate interpolation lines in image to create an
# interpolation field.
