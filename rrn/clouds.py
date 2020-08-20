from ctypes import c_int, c_double, CDLL  # visibility detector
import os.path

import numpy as np
from numpy.ctypeslib import ndpointer  # visibility detector


def load_config():
    """
    Carlo De Franchis's load_config function adapted to this wrapper.
    """

    working_directory_path = os.path.realpath(__file__)
    path = "/".join(working_directory_path.split("/")[:-1])

    paths = {
        "corrsift": os.path.join(path, "corrsift", "libcorrsift.so"),
        "angles": os.path.join(path, "visibility", "visibilitylib.so"),
    }

    for lib in paths.values():
        if not os.path.isfile(lib):
            raise FileNotFoundError(
                "Cannot find lib %s, please check this path" % lib)

    return paths


def grompone_timeseries(ims):
    """Apply Rafael Grompone von Gioi time series cloud (visibility) detector

    This detector takes as input a series of registered images.
    The images are single-channel.

    Args:
        ims (list of 2D numpy arrays): N input images of size X * Y
            (3D numpy array is also ok): idem.

    Returns:
        masks (list of 2D ndarrays, bool): list of masks for the input images
    """

    ims = np.array(ims)
    assert ims.ndim == 3

    N, Y, X = ims.shape  # X for width, Y for height, N for the number of images
    ims = np.ascontiguousarray(ims.astype(np.float64))
    mask = np.empty((N, Y, X), dtype=np.float64, order='C')

    here = os.path.dirname(os.path.abspath(__file__))
    lib_path = os.path.join(here, 'visibility', 'visibilitylib.so')
    lib = CDLL(lib_path)
    lib.visibility.argtypes = (c_int, c_int, c_int,
                               ndpointer(dtype=c_double, shape=(N, Y, X)),
                               ndpointer(dtype=c_double, shape=(N, Y, X)))
    lib.visibility(X, Y, N, ims, mask)

    return [mask[n,:,:] != 0 for n in range(mask.shape[0])]


def angles(im0, im1, lib_path=load_config()['angles']):
    """Wrapper for the C function that measure gradient angles differences

    Args:
        im0 (2D np array): image
        im1 (2D np array): image
        lib_path (str): path to the shared library angles.so file

    Return:
        2D np array measuring the difference of angles between the gradients of
                pixels in the two input images. Output is in [0, 1].
                Some not defined values are set to nan.
    """

    assert im0.ndim == 2 and im1.ndim == 2
    Y, X = im0.shape  # X is widht, Y height
    assert Y == im1.shape[0] and X == im1.shape[1]
    im0 = np.ascontiguousarray(im0.astype(np.float64))
    im1 = np.ascontiguousarray(im1.astype(np.float64))
    adiff = np.empty((Y, X), dtype=np.float64, order='C')
    lib = CDLL(lib_path)
    lib.angles.argtypes = (c_int, c_int,
                           ndpointer(dtype=c_double, shape=(Y, X)),
                           ndpointer(dtype=c_double, shape=(Y, X)),
                           ndpointer(dtype=c_double, shape=(Y, X)))
    lib.angles(X, Y, im0, im1, adiff)

    # NOTDEF values in adiff are set to -1024. Changing them to nan
    adiff[adiff == -1024] = 1  # np.nan

    # check the output; the cases below should never happen
    if np.any(adiff < 0):
        print('WARNING! Found negative values in the output of the angle '
                'function. They were clipped.')
        adiff = np.maximum(adiff, 0)
    if np.any(adiff > 1):
        print('WARNING! Found values superior to 1 in the output of the angle '
                'function. They were clipepd.')
        adiff = np.minimum(adiff, 1)

    return adiff


def cloud_masks(im_list, inputs):
    """
    Compute a ground visibility mask for each image in the provided list.
    The coverage value is returned too (in %).

    This uses R. Grompone von Gioi ground visibility detector

    Args:
        im_list: list of images (list of numpy arrays)
                give LUM images for the visibility-based detector
        inputs: list if inputs' names (used to avoid L8 images)

    Returns:
        cl_masks: (list of 2D boolean ndarrays) cloud masks
        cl_coverage: (list of floats) percentage of the image covered by clouds
    """

    cl_masks = []  # cloud masks (list of np arrays)
    cl_coverage = []  # cloud coverage (list of scalars)

    assert im_list[0].ndim == 2, 'This method needs luminance images'
    cl_masks = grompone_timeseries(im_list)

    cl_coverage = [np.sum(c) / np.size(c) * 100 for c in cl_masks]

    return cl_masks, cl_coverage
