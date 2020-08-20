"""
A python wrapper for T. Dagobert's SIFT correlation.
This wrapper is heavily inspired by one of Carlo de Franchis.
Test for example with:

import numpy as np
import rasterio

from corrsift_wrapper import corrsift

def test_corrsift(im0_path, im1_path):

    with rasterio.open(im0_path, 'r') as f0, \
         rasterio.open(im1_path, 'r') as f1:
        im0 = f0.read().squeeze()
        im1 = f1.read().squeeze()

        cmap = corrsift(im0.astype(np.float32), im1.astype(np.float32))

        iio.write('test.png', (cmap * 255).astype(np.uint8))

"""

import os
import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer


def load_config():
    """
    Carlo De Franchis's load_config function adapted to this wrapper.
    """

    working_directory_path = os.path.realpath(__file__)
    path = "/".join(working_directory_path.split("/")[:-1])

    paths = {
        "corrsift": os.path.join(path, "corrsift", "libcorrsift.so"),
    }

    for lib in paths.values():
        if not os.path.isfile(lib):
            raise FileNotFoundError(
                "Cannot find lib %s, please check this path" % lib)

    return paths


def corrsift(im0, im1, lib_path=load_config()["corrsift"]):
    """
    Python wrapper for Tristan Dagobert's C 'corrsift' function.

    The SIFT correlation gives the correlation coefficient for each pixel
    between the SIFT descriptors.

    Args:
        im0 (2D numpy array): image
        im1 (2D numpy array): image
        lib_path (str): path to the shared library corrsift.so file

    Return:
        2D numpy array representing the correlation map in [0, 1]
    """

    assert im0.ndim == 2 and im1.ndim == 2, 'input images must 2D'
    assert im1.shape == im0.shape, 'input images must have the same size'

    im0 = im0.astype(np.float32, order='C')
    im1 = im1.astype(np.float32, order='C')
    im_h, im_w = im0.shape

    # cmap: correlation map between im0 and im1
    cmap = np.zeros((im_h, im_w), dtype=np.float32)

    lib_corrsift = ctypes.CDLL(lib_path)
    lib_corrsift.corrsift.argtypes = (
            ndpointer(dtype=ctypes.c_float, shape=(im_h, im_w)),
            ndpointer(dtype=ctypes.c_float, shape=(im_h, im_w)),
            ndpointer(dtype=ctypes.c_float, shape=(im_h, im_w)),
            ctypes.c_int,
            ctypes.c_int)

    lib_corrsift.corrsift(
            im0,
            im1,
            cmap,
            im_h,
            im_w)

    # verifications; the cases below should never happen
    cmap_is_finite = np.isfinite(cmap)
    if not np.all(cmap_is_finite):
        print('WARNING! some values in cmap were not finite! '
              'They were replaced by zeros.')
        cmap[np.logical_not(cmap_is_finite)] = 0
    if np.any(cmap < -1):
        print('WARNING: values inferior to -1 corrsift! They were clipped.')
        cmap = np.maximum(cmap, -1)
    if np.any(cmap > 1):
        print('WARNING: values superior to 1 in corrsift! They were clipped.')
        cmap = np.minimum(cmap, 1)

    return cmap

