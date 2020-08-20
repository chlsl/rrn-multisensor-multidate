import os
from scipy.ndimage import uniform_filter  # local std to find key images
import numpy as np


def find_key_images(im_list_sunny, cl_coverage_sunny,
                    local_max, inputs_sunny,
                    cloud_masks_sunny):
    """Find highest quality images in the time series

    Args:
        im_list_sunny: list of ndarrays, sunny input images
        cl_coverage_sunny: list of ndarrays, sunny input images' cloud masks
        local_max: int, temporal neighborhood size
        inputs_sunny: list of str, sunny inputs' paths
        cloud_masks_sunny: list of str, sunny inputs' cloud masks

    Returns:
        list of indices of the key images
    """

    # "intrinsic" quality (resolution and correction level <==> favor S2-L2A)
    assert len(inputs_sunny) == len(im_list_sunny)
    ms_resolution = [0.1 if '_L8_' in i or '_L1C_' in i else 1
                     for i in inputs_sunny]

    # proportion of visible ground
    ms_sunnyness = list(1 - np.array(cl_coverage_sunny) / 100)

    # local contrast
    ms_contrast = []
    assert len(cloud_masks_sunny) == len(im_list_sunny)
    for n, im in enumerate(im_list_sunny):
        assert im.shape[0] == 3, 'image should have three channels'
        lum = np.sum(im, axis=0) / 3
        # print('lum.shape = {}'.format(lum.shape))
        # print('lum[~cl_masks_sunny[n]].shape = {}'.format(lum[~cl_masks_sunny[n]].shape))
        if np.isnan(lum).any():
            print('found nan values in the image.')

        # replace nan and inf by image's mean value (for local_std)
        lum[~np.isfinite(lum)] = np.mean(lum[np.isfinite(lum)])
        local_std = np.sqrt(np.maximum(uniform_filter(lum ** 2, size=15)
                                       - uniform_filter(lum, size=15) ** 2, 0))
        std = np.mean(local_std[cloud_masks_sunny[n]])
        std = std / np.std(lum[cloud_masks_sunny[n]])
        ms_contrast.append(std)

    assert len(cl_coverage_sunny) == len(ms_contrast), \
        'len(cl_coverage_sunny) = {}; len(ms_contrast) = {}'.format(
            len(cl_coverage_sunny), len(ms_contrast))

    # compute the quality score of each image
    scores = np.array(ms_sunnyness) \
             * np.array(ms_contrast) \
             * np.array(ms_resolution)

    # find key images (higher scores in a temporal neighborhood)
    # the list "keys" is the index in the sunny sequences
    keys = [n for n, i in enumerate(scores)
            if (i >= scores[max(0, n-local_max):min(len(scores), n+local_max+1)]
                ).all()]

    # print scores and announce key images
    for n, i in enumerate(inputs_sunny):
        print('image %3d: visibility %.4f; local contrast %.4f; '
              'accuracy %.1f  --  score %.4f  --  key?  %s  (%s)'
              % (n, ms_sunnyness[n], ms_contrast[n],
                 ms_resolution[n], scores[n], 'YES' if n in keys else 'NO ',
                 os.path.basename(i)))

    return keys
