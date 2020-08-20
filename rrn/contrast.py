import iio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import ponomarenko
from scipy.ndimage import uniform_filter

from rrn.clouds import angles
from rrn.corrsift_wrapper import corrsift
from rrn.model_parameters_estimation import affine_coefficients
from rrn.utils import write_over_image


def tonemap_sequence(images, clip=True):
    """Tone-mapping of a sequence, from an arbitrary range to [0, 255].

    Args:
        images: list of 3D ndarrays, sequence of images

    Return:
        tonemapped sequence
    """

    assert images, 'Tone-mapping requires a non-empty list'

    lu_ims = [np.mean(i, axis=0) for i in images]
    percentiles = [(np.percentile(i[2:-2, 2:-2], .1),
                    np.percentile(i[2:-2, 2:-2], 99)) for i in lu_ims]
    global_minper, global_maxper = np.median(np.array(percentiles), axis=0)
    global_dyn = global_maxper - global_minper;

    print(f'Tone-mapping a series of {len(images)} images.\n'
          f'Their median 1 and 99 percentiles are {global_minper:.4f}, '
          f'{global_maxper:.4f} (i.e. dynamic is roughly {global_dyn:.4f}).')

    if clip:
        tm_ims = [(np.clip((i - global_minper) / global_dyn, 0, 1)**(3/4)
                   * 255).astype(np.uint8) for i in images]
    else:
        raise NotImplementedError

    return tm_ims


def stabilize_contrast(im_list, keys, output_list, masks_external=None,
                       acoefs_method=None, usecorrsift=False):
    """Stabilize contrast using key images in the list

    Args:
        im_list: (list of 3D np arrays) color input images
        keys: (list of scalars) index for the key images (used as target for
                contrast corrections)
        output_list: (list of str) names for the output images (one per image)
        masks_external: (list of 2D np arrays) masks corresponding to the input
            image list (same size, same order). Warning: use True for pixels
            that can be used and False for the others.

    Returns:
        (list of 3D np arrays) images with stabilized contrast
    """

    assert len(im_list) == len(output_list)

    if im_list[0].ndim == 2:
        im_size = im_list[0].shape[0:2]
    elif im_list[0].ndim == 3:
        im_size = im_list[0].shape[1:3]
    else:
        raise ValueError('Incorrect number of dimensions for images in list')

    for k in keys:
        iio.write(output_list[k].replace('MSK', 'MSK0_final'),
                write_over_image(np.full(im_size, 255, dtype=np.uint8),
                    'key'))
        iio.write(output_list[k].replace('MSK', 'MSK1_final'),
                write_over_image(np.full(im_size, 255, dtype=np.uint8),
                    'key'))

    for n in range(len(im_list)):
        if n < keys[0]:
            iio.write(output_list[n].replace('MSK', 'MSK1_final'),
                    write_over_image(np.full(im_size, 255, dtype=np.uint8),
                        'placeholder'))
        if n > keys[-1]:
            iio.write(output_list[n].replace('MSK', 'MSK0_final'),
                    write_over_image(np.full(im_size, 255, dtype=np.uint8),
                        'placeholder'))

    # estimate noise using ponomarenko, for ransac
    print('Estimating noise in the sequence using ponomarenko...')
    stdnoise_all = [ponomarenko.estimate_noise(np.sum(i, axis=0) / i.shape[0],
                    num_bins=3)[1][1] for i in im_list]
    stdnoise = np.median(stdnoise_all)
    print(f'The median value of the noise (std) over the {len(im_list)} images '
          f'of the sequence is {stdnoise}.')

    # computing luminance images
    lum = np.empty((len(im_list), im_list[0].shape[1], im_list[0].shape[2]))
    for n in range(lum.shape[0]):
        lum[n, :, :] = np.sum(im_list[n], axis=0) / 3

    # initialize coefficients with nans
    assert im_list[0].ndim == 3 and im_list[0].shape[0] == 3, 'expecting color'
    coefs = np.full((len(im_list), 3, 2), np.nan)

    print('Normalizing the radiometry using the keys images...')

    # 1) compute masks
    # 2) compute affine transformations, with interpolation when needed
    # 3) Apply them

    for m, k in enumerate(keys):
        # prev_key and next_key are the key images except at the beginning and
        # end of the sequence where their value is set so that the sift
        # correlation is correctly computed.
        print('proceeding with key image at index %d' % k)
        prev_key = keys[m - 1] if m > 0 else - 1
        next_key = keys[m + 1] if m < len(keys) - 1 else len(im_list)

        for n in range(prev_key + 1, next_key):

            assert im_list[n].ndim == 3, 'expecting rgb image'
            print('img %3d, pair (%3d, %3d)' % (n, k, n), end=' -- ')

            if n == k:
                print('skipping, do not need stabilization')
                coefs[n, :, :] = [1., 0.]
                continue

            # Step 1/3: Compute mask
            # ----------------------

            # set mask to - True for values that can be used in the regression
            #             - False for non-usable pixels (clouds, water, etc.)

            if usecorrsift:  # default is False
                print('corrsift...', end='', flush=True)
                mk_sift = corrsift(lum[k], lum[n])
                mk_sift = np.maximum(mk_sift, 0)  # correlation in [-1, 1]
                mk_sift = mk_sift > np.percentile(mk_sift, 75)
            else:
                print('angles...', end='', flush=True)
                # angles return values in [0, 1] and some nans.
                # 0 means same angle, so it's the best result one can expect
                # (this is reversed with reference to the sift correlation)
                ang = angles(lum[k].astype(np.float64),
                                 lum[n].astype(np.float64))
                ang_avg = uniform_filter(ang, size=3, mode='nearest')
                threshold = np.percentile(ang_avg, 10)
                mk_sift = ang_avg < threshold  # 0.1

            mask = mk_sift

            # write masks
            # the output_list list contains names like <date>_(...)_L1C_MSK.png
            maskname = output_list[n]
            keyorder = 0 if n < k else 1
            ### if masks_external is not None:
            ###     iio.write(maskname.replace('MSK', f'MSK{keyorder}_ext'),
            ###               (255 * mk_ext).astype(np.uint8))
            ### iio.write(maskname.replace('MSK', f'MSK{keyorder}_sift'),
            ###           (255 * mk_sift).astype(np.uint8))
            iio.write(maskname.replace('MSK', f'MSK{keyorder}_final'),
                      (255 * mask).astype(np.uint8))

            if not mask.any():  # if the mask is empty
                print('WARNING: mask is empty. Impossible to compute coefs... '
                      'Skipping this image. It will remain unchanged.')
                coefs[n, :, :] = [1., 0.]
                continue

            print('done', end=' -- ')

            # Step 2/3: compute affine coefficients
            # -------------------------------------

            print(f'coefs (with {acoefs_method})...', end='', flush=True)
            a = np.ones(3)
            b = np.zeros(3)
            for c in range(3):  # three color channels
                a[c], b[c], inliers = affine_coefficients(target=im_list[k][c],
                                                          image=im_list[n][c],
                                                          mask=mask,
                                                          method=acoefs_method,
                                                          stdnoise=stdnoise)

                # plot coefficients and save masks
                target_init = im_list[k][c]
                source_init = im_list[n][c]
                target_out = target_init[np.logical_not(mask)]
                source_out = source_init[np.logical_not(mask)]
                target_nor = target_init[mask][np.logical_not(inliers)]
                source_nor = source_init[mask][np.logical_not(inliers)]
                target_inl = target_init[mask][np.array(inliers)]
                source_inl = source_init[mask][np.array(inliers)]
                if c == 0: color_nor, color_inl = 'tab:red', 'tab:orange'
                elif c == 1: color_nor, color_inl = 'tab:green', 'tab:olive'
                elif c == 2: color_nor, color_inl = 'tab:blue', 'tab:cyan'
                color_out = 'tab:gray'
                fig, axs = plt.subplots(1, 1)
                axs.scatter(source_out, target_out, marker='.', color=color_out, s=1, alpha=.2)
                axs.scatter(source_nor, target_nor, marker='.', color=color_nor, s=1, alpha=.2)
                axs.scatter(source_inl, target_inl, marker='.', color=color_inl, s=1, alpha=.2)
                axs.plot([0, 2.5], [b[c], a[c]*2.5 + b[c]], color=color_nor, alpha=1)
                axs.axis([0, 2.5, 0, 2.5])
                axs.set_aspect('equal', 'box')
                fig.savefig(maskname.replace('MSK', f'PLT{keyorder}c{c}'),
                            bbox_inches='tight', transparent=True)
                plt.close(fig)
                # same inliers as mask too
                # mask as same size as image
                ransac_pos_inliers = np.nonzero(mask.flatten())[0]
                ransac_pos_inliers = ransac_pos_inliers[np.array(inliers)]
                ransac_mask_inliers = np.zeros((target_init.size,), dtype=np.uint8)
                ransac_mask_inliers[ransac_pos_inliers] = 255
                ransac_mask_inliers = ransac_mask_inliers.reshape(target_init.shape)
                iio.write(maskname.replace('MSK', f'RANSAC{keyorder}c{c}'),
                          ransac_mask_inliers)

            print('done: '
                  f'R=({a[0]:4.2f}, {b[0]:4.2f}); '
                  f'G=({a[1]:4.2f}, {b[1]:4.2f}); '
                  f'B=({a[2]:4.2f}, {b[2]:4.2f}) -- ', end='')

            # interpolate coefficients
            if m == 0 and n < k:
                coefs[n, :, :] = np.array([a, b]).T
                print('no interp for first imgs')
            elif m == len(keys) - 1 and n > k:
                coefs[n, :, :] = np.array([a, b]).T
                print('no interp for last imgs')
            elif n > k:
                d1 = next_key - n
                d2 = next_key - k
                coefs[n, :, :] = d1 / d2 * np.array([a, b]).T
                print('lin interp; waiting for 2nd img')
            elif n < k:
                d1 = n - prev_key
                d2 = k - prev_key
                coefs[n, :, :] += d1 / d2 * np.array([a, b]).T
                print('lin interp; '
                      f'red  =({coefs[n, 0, 0]:6.4f}, {coefs[n, 0, 1]:6.4f}); '
                      f'green=({coefs[n, 1, 0]:6.4f}, {coefs[n, 1, 1]:6.4f}); '
                      f'blue =({coefs[n, 2, 0]:6.4f}, {coefs[n, 2, 1]:6.4f}).')
            else:
                print('[ERROR] In coefs interp: this case should never happen.')

    # Step 3/3: apply computed affine coefficients
    # --------------------------------------------

    im_list_stab = []
    im_stab = np.zeros(im_list[0].shape)
    for n, im in enumerate(im_list):
        for c in range(3):  # three color channels
            im_stab[c, :, :] = coefs[n, c, 0] * im[c, :, :] + coefs[n, c, 1]
            if np.isnan(im_stab[c, :, :]).any():
                print('WARNING! nans in stabilized image... '
                      'idx %3d, coefs a=%4.1f, b=%8.1f' % (n, coefs[n, c, 0],
                                                           coefs[n, c, 1]))
        im_list_stab.append(im_stab.copy())

    return im_list_stab

