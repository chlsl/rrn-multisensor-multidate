#!/usr/bin/env python3
# vim: set fileencoding=utf-8

"""
Relative Radiometric Normalization (RRN)

Copyright (C) 2019-2020, Charles Hessel <charles.hessel@ens-paris-saclay.fr>



Official repository and future versions:

    https://github.com/chlsl/rrn-multisensor-multidate

If you use this code, please cite:

    @article{isprs-annals-V-2-2020-845-2020,
        AUTHOR = {Hessel, C. and Grompone von Gioi, R. and Morel, J. M. and Facciolo, G. and Arias, P. and de Franchis, C.},
        TITLE = {RELATIVE RADIOMETRIC NORMALIZATION USING SEVERAL AUTOMATICALLY CHOSEN REFERENCE IMAGES FOR MULTI-SENSOR, MULTI-TEMPORAL SERIES},
        JOURNAL = {ISPRS Annals of Photogrammetry, Remote Sensing and Spatial Information Sciences},
        VOLUME = {V-2-2020},
        YEAR = {2020},
        PAGES = {845--852},
        URL = {https://www.isprs-ann-photogramm-remote-sens-spatial-inf-sci.net/V-2-2020/845/2020/},
        DOI = {10.5194/isprs-annals-V-2-2020-845-2020}
    }
"""

import argparse
import os
import sys

import iio
import numpy as np

from rrn.clouds import cloud_masks
from rrn.contrast import stabilize_contrast, tonemap_sequence
from rrn.keys import find_key_images
from rrn.utils import (reparse, sp_switch, write_over_image, ndwi,
        rasterio_write, mkdir_p, load_images)


def main(inputs, outdir, bands='RGB',
         local_max=None, cloud_max_area=None, key_images=None,
         load_previous=False, no_overlay=False,
         acoefs_method=None, verbose=False, clouds=None, usecorrsift=False,
         tonemap=False, roughnorm=True):
    """
    Relative Radiometric Normalization (RRN)

    This function normalize the radiometry of image time series.
    A ground visibility detector is first applied. Then, cloudy images are
    removed. Reference images are found and they are used to stabilize the
    color and contrast of the remaining images. Tone-mapping is applied as a
    final and optional step.

    Args:
        inputs: list of str, input images paths
        outdir: str
        bands: list of str
        local_max
        cloud_max_area
        key_images
        load_previous
        no_overlay
        acoefs_method
        verbose
        clouds
        usecorrsift: bool, use sift correlation instead of angles for masks
        tonemap
        roughnorm

    Returns:
        None. The images are saved to the disk in 'outdir'
    """

    # re-parse "band" to avoid duplicates or incompatible combinations
    bands = reparse(bands)  # band is a list

    # input list (inpl): separate and sort the input files (paths) given
    inpl = {b: sorted([i for i in inputs if '_' + b in i]) for b in bands}
    assert all(len(inpl[b]) == len(inpl[bands[0]]) for b in bands), \
        'expecting equal numbers of images in all input bands.'

    # sort images by date, whatever the directory name there're in
    for b in inpl:
        inputs_basename = [os.path.basename(i) for i in inpl[b]]
        inpl.update({b: [p for i, p in sorted(zip(inputs_basename, inpl[b]))]})
        # remove duplicates
        inputs_basename = [os.path.basename(i) for i in inpl[b]]
        dumpdup = []
        for n, (p, fullp) in enumerate(zip(inputs_basename, inpl[b])):
            if n == 0:
                dumpdup.append(fullp)
                continue
            else:
                curr_path = p
                print(p)
                prev_path = inputs_basename[n-1]
                if prev_path[:20] in curr_path:
                    if '_L1C_' in curr_path and '_L2A_' in prev_path:
                        continue
                        # inpl[b].pop(n)
                        print(f'{curr_path} removed, cause duplicated image')
                    elif '_L2A_' in curr_path and '_L1C_' in prev_path:
                        dumpdup.pop(-1)
                        dumpdup.append(fullp)
                        # inpl[b].pop(n+1)
                        print(f'{prev_path} removed, cause duplicated image')
                else:
                    dumpdup.append(fullp)
        inpl.update({b: dumpdup.copy()})

    print('List of sorted input images (first band only):')
    for p in inpl[bands[0]]: print(p)

    # read input images
    imgl = {}
    for b in bands:
        img_list, msk_list = load_images(inpl[b])
        imgl.update({b: img_list, b+'msk':msk_list})

    # stop when no image was loaded
    if len(imgl[bands[0]]) == 0:
        raise ValueError(
                'Could not read the input images.\n'
                'The image names are supposed to contain "_XXX" where XXX '
                'is RGB or B02, B03 etc. Use the argument "--band" or "-b" to '
                'specify which bands to load.')

    # create output directory and filenames
    print('creating directory %s' % outdir)
    mkdir_p(outdir)

    # create RGB images if needed
    if 'RGB' not in imgl:
        imgl.update({'RGB': [np.array([r, g, b]) for b, g, r
                             in zip(imgl['B02'], imgl['B03'], imgl['B04'])]})
        # files in inpl['RGB'] may very well not exist but this is useful anyway
        inpl.update({'RGB': [i.replace('_B02', '_RGB') for i in inpl['B02']]})

    # create dictionnary containing output images' names. Its name is outl
    outl = {'RGB': [outdir + '/' + os.path.basename(i) for i in inpl['RGB']]}

    # create luminance images
    assert 'RGB' in imgl, 'at this point the RGB images should exist!'
    imgl.update({'LUM': [np.sum(i,axis=0)/3 for i in imgl['RGB']]})

    # normalize the input images
    # This is useful for visualization purposes, but not required in our method.
    # after this, the RGB images are roughly in the range [0, 1]. Beware: there
    # is no clipping, so many values are still out of this range!
    if roughnorm:
        rough_norm_params_S2 = [64, 1448-64]
        rough_norm_params_L8 = [4776, 14671-4776]
        rough_norm_params_L4 = [0, 255]
        possible_bands = {'B02', 'B03', 'B04', 'B08', 'RGB', 'LUM'}
        for b in possible_bands:
            if not imgl.get(b): continue
            print(f'Rough normalization (using constants) for band {b}')
            roughnormlist = []
            for i, p in zip(imgl[b], outl['RGB']):
                if '_L4_' in p or '_L5_' in p or '_L7_' in p:
                    roughnormlist.append((i - rough_norm_params_L4[0]) / rough_norm_params_L4[1])
                elif '_L8_' in p:
                    roughnormlist.append((i - rough_norm_params_L8[0]) / rough_norm_params_L8[1])
                elif '_S2A_' in p or '_S2B_' in p:
                    roughnormlist.append((i - rough_norm_params_S2[0]) / rough_norm_params_S2[1])
                else:
                    raise ValueError('Unable to perform the rough '
                            'normalization. This is probably because the '
                            'filenames do not contain the satellite name, or '
                            'because this satellite is not handled yet.\n'
                            'Try again with option --no-roughnorm.')
            imgl.update({b: roughnormlist})


    # compute water coefficient if possible
    if 'B08' in imgl:
        print('computing the normalized difference water index...')
        if 'B02' in imgl:
            imgl.update({'WTR': [ndwi(g, ir) > 0
                                 for g, ir in zip(imgl['B03'], imgl['B08'])]})
        elif 'RGB' in imgl:
            imgl.update({'WTR': [ndwi(rgb[1], ir) > 0
                                 for rgb, ir in zip(imgl['RGB'], imgl['B08'])]})
        # also add the file names into the inputs_list (inpl) dictionary
        inpl.update({'WTR': [i.replace('_B02', '_WTR')
                             for i in inpl['B08']]})
        outl.update({'WTR': [outdir + '/' + os.path.splitext(
                             os.path.basename(i))[0].replace('_B08',
                             '_WTR') + '.png' for i in inpl['B08']]})
        for path, imag in zip(outl['WTR'], imgl['WTR']):
            iio.write(path, 255 * imag.astype(np.uint8))

    # find corresponding cloud masks if --load-previous flag is enabled
    if load_previous:
        clouds = [outdir + '/'
                + os.path.splitext(os.path.basename(path))[0].replace(
                '_RGB', '_CLD') + '.png'
                for path in inpl['RGB']]
        if not all(os.path.isfile(f) for f in clouds):
            print('NOTE: the option --load-previous was used but some visibility'
                  ' masks could not be found (could be one or all). '
                  'Recomputing them all.')
            clouds = None

    # load cloud masks
    if clouds:  # clouds (path to cloud masks) can also be passed in option
        inpl.update({'CLD': clouds})

        # load cloud masks
        print('loading already computed visibility masks...')
        imgl.update({'CLD': load_images(inpl['CLD'])[0]}, option=True)

        # ensure that they are boolean
        imgl['CLD'] = [(cld != 0) if (cld.ndim == 2)  # for black/white masks
                       else (cld[0,:,:] > cld[1,:,:])  # for green/red masks
                       for cld in imgl['CLD']]
        cl_coverage = [np.sum(m) / np.size(m) * 100 for m in imgl['CLD']]

    # compute the cloud masks, unless they already exist
    if not imgl.get('CLD'):  # if they are not already there (see above)
        print('computing visibility masks...')
        if not imgl.get('LUM'):
            imgl.update({'LUM': [np.sum(i,axis=0)/3 for i in imgl['RGB']]})
        cl_masks, cl_coverage = cloud_masks(imgl['LUM'], inpl['RGB'])
        imgl.update({'CLD': cl_masks})

    assert isinstance(imgl['CLD'], list)
    assert isinstance(cl_coverage, list)
    assert all(isinstance(i, (np.ndarray, np.generic)) for i in imgl['CLD'])

    # The lists are then split between sunny and cloudy images
    im_list_sunny = []
    im_list_cloudy = []
    cl_masks_sunny = []
    cl_masks_cloudy = []
    inputs_sunny = []
    inputs_cloudy = []
    cl_coverage_sunny = []
    cl_coverage_cloudy = []

    for im, cl, path, cover in zip(imgl['RGB'], imgl['CLD'], inpl['RGB'], cl_coverage):
        if cover < cloud_max_area:
            im_list_sunny.append(im)
            cl_masks_sunny.append(cl)
            inputs_sunny.append(path)
            cl_coverage_sunny.append(cover)
        else:
            im_list_cloudy.append(im)
            cl_masks_cloudy.append(cl)
            inputs_cloudy.append(path)
            cl_coverage_cloudy.append(cover)

    imgl.update({'CLD-SUN':cl_masks_sunny})  # cloud masks in dict

    if 'WTR' in imgl:  # water masks in two separated dict, sunny and cloudy
        imgl.update({'WTR-SUN':[w for w, c in zip(imgl['WTR'], cl_coverage)
                                if c < cloud_max_area]})
        imgl.update({'WTR-CLD':[w for w, c in zip(imgl['WTR'], cl_coverage)
                                if c >= cloud_max_area]})

    # input file names
    inpl.update({'RGB-SUN': inputs_sunny})

    # output file names
    outl.update({'RGB-SUN': [outdir + '/' +
        os.path.splitext(os.path.basename(i))[0] + '_sunny' +
        os.path.splitext(os.path.basename(i))[1] for i in inputs_sunny]})
    outl.update({'RGB-CLD': [outdir + '/' +
        os.path.splitext(os.path.basename(i))[0] + '_cloudy' +
        os.path.splitext(os.path.basename(i))[1] for i in inputs_cloudy]})

    # writing sunny and cloudy images in output directory for comparisons
    # also write the cloud masks
    for t, img_t, msk_t, inp_t in zip(['sunny', 'cloudy'],
                                      [im_list_sunny, im_list_cloudy],
                                      [cl_masks_sunny, cl_masks_cloudy],
                                      [inputs_sunny, inputs_cloudy]):
        if not img_t:  # there is no sunny/cloudy images
            continue
        print('writing %s images (and their visibility masks) in %s' % (t, outdir))
        for im, ma, im_path in zip(img_t, msk_t, inp_t):
            im_name, im_ext = os.path.splitext(os.path.basename(im_path))
            im_path_out = os.path.join(outdir, im_name + '_' + t + im_ext)
            rasterio_write(im_path_out, im.transpose((1, 2, 0)))
            ma_path_out = os.path.join(outdir, im_name.replace('_RGB', '_CLD')
                    + '_' + t + '.png')
            iio.write(ma_path_out, ma.astype(np.uint8) * 255)  # masks are bool
        if tonemap:
            img_t = tonemap_sequence(img_t)
            for im, im_path in zip(img_t, inp_t):
                im_name = os.path.splitext(os.path.basename(im_path))[0]
                im_path_out = os.path.join(outdir, im_name+'_'+t+'_tonemap.png')
                iio.write(im_path_out, im.transpose((1, 2, 0)))

    # print infos about clouds coverage and removed images
    if len(inputs_cloudy) > 0:
        print('%d image%s removed (over %d):'
              % (len(inputs_cloudy), 's' if len(inputs_cloudy) > 1 else '',
                 len(inpl['RGB'])))
        for n, impath in enumerate(inputs_cloudy):
            print('%s -- invisible ground %5.2f%%' % (impath,
                                                    cl_coverage_cloudy[n]))
        if len(im_list_sunny) > 1:
            print('continuing with %d images:' % len(im_list_sunny))
            for n, impath in enumerate(inputs_sunny):
                print('%s -- invisible ground %5.2f%%' % (impath,
                                                        cl_coverage_sunny[n]))

    # exit if all images are cloudy
    if not imgl['CLD-SUN']:
        print('It seems that all images were detected with too little visibility. They were '
              'removed from the sequence, which is now empty and thus cannot '
              'be stabilized.\nYou may want to try to set the '
              '"--no_visible_ground_area" option to a higher value. For now it is set '
              f'to {cloud_max_area}%.')
        sys.exit(0)

    # exit if number of remaining images is 0 or 1
    if len(im_list_sunny) == 0:
        sys.exit('All images have a too small visible area. Stopping now. \n'
                 'You might want to try a higher no-visible-ground threshold.')
    if len(im_list_sunny) == 1:
        sys.exit('Only one image left. '
                 'This is not enough for contrast stabilization. '
                 'Stopping now. \n'
                 'You might want to try a higher no-visible-ground threshold.')

    if key_images is None:  # then choose the key images
        print('Finding key images:')
        keys = find_key_images(im_list_sunny, cl_coverage_sunny,
                               local_max, inputs_sunny, imgl['CLD-SUN'])
    else:
        keys = key_images

    assert len(keys) > 0, 'the stabilization requires at least one key image'

    # print key images names
    print('found %d key image%s:' % (len(keys), sp_switch(keys)))
    for n, k in enumerate(keys):
        print('key image %3d, index %3d, file %s' % (n, k, inputs_sunny[k]))

    # create output images' names
    output_list = [outdir + '/'
                   + os.path.splitext(
                       os.path.basename(p))[0].replace('_RGB', '_MSK')
                   + '.png' for p in inputs_sunny]

    # get masks, to be handed to the stabilization function
    if 'WTR-SUN' in imgl:
        masks = [np.logical_not(np.logical_or(c, w)) for c in imgl['CLD-SUN']
                                                     for w in imgl['WTR-SUN']]
    else:
        masks = [np.logical_not(c) for c in imgl['CLD-SUN']]

    # stabilize image list using keys
    im_list_stab = stabilize_contrast(im_list_sunny, keys,
                                      output_list, masks,
                                      acoefs_method=acoefs_method,
                                      usecorrsift=usecorrsift)

    # write stabilized images in output directory
    assert len(inputs_sunny) == len(im_list_stab)
    for n, (im_path, im_stab) in enumerate(zip(inputs_sunny, im_list_stab)):
        imname, imext = os.path.splitext(os.path.basename(im_path))
        im_path_out = outdir + '/' + imname + '_stab' + imext
        print('writing %s' % im_path_out)
        # overlay with date
        if not no_overlay:
            im_stab = write_over_image(im_stab,
                                       imname[:10]
                                       + (' *' if n in keys else '')
                                       + (' L8' if '_L8_' in imname else '')
                                       + (' S2A' if '_S2A_' in imname else '')
                                       + (' S2B' if '_S2B_' in imname else ''),
                                       copy=True)
        rasterio_write(im_path_out, im_stab.transpose((1, 2, 0)))

    # im_list_stab: list of stabilized sunny images
    # inputs_sunny: list of pathnames, restricted to the sunny images
    if tonemap:
        # tonemap sequence: from whatever range to [0, 255] (clipped)
        im_list_tonemap = tonemap_sequence(im_list_stab)
        for n, (im, path) in enumerate(zip(im_list_tonemap, inputs_sunny)):
            imname = os.path.splitext(os.path.basename(path))[0]
            path = os.path.join(outdir, imname+'_stab_tonemap.png')
            if not no_overlay:
                im = write_over_image(im, imname[:10]
                        + (' L8 (pansharpened and upsampled)'
                           if '_L8_' in imname else '')
                        + (' S2A' if '_S2A_' in imname else '')
                        + (' S2B' if '_S2B_' in imname else '')
                        + (' L1C' if '_L1C_' in imname else '')
                        + (' L2A' if '_L2A_' in imname else '')
                        + (' -- KEY' if n in keys else ''), copy=True)
            iio.write(path, im.transpose((1, 2, 0)))

    # move non-essential saved files in a subdirectory
    mkdir_p(os.path.join(outdir, 'supplement'))
    files_to_move = [f for f in os.listdir(outdir) if
                     any(pattern in f for pattern in ('_MSK0_final.png',
                                                      '_MSK1_final.png',
                                                      '_CLD_sunny.png',
                                                      '_RGB_sunny.tif',
                                                      '_RGB_sunny_tonemap.tif',
                                                      '_RGB_sunny_stab.png',
                                                      '_PLT0c0.png',
                                                      '_PLT0c1.png',
                                                      '_PLT0c2.png',
                                                      '_PLT1c0.png',
                                                      '_PLT1c1.png',
                                                      '_PLT1c2.png',
                                                      '_RANSAC0c0.png',
                                                      '_RANSAC0c1.png',
                                                      '_RANSAC0c2.png',
                                                      '_RANSAC1c0.png',
                                                      '_RANSAC1c1.png',
                                                      '_RANSAC1c2.png'))]
    for f in files_to_move:
        os.rename(os.path.join(outdir,f), os.path.join(outdir, 'supplement', f))


def cli():
    """
    Command line interface
    """

    # Program description
    parser = argparse.ArgumentParser(description='Relative Radiometric '
                                     'Normalization')

    # Mandatory, positional
    pos = parser.add_argument_group('positional arguments (mandatory)')
    pos.add_argument('imlist', nargs='+',
                     help='list of images')

    # Mandatory, not positional
    man = parser.add_argument_group('mandatory arguments')
    man.add_argument('-b', '--band', nargs='+', type=str, default='RGB',
                     choices={'B02', 'B03', 'B04', 'B08', 'RGB'},
                     help='space-separated list of bands to load (in the '
                     'provided files, see "imlist"). Specifying RGB implies '
                     'that the B02, B03 and B04 bands were merged '
                     'beforehand. Default is RGB.')
    man.add_argument('-o', '--outdir', required=True, type=str,
                     help='output directory (created if nonexistent)')

    # Method's parameters
    prm = parser.add_argument_group('parameters of the method')
    prm.add_argument('-m', '--local-max', type=int, default=9,
                     help='a notion of how "local" the local maxima are. '
                          'Default is 9.')
    prm.add_argument('--no-visible-ground-area', type=int, default=25,
                     help='threshold on the proportion (percentage) of not '
                     'visible ground pixels. Images whose proportion of not '
                     'visible ground pixels is above this threshold are '
                     'rejected.  Default is 25.')

    # By-pass some of the method's steps
    byp = parser.add_argument_group("by-pass some of the method's steps")
    byp.add_argument('--key-images', nargs='+', type=int, default=None,
                     help='list of index for key images.')
    byp.add_argument('-c', '--clouds', nargs='+', type=str, default=None,
                     help='input cloud masks. If not provided, they are '
                             'computed using R. Grompone von Gioi '
                             'ground visibility detector.')
    byp.add_argument('--load-previous', action='store_true', default=False,
                     help='use this flag to load previously computed '
                          'intermediary results rather than computing them. '
                          'Useful for debugging mainly.')

    # Alternatives for some of the method's steps
    alt = parser.add_argument_group('alternatives for some steps')
    alt.add_argument('--acoefs-method', type=str, default='ransac',
                     choices={'L1', 'L1constrained', 'median', 'ransac'},
                     help='computation method for the affine coefficients. '
                          'Default is ransac.')
    alt.add_argument('--usecorrsift',
                     dest='usecorrsift', action='store_true',
                     help='flag to replace corrsift (slow) by the angles '
                          '(faster). The angle error is less accurate. '
                          'Default is True.')
    alt.add_argument('--dont-usecorrsift',
                     dest='usecorrsift', action='store_false',
                     help='set to False the usecorrsift flag.')
    alt.set_defaults(usecorrsift=False)

    # Options
    opt = parser.add_argument_group('more options!')
    opt.add_argument('--verbose', action='store_true', default=False,
                     help='print non-critical information')
    opt.add_argument('--no-overlay', action='store_true', default=False,
                     help='use this flag to remove the bottom overlay '
                          'with the date added on the stabilized images.')
    opt.add_argument('--roughnorm', dest='roughnorm', action='store_true',
                     help='flag to apply a rough normalization to the input '
                     'images. Default is True.')
    opt.add_argument('--no-roughnorm', dest='roughnorm', action='store_false',
                     help='Set to False the roughnorm flag.')
    opt.set_defaults(roughnorm=True)
    opt.add_argument('--tonemap', dest='tonemap', action='store_true',
                     help='flag to apply a tonemapping step on the output '
                     'images. Default is True. (use --no-tonemap to deactivate)')
    opt.add_argument('--no-tonemap', dest='tonemap', action='store_false',
                     help='set to False the tonemap flag.')
    opt.set_defaults(tonemap=True)

    # Note: argparse infers var name from arg name and replaces any '-' by '_'
    args = parser.parse_args()

    # print('RECAP OF OPTIONS:')
    # print(args)

    main(args.imlist, args.outdir, args.band,
         args.local_max, args.no_visible_ground_area, args.key_images,
         args.load_previous, args.no_overlay, args.acoefs_method,
         args.verbose, args.clouds, args.usecorrsift,
         args.tonemap, args.roughnorm)


if __name__ == '__main__':
    cli()
