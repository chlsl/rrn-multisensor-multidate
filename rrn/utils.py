import errno  # for mkdir_p
import warnings
import os

import numpy as np
from PIL import Image, ImageDraw, ImageFont
import rasterio


warnings.filterwarnings("ignore",
                        category=rasterio.errors.NotGeoreferencedWarning)


def reparse(band):
    band = [band] if type(band) is not list else band  # ensure list
    band = set(band)
    if 'RGB' in band:
        band -= {'B02', 'B03', 'B04'}
    elif 'B02' in band:
        band |= {'B02', 'B03', 'B04'}
    return sorted(list(band))


def sp_switch(l, s='', p='s'):
    """Switch between singular and plural in function of l's length

    Args:
    l: list
    s: singular string (default value '')
    p: plural string (default value 's')

    Return:
    correct string
    """
    return p if len(l) > 1 else s


def write_over_image(image, text, copy=False):
    """
    Write text on the bottom left corner of an image (numpy array)
    :return:

    Thanks Giles Booth for the font ChiKareGo2.
    See his blog post: http://www.suppertime.co.uk/blogmywiki/2017/04/chicago/
    """
    # we cannot directly use the image, because PIL only handle uint8 images
    # instead we write the text in a "strip" and deduce a mask from it

    assert image.ndim in (2, 3), 'Can write text only on gray or color images'

    strip_h = 17
    strip_w = image.shape[2] if image.ndim == 3 else image.shape[1]
    strip = Image.fromarray(np.zeros((strip_h, strip_w), dtype=np.uint8))

    strip_labeled = ImageDraw.Draw(strip)
    font = ImageFont.truetype(os.path.dirname(os.path.abspath(__file__))
                              + '/fonts/ChiKareGo2.ttf', size=16)
    strip_labeled.text((4, 1), text, fill=255, font=font)

    if copy:
        imout = image.copy()

    if image.ndim == 3:
        text_mask = np.broadcast_to(np.array(strip) == 255,
                                    (image.shape[0], strip_h, strip_w))
        image_crop = imout[:, -strip_h:, :] if copy else image[:, -strip_h:, :]
    else:
        text_mask = np.array(strip) == 255
        image_crop = imout[-strip_h:, :] if copy else image[-strip_h:, :]

    bkgd_mask = np.logical_not(text_mask)
    image_crop[text_mask] = 2**16 - 1
    image_crop[bkgd_mask] = (.6 * image_crop[bkgd_mask]).astype(image.dtype)

    return imout if copy else image


def ndwi(green, infrared):
    """
    McFeeters’s Normalized Difference Water Index [1]
    NDWI = (Green - NIR) / (Green + NIR)
    Positive values depict water bodies feature.
    Green and NIR values are top-of-atmosphere reflectances.
    [1] McFeeters, S.K. The use of the normalized difference water index
        (NDWI) in the delineation of open water features. Int. J. Remote
        Sens. 1996, 17, 1425–1432.
    """
    epsilon = 7/3 - 4/3 - 1  # gives machine's epsilon
    return (np.nan_to_num(green) - np.nan_to_num(infrared)) / (
            np.nan_to_num(green) + np.nan_to_num(infrared) + epsilon)


# Function by Rafael Grompone von Gioi
def rasterio_write(path, array, profile={}, tags={}):
    """
    Write a numpy array in a tiff or png file with rasterio.

    Args:
        path (str): path to the output tiff/png file
        array (numpy array): 2D or 3D array containing the image to write.
        profile (dict): rasterio profile (ie dictionary of metadata)
        tags (dict): dictionary with additional geotiff tags
    """
    # determine the driver based on the file extension
    extension = os.path.splitext(path)[1].lower()
    if extension in ['.tif', '.tiff']:
        driver = 'GTiff'
    elif extension in ['.png']:
        driver = 'png'
    else:
        raise NotImplementedError('format {} not supported'.format(extension))

    # read image size and number of bands
    array = np.atleast_3d(array)
    height, width, nbands = array.shape

    # define image metadata dict
    profile.update(driver=driver, count=nbands, width=width, height=height,
                   dtype=array.dtype)

    # write to file
    with rasterio.Env():
        with rasterio.open(path, 'w', **profile) as dst:
            dst.write(np.transpose(array, (2, 0, 1)))
            dst.update_tags(**tags)


# Function mkdir_p copyright: Carlo de Franchis
def mkdir_p(path):
    """
    Create a directory without complaining if it already exists.
    """
    if path:
        try:
            os.makedirs(path)
        except OSError as exc:  # requires Python > 2.5
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise


def load_images(inputs, option=None):
    """
    Read the input images (and their masks) using rasterio.

    Args:
        inputs

    Return:
        a first list of numpy arrays for the images
        a second list of numpy arrays for their associated masks
    """

    # The read_masks() function assumes that the 'nodata' property is set. If
    # not set, the mask is defined as valid for all pixels (no warning). In
    # case a sidecar GeoTIFF exists alongside the input file (looks like
    # <filename.ext>.msk), rasterio ignores the nodata metadata and returns the
    # mask based on the sidecar.

    im_list = []
    mask_list = []
    for impath in inputs:
        print(f'reading {impath}', end='')
        with rasterio.open(impath, 'r') as f:
            # read image
            imrgb = f.read().squeeze()
            imrgb = imrgb.astype(np.float32)
            np.nan_to_num(imrgb, copy=False)
            assert ~np.isnan(imrgb).any()
            im_list.append(imrgb)

            # read image's mask
            mask = f.read_masks().squeeze()  # TODO: is squeeze needed really?
            mask_list.append(mask)
        print('')

    return im_list, mask_list
