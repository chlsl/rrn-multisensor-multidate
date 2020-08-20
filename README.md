# Relative Radiometric Normalization (RRN)

This program normalizes the radiometric values of multi-sensor and
multi-temporal overhead image time-series.

This is the implementation for the ISPRS 2020 paper: [Relative radiometric
normalization using several automatically chosen reference images for
multi-sensor, multi-temporal
series](https://doi.org/10.5194/isprs-annals-V-2-2020-845-2020)
> Hessel, C., Grompone von Gioi, R., Morel, J. M., Facciolo, G., Arias, P., and
> de Franchis, C.: RELATIVE RADIOMETRIC NORMALIZATION USING SEVERAL
> AUTOMATICALLY CHOSEN REFERENCE IMAGES FOR MULTI-SENSOR, MULTI-TEMPORAL
> SERIES, ISPRS Ann.  Photogramm. Remote Sens.  Spatial Inf. Sci., V-2-2020,
> 845–852, https://doi.org/10.5194/isprs-annals-V-2-2020-845-2020, 2020.

Version 0.0.1, released on August 20, 2020.

By [Charles Hessel](mailto:charles.hessel@ens-paris-saclay.fr), CMLA, ENS
Paris-Saclay.

The code contains contributions from Rafael Grompone von Gioi, Tristan
Dagobert, Carlo de Franchis.

Main source code repository for this software:
https://github.com/chlsl/rrn-multisensor-multidate


## Installation

The simplest way to install `rrn` is to run
```
pip install .
```
from the directory in which the `setup.py` file is located.

This will compile the two C modules that are included in the project, install
`rrn` as a package available in your python projects, and also make the
command-line interface `rrn` available from anywhere on your machine.

This project uses `libtiff`, `libpng` and `libjpeg`.

The required python packages will be installed by pip. They are listed in the
`requirements.txt` file.

#### Testing

Execute the script `tests/test.sh` to verify your installation. This should
create a directory `out`, which contains the output normalized time series.


## Usage

The command
```
rrn sequence-registered/*{B02,B02,B04}.tif --band B02 B03 B04 --outdir sequence-registered-stabilized
```
will take all images given in argument and save their stabilized version in the
directory "sequence-registered-stabilized". The stabilized images are floating
point tif RGB, and to their names is appended `_stab`. Tone-mapped versions of
the same images are saved too; they are suffixed `_stab_tonemap` and saved as
png.  The input images can be either RGB or separate bands. Their type must be
specified with the argument `--band` or `-b`.


#### Input images

The script takes as input a registered time series.

For RNN to correctly load the input images, their name must contain indications
concerning the band and the satellite (and correction level for S2).
  - Concerning the bands, the filenames must contain `_BBB`, where `BBB` stands
    for B02, B03, B04 or RGB. The `B02` band is supposed to be the blue, `B03`
    the green and `B04` the red. You can also give a 3-channels image with the
    3 bands already merged, in which case the input names must contain `_RGB`.
    Rename your images if needed.
  - Likewise, the satellite and correction level, are collected from the filenames.
    - For Sentinel 2 we're looking for `_S2A_` and `_S2B_`; and `_L1C_` or
      `_L2A_` for the correction level.
    - For Landsat 8 we look for `_L8_`.

**Example:**
The image `2019-07-03_S2B_orbit_094_tile_31UDQ_L2A_band_B02.tif` is detected as
a Sentinel-2 image, correction level 2A and band B02 (blue).
You can see more examples in the test directory.


#### Temporal consistency

RRN uses several key images so as to handle seasonal color variations.
To obtain more stable sequences, use the `--local-max` option.
Setting it to a higher value will cause RRN to use fewer key images and
therefore to allow less color change. If you set it to the length of the input
sequence or above, only one key image will be used.


#### Overlay

To remove the overlay with the date that appears on the bottom of the stabilized
images, use the flag `--no-overlay`.


#### Options

Use `rrn --help` to get help and a list of available options (reported below).

```
Relative Radiometric Normalization

optional arguments:
  -h, --help            show this help message and exit

positional arguments (mandatory):
  imlist                list of images

mandatory arguments:
  -b {RGB,B08,B03,B04,B02} [{RGB,B08,B03,B04,B02} ...], --band {RGB,B08,B03,B04,B02} [{RGB,B08,B03,B04,B02} ...]
                        space-separated list of bands to load (in the provided
                        files, see "imlist"). Specifying RGB implies that the
                        B02, B03 and B04 bands were merged beforehand. Default
                        is RGB.
  -o OUTDIR, --outdir OUTDIR
                        output directory (created if nonexistent)

parameters of the method:
  -m LOCAL_MAX, --local-max LOCAL_MAX
                        a notion of how "local" the local maxima are. Default
                        is 9.
  --no-visible-ground-area NO_VISIBLE_GROUND_AREA
                        threshold on the proportion (percentage) of not visible
                        ground pixels. Images whose proportion of not visible
                        ground pixels is above this threshold are rejected.
                        Default is 25.

by-pass some of the method's steps:
  --key-images KEY_IMAGES [KEY_IMAGES ...]
                        list of index for key images.
  -c CLOUDS [CLOUDS ...], --clouds CLOUDS [CLOUDS ...]
                        input cloud masks. If not provided, they are computed
                        using R. Grompone von Gioi ground visibility detector.
  --load-previous       use this flag to load previously computed intermediary
                        results rather than computing them. Useful for
                        debugging mainly.

alternatives for some steps:
  --acoefs-method {ransac,L1constrained,L1,median}
                        computation method for the affine coefficients. Default
                        is ransac.
  --usecorrsift         flag to replace corrsift (slow) by the angles (faster).
                        The angle error is less accurate. Default is True.
  --dont-usecorrsift    set to False the usecorrsift flag.

more options!:
  --verbose             print non-critical information
  --no-overlay          use this flag to remove the bottom overlay with the
                        date added on the stabilized images.
  --roughnorm           flag to apply a rough normalization to the input
                        images. Default is True.
  --no-roughnorm        Set to False the roughnorm flag.
  --tonemap             flag to apply a tonemapping step on the output images.
                        Default is True. (use --no-tonemap to deactivate)
  --no-tonemap          set to False the tonemap flag.
```


## Licence

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.

