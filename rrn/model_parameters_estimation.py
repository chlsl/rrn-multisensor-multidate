import iio
import numpy as np
from scipy.optimize import minimize


def affine_coefficients(target, image, mask, method=None, stdnoise=None):
    assert target.ndim < 3 and image.ndim < 3, 'vector or 1-channel image'
    if method is None:
        method = 'ransac'
    if method == 'L1':
        a, b = affine_coefficients_opt(target, image, mask, force_a_ge_1=False)
        return a, b, None
    elif method == 'L1constrained':
        a, b = affine_coefficients_opt(target, image, mask, force_a_ge_1=True)
        return a, b, None
    elif method == 'median':
        a, b = affine_coefficients_med(target, image, mask)
        return a, b, None
    elif method == 'ransac':
        threshold = 20 * stdnoise
        a, b, inliers = ransac(image[mask], target[mask], threshold)
        return a, b, inliers
    else:
        raise ValueError('Incorrect acoefs_method.')


def affine_coefficients_opt(target, image, mask, force_a_ge_1=True):
    """
    Compute the affine coefficients (a, b) so that a x image + b is as close as
    possible to the target image (for the pixels in the provided mask).
    To avoid contrast compression, it can be usefull to force a >= 1.

    :param target:
    :param image:
    :param mask:
    :param force_a_ge_1:
    :return: affine coefficients (a, b) so that a * image(p) + b - target(p) is
             minimal for pixels p in the mask.

    This function is inspired by Axel's function "correct_to_target_images"
    """
    assert target.ndim == 2, 'target image must be gray'
    assert image.ndim == 2, 'input image must be gray'
    assert target.shape == image.shape

    # make sure the data is np.float64 (with float32 the method 'L-BFGD-B' does
    # not converge).
    target = target.astype(np.float64)  # float64 is required
    image = image.astype(np.float64)  # not sure flatten is usefull
    if force_a_ge_1:
        res = minimize(
            lambda w: np.sum(np.abs((target - w[0]*image - w[1])[mask])),
            x0=np.array([1, 0]), method='L-BFGS-B',
            bounds=((1, None), (None, None)))
    else:
        res = minimize(
            lambda w: (np.abs((target - w[0]*image - w[1])[mask])).sum(),
            np.array([1, 0]), method='Powell')
    a, b = res.x
    return a, b


def affine_coefficients_med(target, image, mask):
    """Compute affine coefficients using median and mad
    The coefficients (a, b) are such that
    (a * image + b - target) is minimal.
    In this function the median and median of absolute deviation to the median
    are used.
    """
    assert target.ndim == 2, 'target image must be gray'
    assert image.ndim == 2, 'input image must be gray'
    assert target.shape == image.shape

    med_target = np.median(target[mask])
    mad_target = np.median(np.abs(target[mask] - med_target))
    med_image = np.median(image[mask])
    mad_image = np.median(np.abs(image[mask] - med_image))

    a = mad_target / mad_image
    b = med_target - a * med_image
    return a, b


def ransac(source, target, threshold, tolerance=None, hourglass=2048):
    """RANSAC Random Sample Consensus for image values transformations

    Args:
        source (1D ndarray): source image (vectorized)
        target (1D ndarray): target image (vectorized)
        threshold (float): maximal (orthogonal) distance to line to count the
                point as inlier. Perhaps the most important parameter.
        tolerance (float): stopping criterion: if error of model is below, it
                is considered as sufficient and the algorithm stops.
                Default is threshold / 3
        hourglass (int): maximal number of iterations. Default is 1024

    Return:
        model parameters ([a, b] so that a * source + b ~= target), inliers set
    """

    assert source.ndim == 1 and target.ndim == 1
    assert target.shape[0] == target.shape[0]

    if not tolerance: tolerance = threshold / 3
    assert tolerance < threshold

    initial_size = source.size  # store this for later

    def reduce_points_density(img):
        assert img.size > 0

        permin = np.percentile(img, .1)
        permax = np.percentile(img, 99.9)

        # values before permin will be assigned to bin 0
        # values after permax will be assigned to bin 99
        ibins = np.linspace(permin, permax, 99)
        iquant = np.digitize(img, ibins)

        # number of values in bins if all equals -- multiplied for threshold
        thresh = np.round(3 * source.size / 100).astype(np.int)

        # a non-standard histogram, which contains the position of the pixels;
        # when no element, empty array
        histog = [(iquant == b).nonzero()[0] for b in range(101)]
        histog = [h for h in histog if h.size > 0]  # remove empty arrays

        hiscut = []
        for h in histog:
            if h.size > thresh:
                np.random.shuffle(h)
                hiscut.append(h[:thresh])  # take only the first elmts
            else:
                hiscut.append(h)

        return [i for h in hiscut for i in h]

    # remove pixels at extremities
    permin = np.percentile(source, .1)
    permax = np.percentile(source, 99.99)
    if permax <= permin:
        print('WARNING, permax < permin! permax = {permax}, permin = {permin}')
        permax = +np.inf
        permin = 0
    idx_init = np.nonzero(np.logical_and(source > permin, source < permax))[0]
    permin = np.percentile(target, .1)
    permax = np.percentile(target, 99.99)
    if permax <= permin:
        print('WARNING, permax < permin! permax = {permax}, permin = {permin}')
        permax = +np.inf
        permin = 0
    idx_init_2 = np.nonzero(np.logical_and(target[idx_init] > permin,
        target[idx_init] < permax))[0]
    idx_init = idx_init[idx_init_2]
    source = source[idx_init]
    target = target[idx_init]

    # remove pixel in too dense areas, source
    idx_source = reduce_points_density(source)
    source, target = source[idx_source], target[idx_source]

    # remove pixel in too dense areas, target
    idx_target = reduce_points_density(target)
    source, target = source[idx_target], target[idx_target]

    # print(f'Running RANSAC with threshold set to {threshold:.4f}; '
    #         f'Removed {100 * (initial_size - source.size) / initial_size:.1f}% '
    #         f'of the {initial_size} pixels in the source and target images.')

    # define a function for the regression using all found inliers
    def total_least_squares(x, y, d=1):
        """Total Least Squares, an alternative to the ordinary least squares.
        See also "Deming Regression": en.wikipedia.org/wiki/Deming_regression
        """

        assert type(x) is np.ndarray and type(y) is np.ndarray
        assert x.shape == y.shape and x.ndim == 1

        x_m, y_m = np.mean(x), np.mean(y)
        xx_v, yy_v = np.mean((x - x_m)**2), np.mean((y - y_m)**2)
        xy_c = np.mean((x - x_m) * (y - y_m))

        a = (yy_v - d * xx_v + np.sqrt((yy_v - d * xx_v)**2 + 4 * d * xy_c**2)
                ) / (2 * xy_c)
        b = y_m - a * x_m

        return a, b


    indices = np.arange(source.shape[0], dtype=np.int)
    best_model = (1, 0)  # best model
    best_nbinliers = 0  # best number of inliers
    best_error = +np.inf  # best error (smallest one)
    best_inliers = None  # best inliers (bool array)

    while best_error > tolerance:

        # estimate model parameters from two random points
        np.random.shuffle(indices)
        x = source[indices[1]] - source[indices[0]]
        y = target[indices[1]] - target[indices[0]]
        # if x is zero, we can't continue (division by zero)
        # if y is zero, we already now that the model is very very bad. Skip it.
        if x == 0 or y == 0: continue
        a = y / x
        b = target[indices[0]] - a * source[indices[0]]

        # distance of the points to the line (shortest distance)
        distances = np.abs(target - a * source - b) / np.sqrt(a**2 + 1)

        # for the model estimated, find all other points that match (= inliers)
        inliers = distances < threshold
        nbinliers = np.sum(inliers)  # number of inliers

        # If the number of inliers is sufficiently large, use all the inliers
        # to make a better estimation of the model's parameters
        if nbinliers > best_nbinliers:

            while hourglass > 0:

                aa, bb = total_least_squares(source[inliers], target[inliers])
                distances_refined = np.abs(target - aa * source - bb) \
                                    / np.sqrt(aa**2 + 1)
                inliers_refined = distances_refined < threshold
                nbinliers_refined = np.sum(inliers_refined)

                # if the nb of inliers after refinement is larger than
                # nbinliers, then this is an even better model
                if nbinliers_refined > nbinliers:

                    a, b = aa, bb
                    distances = distances_refined
                    inliers = inliers_refined
                    nbinliers = nbinliers_refined
                    hourglass -= 1

                    # average error for the inliers
                    error = np.mean(distances[inliers])

                    best_nbinliers = nbinliers
                    best_model = (a, b)
                    best_error = error
                    best_inliers = inliers

                else:
                    break

        hourglass -= 1

        if hourglass <= 0:
            break

    best_inliers_full = np.zeros((initial_size,), dtype=np.bool)
    idx_init = np.array(idx_init)
    idx_source = np.array(idx_source)
    idx_target = np.array(idx_target)
    best_inliers = np.array(best_inliers)
    best_inliers_full[idx_init[idx_source[idx_target[best_inliers]]]] = True

    return (*best_model, best_inliers_full)

