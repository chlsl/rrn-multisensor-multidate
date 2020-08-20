/* Interface for python, using ctypes module
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "corrsift.h"
// #include "iio.h"

void corrsift(float *im0, float *im1, float *cmap, int im_h, int im_w);

void corrsift(float *im0, float *im1, float *cmap, int im_h, int im_w) {

    int nlig=im_h, ncol=im_w;
    int i;

    int nb_ims=2; 
    int nb_pts_sift, nb_pts_lig, nb_pts_col;
    sift_point_t ** pts_sift=NULL;
    sift_point_t * pts_sift_ref=NULL;
    float ** corr_ims_sift=NULL;
    float ** ims_sift=NULL;

    // nb_imgs is always 2
    ims_sift = (float**) malloc(nb_ims * sizeof(float*));
    ims_sift[0] = im0;
    ims_sift[1] = im1;

    /* calcul des descripteurs SIFT */
    pts_sift = (sift_point_t **) malloc(nb_ims * sizeof(sift_point_t*));
    for (i = 0; i < nb_ims; i++) 
    {
        // printf("ims_sift[%ld]=%p\n", i, &(ims_sift[i]));
        pts_sift[i] = create_sift_descriptors(ims_sift[i], nlig, ncol,
                &nb_pts_sift, &nb_pts_lig, &nb_pts_col);
    }
    // printf("[OK] create_sift_descriptors pts_ligne=%d pts_col=%d\n",
    //         nb_pts_lig, nb_pts_col);
    // fflush(stdout);

    /* calcul des descripteurs médians */
    pts_sift_ref = pts_sift[0];  // I use the first image as reference, always.

    /* calcul des corrélations */
    corr_ims_sift = (float **) malloc(nb_ims * sizeof(float*));
    for (i = 0; i < nb_ims; i++) 
    {
        corr_ims_sift[i] = compute_sift_correlation(pts_sift[i], pts_sift_ref,
                nb_pts_lig, nb_pts_col, nlig, ncol);
        // sprintf(nom, "correla_sift_%02d.tif", i);
        // iio_save_image_float((char *)nom, corr_ims_sift[i], ncol, nlig);
        if (i == 1) {
            for (int j = 0; j < ncol * nlig; j++) {
                cmap[j] = corr_ims_sift[i][j];
            }
        }
    }
    // iio_save_image_float("iio_out.tif", cmap, im_w, im_h);

    for (i = 0; i < nb_ims; i++) {
        free(corr_ims_sift[i]);
        free(pts_sift[i]);
    }
    free(corr_ims_sift);
    free(pts_sift);
    // free(pts_sift_ref);
}
