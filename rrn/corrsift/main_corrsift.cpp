#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern "C"
{
#include "iio.h"
}

#include "corrsift.h"

#define VERSION "1.0"
#define NB_ITERATIONS 1
/*****************************************************************************/
    void 
traiter_commande(int argc, char **argv, 
        char *** fics_sift, int * nb_fics_sift)
{
    int i;

    if (argc == 1)
    {
        printf("corrsift.exe %s\n", VERSION);

        printf("%s S(1)… S(n)\n", argv[0]);
        exit(EXIT_SUCCESS);
    }
    *nb_fics_sift = argc-1;
    *fics_sift = (char **) malloc(*nb_fics_sift * sizeof(char *));

    for (i=1; i<argc; i++)
        (*fics_sift)[i-1] = argv[i];

    return;
}
/*****************************************************************************/
int main(int argc, char *argv[])
{
    int nlig, ncol;
    int i;

    char ** fics_sift; 
    int nb_ims=0; 

    char nom[60];

    int nb_pts_sift, nb_pts_lig, nb_pts_col;
    sift_point_t ** pts_sift=NULL;
    sift_point_t * pts_sift_ref=NULL;
    float ** corr_ims_sift=NULL;

    traiter_commande(argc, argv, 
            &fics_sift, &nb_ims);

    printf("Fic  SIFT %d:\n", nb_ims);
    for (i = 0; i < nb_ims; i++) 
    {
        printf("%s\n", fics_sift[i]);
    }

    float ** ims_sift=NULL;

    /* [1] lecture des descripteurs SIFT */
    if (nb_ims > 0)
    {
        ims_sift = (float **) malloc(nb_ims * sizeof(float *));
        for (i = 0; i < nb_ims; i++)
        {
            ims_sift[i] =  iio_read_image_float(fics_sift[i], &ncol, &nlig); 
        }
    }

    const long int npix = nlig * ncol;

    /* calcul des descripteurs SIFT */
    pts_sift = (sift_point_t **) malloc(nb_ims * sizeof(sift_point_t*));
    for (i = 0; i < nb_ims; i++) 
    {
        printf("ims_sift[%ld]=%p\n", i, &(ims_sift[i]));
        pts_sift[i] = create_sift_descriptors(ims_sift[i], nlig, ncol,
                &nb_pts_sift, &nb_pts_lig, &nb_pts_col);
    }
    printf("[OK] create_sift_descriptors pts_ligne=%d pts_col=%d\n", nb_pts_lig, nb_pts_col);
    fflush(stdout);

    /* calcul des descripteurs médians */
    pts_sift_ref = compute_median_descriptors(pts_sift, nb_pts_sift, nb_ims);
    printf("[OK] compute_median_descriptors\n");
    fflush(stdout);

    /* calcul des corrélations */
    corr_ims_sift = (float **) malloc(nb_ims * sizeof(float*));
    for (i = 0; i < nb_ims; i++) 
    {
        corr_ims_sift[i] = compute_sift_correlation(pts_sift[i], pts_sift_ref,
                nb_pts_lig, nb_pts_col, nlig, ncol);
        sprintf(nom, "correla_sift_%02d.tif", i);
        iio_save_image_float((char *)nom, corr_ims_sift[i], ncol, nlig);
        free(corr_ims_sift[i]);
        free(pts_sift[i]);
    }

    free(corr_ims_sift);
    free(pts_sift);
    free(pts_sift_ref);
}
