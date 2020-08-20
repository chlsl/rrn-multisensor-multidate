
/*
* nimbus.cpp
* 
* Copyright (C) 2016, Tristan Dagobert, CMLA, École Normale Supérieure de Cachan.
* 
* This software is a computer program.[describe
* functionalities and technical features of your software].
* 
* This software is governed by the CeCILL-C license under French law and
* abiding by the rules of distribution of free software.  You can  use, 
* modify and/ or redistribute the software under the terms of the CeCILL-C
* license as circulated by CEA, CNRS and INRIA at the following URL
* "http://www.cecill.info". 
* 
* As a counterpart to the access to the source code and  rights to copy,
* modify and redistribute granted by the license, users are provided only
* with a limited warranty  and the software's author,  the holder of the
* economic rights,  and the successive licensors  have only  limited
* liability. 
* 
* In this respect, the user's attention is drawn to the risks associated
* with loading,  using,  modifying and/or developing or reproducing the
* software by the user in light of its specific status of free software,
* that may mean  that it is complicated to manipulate,  and  that  also
* therefore means  that it is reserved for developers  and  experienced
* professionals having in-depth computer knowledge. Users are therefore
* encouraged to load and test the software's suitability as regards their
* requirements in conditions enabling the security of their systems and/or 
* data to be ensured and,  more generally, to use and operate it in the 
* same conditions as regards security. 
* 
* The fact that you are presently reading this means that you have had
* knowledge of the CeCILL-C license and that you accept its terms.
* 
*/
    
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// extern "C"
// {
#include "iio.h"
#include "lib_sift.h"
#include "bilinear_interpolation.h"
// }

#include "corrsift.h"
// #include "gaussian.h"
// #include "nfa_multi_uniform.h"

#define ETA 1e-6

#define SAMPLE_FACTOR 10
#define HISTO_NB_CLASSES 2000

#define SIFT_NB_CLASSES 128

/*****************************************************************************/
float * 
extract_partition(int nl, int nc, float * im, int ncol, int dl, int dc)
		  
{
  int l, c;
  float * part = (float *) malloc(nl * nc * sizeof(float));

  for (l=0; l<nl; l++) 
    {
      for (c=0; c<nc; c++) 
	{
	  part[l*nc + c] = im[(l + dl) * ncol + (c + dc)];
	}
    }
  return part;
}
/*****************************************************************************/
void 
merge_partition(float * im, int ncol, float * part, int nl, int nc,
		int dl, int dc)
{
  int l, c;
  for (l=0; l < nl; l++) 
    {
      for (c=0; c < nc; c++) 
	{
	  im[(l+ dl) * ncol + (c+ dc)] = part[l* nc + c];
	}
    }

  return;
}
/*****************************************************************************/
float * zoom_sift(float * im, int ncol, int nlig, int factor, int zoom_ncol,
		  int zoom_nlig)
{
  /*     x0−…−x1−…−x2−−x3
   *     |         |   |
   *     …         …   …
   *     |   a     | b |
   *     x4−…−x5−…−x6−−x7
   *     |   c     | d |
   *     x8−…−x9−…−xA−−xB
   *
   * L'image à zoomer a été échantillonnée sur les «x» régulièrement, sauf
   * aux extrémités si le facteur d'échantillonnage n'était pas un diviseur de
   * la taille de l'image. Il existe 4 partitions possibles de l'image à 
   * ré-échantillonner: {a}, {a, b}, {a, c} et {a, b, c, d}. Il faut donc 
   * interpoler adéquatement.
   */

  int nl, nc;
  int dl, dc; /* biais */

  int a_ncol, a_nlig;
  int A_ncol, A_nlig;
  int B_ncol, B_nlig;
  int C_ncol, C_nlig;
  int D_ncol, D_nlig;

  a_ncol = (zoom_ncol%factor == 0) ? zoom_ncol/factor : zoom_ncol/factor + 1;
  a_nlig = (zoom_nlig%factor == 0) ? zoom_nlig/factor : zoom_nlig/factor + 1;

  /* agrégation des partitions {A, B, C, D} */
  float * neo;
  neo = (float *) malloc(zoom_nlig * zoom_ncol * sizeof(float));
 
  /* traitement partition a */
  float * part_a = NULL, * part_A = NULL;
  nl = a_nlig;
  nc = a_ncol;
  /* copie de {x0, …, x1, …, x2, …, …, …, x4, …, x5, …, x6} */
  dl = 0;
  dc = 0;
  part_a = extract_partition(nl, nc, im, ncol, dl, dc);
  /* zoom de {a} −−> {A} */

  A_ncol = factor * (a_ncol-1) + 1;
  A_nlig = factor * (a_nlig-1) + 1;
  part_A = (float *) malloc(A_nlig * A_ncol * sizeof(float));
  bilinear_interpolant_vec(part_A, A_ncol, A_nlig, part_a, nc, nl, 1);

  /* copie de {A} */
  dl = 0;
  dc = 0;
  merge_partition(neo, zoom_ncol, part_A, A_nlig, A_ncol, dl, dc);

  /* traitement partition b */
  float * part_b = NULL, * part_B = NULL; 
  nl = a_nlig;
  nc = ncol - a_ncol + 1;
  if (nc > 1)
    {
      /* copie de {x2, x3, …, …, x6, x7} */
      dl = 0;
      dc = a_ncol - 1;
      part_b = extract_partition(nl, nc, im, ncol, dl, dc);
      /* zoom de {b} −−> {B} */

      B_nlig = A_nlig;
      B_ncol = zoom_ncol - A_ncol + 1;
      part_B = (float *) malloc(B_nlig * B_ncol * sizeof(float));
      bilinear_interpolant_vec(part_B, B_ncol, B_nlig, part_b, nc, nl, 1);

      /* copie de {B} */
      dl = 0;
      dc = A_ncol-1;
      merge_partition(neo, zoom_ncol, part_B, B_nlig, B_ncol, dl, dc);

    }

  /* traitement partition c */
  float * part_c = NULL, * part_C = NULL; 
  nl = nlig - a_nlig + 1;
  nc = a_ncol;
  if (nl > 1)
    {
      /* copie de {x4, …, x5, …, x6, …, …, …, x8, …, x9, …, xA} */
      dl = a_nlig - 1;
      dc = 0;
      part_c = extract_partition(nl, nc, im, ncol, dl, dc);
      /* zoom de {c} −−> {C} */
      C_nlig = zoom_nlig - A_nlig + 1;
      C_ncol = A_ncol;
      part_C = (float *) malloc(C_nlig * C_ncol * sizeof(float));
      bilinear_interpolant_vec(part_C, C_ncol, C_nlig, part_c, nc, nl, 1);

      /* copie de {C} */
      dl = A_nlig-1;
      dc = 0;
      merge_partition(neo, zoom_ncol, part_C, C_nlig, C_ncol, dl, dc);

    }

  /* traitement partition d */
  float * part_d = NULL, * part_D = NULL; 
  nl = nlig - a_nlig + 1;
  nc = ncol - a_ncol + 1;
  if (nl > 1 && nc > 1)
    {
      /* copie de {x6, x7, xA, xB} */
      dl = a_nlig - 1;
      dc = a_ncol - 1;
      part_d = extract_partition(nl, nc, im, ncol, dl, dc);
      /* zoom de {d} −−> {D} */
      D_nlig = zoom_nlig - A_nlig + 1;
      D_ncol = zoom_ncol - A_ncol + 1;
      part_D = (float *) malloc(D_nlig * D_ncol * sizeof(float));
      bilinear_interpolant_vec(part_D, D_ncol, D_nlig, part_d, nc, nl, 1);

      /* copie de {D} */
      dl = A_nlig-1;
      dc = A_ncol-1;
      merge_partition(neo, zoom_ncol, part_D, D_nlig, D_ncol, dl, dc);

    }

  free(part_a);
  free(part_b);
  free(part_c);
  free(part_d);
  free(part_A);
  free(part_B);
  free(part_C);
  free(part_D);
  return neo;
}
/*****************************************************************************/
sift_point_t * 
compute_median_descriptors(sift_point_t ** pts_sift, int nb_pts_clef, int nb_ims)
{
  int i, j;
  sift_point_t * pts_sift_ref;
  pts_sift_ref = (sift_point_t *) malloc(nb_pts_clef * sizeof(sift_point_t));

  unsigned char ** sift_descriptors = 
    (unsigned char **) malloc(nb_ims * sizeof(unsigned char *));
  for (i = 0; i < nb_pts_clef; i++) 
    {
      for (j = 0; j < nb_ims; j++) 
	{
	  sift_descriptors[j] = pts_sift[j][i].descriptor;
	}
      j = compute_the_median_descriptor(sift_descriptors, nb_ims, SIFT_NB_CLASSES);
      memcpy(pts_sift_ref[i].descriptor, sift_descriptors[j],  
	     SIFT_NB_CLASSES * sizeof(unsigned char));
    }

  free(sift_descriptors);
  return pts_sift_ref;
}
/*****************************************************************************/
int 
compute_the_median_descriptor(unsigned char ** descriptors, int nb_desc, int nb_classes)
{
  int i, j, c;
  float d, d2;
  float sum_of_distances[nb_desc];
  
  memset(sum_of_distances, 0, nb_desc * sizeof(float));

  /* algorithme de Weiszfeld (1938)… */
  for (i = 0; i < nb_desc; i++) 
    {
      for (j = i+1; j < nb_desc; j++) 
	{
	  for (d2=0.f, c = 0; c < nb_classes; c++) 
	    {
	      d = (descriptors[i][c] - descriptors[j][c]);
	      d2 += d * d; 
	    }
	  sum_of_distances[i] += d2;
	  sum_of_distances[j] += d2;
	}
    }

  /* recherche du minimum */
  float mini = sum_of_distances[0];
  int idx_mini = 0;

  for (i = 0; i < nb_desc; i++) 
    {
      if (sum_of_distances[i] < mini)
	{
	  mini = sum_of_distances[i];
	  idx_mini = i;
	}
    }
  
  return idx_mini;
}

/*****************************************************************************/
sift_point_t * 
create_sift_descriptors(float * im, int nlig, int ncol,	int * nb_pts_clef, 
			int * nb_pts_lig, int * nb_pts_col)
{
  int i, j, l, c;
  long int p;

  *nb_pts_col = (ncol%SAMPLE_FACTOR == 0) ? ncol/SAMPLE_FACTOR + 1: ncol/SAMPLE_FACTOR + 2;
  *nb_pts_lig = (nlig%SAMPLE_FACTOR == 0) ? nlig/SAMPLE_FACTOR + 1: nlig/SAMPLE_FACTOR + 2;
  *nb_pts_clef = *nb_pts_col * *nb_pts_lig;

  // TODO: vérifier pour le dernier point…

  sift_point_t * pts_sift;

  pts_sift = (sift_point_t *) malloc(*nb_pts_clef * sizeof(sift_point_t));
  
  /* initialisation : les points SITF sont calculés de façon régulièrement
   * espacée en «x», sauf pour les derniers points des lignes et colonnes :
  /*                     x−−−x−−−x−x
   *                     |       | |
   *                     |   a   |b|
   *                     x−−−x−−−x−x
   *                     |   c   |d|
   *                     x−−−x−−−x−x
   */

  for (p=0, l=0, i=0; i<*nb_pts_lig; i++, l = (l+SAMPLE_FACTOR < nlig-1) ? l+SAMPLE_FACTOR : nlig-1)
    {
      for (c=0, j=0; j<*nb_pts_col; j++, p++, c = (c+SAMPLE_FACTOR < ncol-1) ? c+SAMPLE_FACTOR : ncol-1)
	{
	  //	  printf("(%d %d) ", l, c);
	  pts_sift[p].x = l;
	  pts_sift[p].y = c;
	  pts_sift[p].scale = 1;
	  pts_sift[p].orientation = 0;
	}
      //printf("\n");
    }  

  /* calcul des descripteurs SIFT de l'image à tester */
  sift_fill_descriptors(im, ncol, nlig, pts_sift, *nb_pts_clef);

  return pts_sift;
}
/*****************************************************************************/
float * 
compute_sift_correlation(sift_point_t * pts_sift, sift_point_t * pts_sift_ref,
			 int nb_pts_lig, int nb_pts_col, int nlig, int ncol)
{
  int FAC = SAMPLE_FACTOR;
  const int nb_classes = 128;
  int i, j, k, l, c;
  long int p;

  /* calcul de la corrélation entre les descripteurs de référence et ceux de 
   * l'image à tester */
  unsigned char * histo, * histo_ref;
  float x, y, sxy=0, sx=0, sy=0, sx2=0, sy2=0;

  float * neo;
  neo = (float *) calloc(nlig * ncol, sizeof(float));

  float * neo_p;
  neo_p = (float *) calloc(nb_pts_lig * nb_pts_col, sizeof(float));

  // printf("neo %d %d\n", nlig, ncol); fflush(stdout);
  for (p=0, l=0, k=0; k<nb_pts_lig; k++, l+=FAC)
    {
      for (c=0, j=0; j<nb_pts_col; j++, c+=FAC, p++)
	{
	  histo = pts_sift[p].descriptor;
	  histo_ref = pts_sift_ref[p].descriptor;

	  for (sxy=0, sx=0, sy=0, sx2=0, sy2=0, i=0; i<nb_classes; i++)
	    {
	      x = histo[i];
	      y = histo_ref[i];
	      sxy += x * y;
	      sx += x;
	      sy += y;
	      sx2 += x * x;
	      sy2 += y * y;
	    }

	  neo_p[p] = (nb_classes * sxy - sx * sy) / 
	    (sqrt(nb_classes * sx2 - sx * sx) * sqrt(nb_classes * sy2 - sy * sy));
	}
    }

  /* interpolation intelligente… */
  neo = zoom_sift(neo_p, nb_pts_col, nb_pts_lig, FAC, ncol, nlig);

  free(neo_p);
  return neo;
}

/*****************************************************************************/

