
/*
* nimbus.h
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
    
typedef struct sift_keypoint_std sift_point_t;

float * zoom_sift(float * im, int ncol, int nlig, int factor, int zoom_ncol,
		  int zoom_nlig);
sift_point_t * 
compute_median_descriptors(sift_point_t ** pts_sift, int nb_pts_clef, int nb_ims);
int 
compute_the_median_descriptor(unsigned char ** descriptors, int nb_desc, int nb_classes);
sift_point_t * 
create_sift_descriptors(float * im, int nlig, int ncol,	int * nb_pts_clef, 
			int * nb_pts_lig, int * nb_pts_col);
float * 
extract_partition(int nl, int nc, float * im, int ncol, int dl, int dc);
void 
merge_partition(float * im, int ncol, float * part, int nl, int nc,
		int dl, int dc);
float * 
compute_sift_correlation(sift_point_t * pts_sift, sift_point_t * pts_sift_ref,
			 int nb_pts_lig, int nb_pts_col, int nlig, int ncol);



