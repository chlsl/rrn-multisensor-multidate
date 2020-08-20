// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

// #include <cmath>
// #include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gaussian.h"

/**
 *
 * Convolution with a Gaussian
 *
 */
 void gaussian1D(
    float *I,             //input/output image
    const int xdim,       //image width
    const int ydim,       //image height
    const double sigma,   //Gaussian sigma
    const int bc,       //boundary condition
    const int precision //defines the size of the window
)
{
    int i, j, k;

    const double den  = 2*sigma*sigma;
    const int   size = (int) (precision * sigma) + 1 ;
    const int   bdx  = xdim + size;
    const int   bdy  = ydim + size;

    if ( bc && size > xdim ) {
    printf("GaussianSmooth: sigma too large for this bc\n");
	// std::cerr << "GaussianSmooth: sigma too large for this bc\n" << std::endl;
	// throw 1;
    exit(1);
    }

    // compute the coefficients of the 1D convolution kernel
    // double *B = new double[size];
    double *B = (double *) malloc(size * sizeof(double));
    for(int i = 0; i < size; i++)
	B[i] = 1 / (sigma * sqrt(2.0 * 3.1415926)) * exp(-i * i / den);

    double norm = 0;

    // normalize the 1D convolution kernel
    for(int i = 0; i < size; i++)
	norm += B[i];

    norm *= 2;

    norm -= B[0];

    for(int i = 0; i < size; i++)
	B[i] /= norm;

//BBB    // convolution of each line of the input image
//BBB    double *R = new double[size + xdim + size];
//BBB    for ( k = 0; k < ydim; k ++ )
//BBB    {
//BBB	for (i = size; i < bdx; i++)
//BBB	    R[i] = I[k * xdim + i - size];
//BBB
//BBB	switch ( bc )
//BBB	{
//BBB	    case 0:   // Dirichlet boundary conditions
//BBB
//BBB		for( i = 0, j = bdx; i < size; i++,j++) R[i] = R[j] = 0; break;
//BBB
//BBB	    case 1:   // Reflecting boundary conditions
//BBB
//BBB		for(i = 0, j = bdx; i < size; i++, j++) 
//BBB		{
//BBB		    R[i] = I[k * xdim + size-i];
//BBB		    R[j] = I[k * xdim + xdim-i-1];
//BBB		}
//BBB		break;
//BBB
//BBB	    case 2:   // Periodic boundary conditions
//BBB
//BBB		for(i=0,j=bdx;i<size;i++,j++) 
//BBB		{
//BBB		  R[i] = I[k * xdim + xdim-size+i];
//BBB		  R[j] = I[k * xdim + i];
//BBB		}
//BBB		break;
//BBB	}
//BBB
//BBB	for ( i=size; i<bdx; i++ )
//BBB	{
//BBB
//BBB	    double sum = B[0] * R[i];
//BBB
//BBB	    for ( int j = 1; j < size; j ++ )
//BBB	      sum += B[j] * ( R[i-j] + R[i+j] );
//BBB
//BBB	    I[k * xdim + i - size] = sum;
//BBB	}
//BBB    }
//BBB
    // convolution of each column of the input image
    // double *T = new double[size + ydim + size];
    double *T = (double *) malloc((size + ydim + size) * sizeof(double));
    for ( k = 0; k < xdim; k ++ )
    {
	for ( i=size; i<bdy; i++ ) T[i] = I[(i - size) * xdim + k];

	switch ( bc )
	{
	    case 0:   // Dirichlet boundary conditions

		for ( i = 0, j = bdy; i < size; i++, j++) T[i] = T[j] = 0; break;

	    case 1:   // Reflecting boundary conditions

		for(i=0,j=bdy;i<size;i++,j++) 
		{
		    T[i] = I[(size-i) * xdim + k];
		    T[j] = I[(ydim-i-1) * xdim + k];
		}
		break;

	    case 2:   // Periodic boundary conditions

		for(i=0,j=bdx;i<size;i++,j++) 
		{
		    T[i] = I[(ydim-size+i) * xdim + k];
		    T[j] = I[i * xdim + k];
		}
		break;
	}

	for ( i=size;i<bdy; i++ )
	{
	    double sum = B[0] * T[i];

	    for ( j = 1; j < size; j ++ )
		sum += B[j] * (T[i-j] + T[i+j]);

	    I[(i - size) * xdim + k] = sum;
	}
    }

    free(B);
    free(T);
    // delete []B;
    // //    delete []R;
    // delete []T;
}




/**
 *
 * Convolution with a Gaussian
 *
 */
void smooth_with_gaussian(
    const float *in,      //input image
    float *out,           //output image
    const int xdim,       //image width
    const int ydim,       //image height
    const double sigma,   //Gaussian sigma
    const int bc,       //boundary condition
    const int precision //defines the size of the window
)
{
    for(int i = 0; i < xdim*ydim; i++) out[i] = in[i];

    gaussian1D(out, xdim, ydim, sigma, bc, precision);
}

