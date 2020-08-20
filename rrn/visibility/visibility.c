/*----------------------------------------------------------------------------

  Copyright (c) 2019 Charles Hessel <hesselcharles@gmail.com>
  Copyright (c) 2018-2019 rafael grompone von gioi <grompone@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.

  ----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------

  CHANGELOG:

  - december 2019: Charles Hessel
    Added the ASLIB horror
  - january 2020: Charles Hessel
    Added function "angles" which reuses Rafael's code

  ----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iio.h"

/*----------------------------------------------------------------------------*/
#ifndef FALSE
#define FALSE 0
#endif /* !FALSE */

#ifndef TRUE
#define TRUE 1
#endif /* !TRUE */

/*----------------------------------------------------------------------------*/
/* PI */
#ifndef M_PI
#define M_PI   3.14159265358979323846
#endif /* !M_PI */

/*----------------------------------------------------------------------------*/
/* Label for pixels with undefined gradient. */
#define NOTDEF -1024.0

/*----------------------------------------------------------------------------*/
/* fatal error, print a message to standard error and exit
 */
static void error(char * msg)
{
  fprintf(stderr,"error: %s\n",msg);
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*/
/* memory allocation, print an error and exit if fail
 */
static void * xmalloc(size_t size)
{
  void * p;
  if( size == 0 ) error("xmalloc input: zero size");
  p = malloc(size);
  if( p == NULL ) error("out of memory");
  return p;
}

/*----------------------------------------------------------------------------*/
/* Euclidean distance between x1,y1 and x2,y2
 */
static double dist(double x1, double y1, double x2, double y2)
{
  return sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );
}

/*----------------------------------------------------------------------------*/
/* Normalized angle difference between 'a' and the symmetric of 'b'
   relative to a vertical axis.
 */
static double norm_angle(double a, double b)
{
  a -= b;
  while( a <= -M_PI ) a += 2.0*M_PI;
  while( a >   M_PI ) a -= 2.0*M_PI;

  return fabs(a) / M_PI;
}

/*----------------------------------------------------------------------------*/
/*                                    main                                    */
/*             (or python interface, depending on ASLIB variable)             */
/*----------------------------------------------------------------------------*/
#ifdef ASLIB
int visibility(int X, /* width of image */
               int Y, /* height of image */
               int N, /* number of images */
               double *images, /* all images (input) */
               double *mask) /* output masks, same shape as 'images' */
#else
int main(int argc, char ** argv)
#endif
{
#ifdef ASLIB
#else
  double * mask;
  int X,Y,N;
  int XX,YY,CC;
#endif
  int S = 500;
  double rho = 0.25;
  double * image;
  double * angles;
  int * done;
  double * adiff;
  int * reg_x;
  int * reg_y;
  int reg_n;
  int n,x,y,k,i,j;
  int W = 14;
  double logNT,logNFA;
  int neighbor_x[8] = {  0,  0, -1,  1, -1, -1,  1,  1 };
  int neighbor_y[8] = { -1,  1,  0,  0, -1,  1, -1,  1 };

#ifdef ASLIB
#else
  /* usage */
  if( argc < 3 ) error("usage: visibility image1 image2 [image3 ... imageN]");

  /* read input */
  N = argc - 1;
#endif
  for(n=0; n<N; n++)
    {
#ifdef ASLIB
      image = images + n*X*Y;
#else
      image = iio_read_image_double_split(argv[n+1], &XX, &YY, &CC);

      if( XX < 1 || YY < 1 ) error("invalid image size");
      if( CC > 1 ) error("only single channel images handled");
      if( n>0 && ( XX!=X || YY!=Y ) ) error("images must have the same size");
#endif

      /* get memory */
      if( n==0 )
        {
#ifdef ASLIB
#else
          X = XX;
          Y = YY;
          mask = (double *) xmalloc( X * Y * N * sizeof(double) );
#endif
          angles = (double *) xmalloc( X * Y * N * sizeof(double) );
          adiff = (double *) xmalloc( X * Y * sizeof(double) );
          reg_x = (int *) xmalloc( X * Y * sizeof(int) );
          reg_y = (int *) xmalloc( X * Y * sizeof(int) );
          done = (int *) xmalloc( X * Y * sizeof(int) );
          for(k=0; k<X*Y*N; k++) mask[k] = 255.0; /* initialize as cloud */
        }

      /* compute gradient angles */
      for(x=1; x<X-1; x++)
      for(y=1; y<Y-1; y++)
        {
          double dx = image[(x+1)+y*X] - image[(x-1)+y*X];
          double dy = image[x+(y+1)*X] - image[x+(y-1)*X];
          double mod = sqrt(dx*dx + dy*dy);

          if( mod <= 0.0 ) angles[x+y*X + n*X*Y] = NOTDEF;
          else             angles[x+y*X + n*X*Y] = atan2(dy,dx);
        }

#ifdef ASLIB
#else
      free( (void *) image );
#endif
    }

  /* find matches */
  for(n=0; n<N; n++)
    for(k=n+1; k<N; k++)
      {
        double s;

        /* compute gradient angle difference */
        for(i=0; i<X*Y; i++) adiff[i] = NOTDEF;
        for(x=1; x<X-1; x++)
        for(y=1; y<Y-1; y++)
          {
            double a = angles[x+y*X + n*X*Y];
            double b = angles[x+y*X + k*X*Y];
            if( a != NOTDEF && b != NOTDEF ) adiff[x+y*X] = norm_angle(a,b);
            else                             adiff[x+y*X] = 1.0;
          }

        /* grow regions as potential match candidates */
        for(x=0; x<X; x++)
        for(y=0; y<Y; y++)
          if( adiff[x+y*X] != NOTDEF && adiff[x+y*X] <= rho )
            {
              /* initialize region */
              s = adiff[x+y*X];
              adiff[x+y*X] = NOTDEF;
              reg_x[0] = x;
              reg_y[0] = y;
              reg_n = 1;

              /* recursively add neighbors of the that match the criterion */
              for(i=0; i<reg_n; i++)
                for(j=0; j<4; j++)
                  {
                    int xx = reg_x[i] + neighbor_x[j];
                    int yy = reg_y[i] + neighbor_y[j];

                    if( xx >= 0 && xx < X && yy >= 0 && yy < Y &&
                        adiff[xx+yy*X] != NOTDEF && adiff[xx+yy*X] <= rho )
                      {
                        s += adiff[xx+yy*X];
                        adiff[xx+yy*X] = NOTDEF;
                        reg_x[reg_n] = xx;
                        reg_y[reg_n] = yy;
                        ++reg_n;
                      }
                  }

              /* Number of test:
                 all pair of images ~ N^2                    ->  N^2
                 each pixel is the center of possible region ->  X*Y
                 we need to test regions of size 1 to X*Y    ->  X*Y
                 the number of polyominos of size n is about
                     0.316915 * 4.0625696^n / n              -> 0.3*4.06^n / n
               */
              logNT = 2.0 * log10(N) + 2.0 * log10(X) + 2.0 * log10(Y)
                    + log10(0.316915) + (double) reg_n * log10(4.0625696)
                    - log10( (double) reg_n );

              /* NFA = NT * s^n / n!
                   log(n!) is bounded by Stirling's approximation:
                   n! >= sqrt(2pi) * n^(n+0.5) * exp(-n)
                   then, log10(NFA) <= log10(NT) + n*log10(s)
                                                 - log10(latter expansion) */
              logNFA = logNT + (double) reg_n * log10(s)
                     - 0.5 * log10(2.0 * M_PI)
                     - ( (double) reg_n + 0.5 ) * log10( (double) reg_n )
                     + (double) reg_n * log10(exp(1.0));

              /* if a meaningful match is found, no cloud could be present */
              if( reg_n > 1 && logNFA < 0.0 )
                for(i=0; i<reg_n; i++)
                  {
                    int xx = reg_x[i];
                    int yy = reg_y[i];
                    mask[xx+yy*X + n*X*Y] = mask[xx+yy*X + k*X*Y] = 0.0;
                  }
            }
      }

  /* grain filter: remove connected regions in mask with less than S pixels */
  for(n=0; n<N; n++)
    {
      for(i=0; i<X*Y; i++) done[i] = FALSE;

      for(x=1; x<X-1; x++)
      for(y=1; y<Y-1; y++)
        if( ! done[x+y*X] )
          {
            /* initialize region */
            double val = mask[x+y*X + n*X*Y];
            if( val == 0.0 ) continue; /* Only remove small clouds */
            done[x+y*X] = TRUE;
            reg_x[0] = x;
            reg_y[0] = y;
            reg_n = 1;

            /* recursively add neighbors of the that match the criterion */
            for(i=0; i<reg_n; i++)
            for(j=0; j<4; j++)
              {
                int xx = reg_x[i] + neighbor_x[j];
                int yy = reg_y[i] + neighbor_y[j];

                if( xx > 0 && xx < X-1 && yy > 0 && yy < Y-1 &&
                    ! done[xx+yy*X] && mask[xx+yy*X + n*X*Y] == val )
                  {
                    done[xx+yy*X] = TRUE;
                    reg_x[reg_n] = xx;
                    reg_y[reg_n] = yy;
                    ++reg_n;
                  }
              }

            /* remove connected region if smaller than S */
            if( reg_n < S )
              for(i=0; i<reg_n; i++)
                mask[reg_x[i] + reg_y[i]*X + n*X*Y]
                  = val == 255.0 ? 0.0 : 255.0;
          }
    }

#ifdef ASLIB
#else
  /* write cloud masks */
  for(n=0; n<N; n++)
    {
      char filename[512];
      sprintf(filename,"%03d.png",n);
      iio_write_image_double_vec(filename, mask + n*X*Y, X, Y, 1);
    }

  /* free memory */
  free( (void *) mask );
#endif
  free( (void *) angles );
  free( (void *) adiff );
  free( (void *) reg_x );
  free( (void *) reg_y );
  free( (void *) done );

  return EXIT_SUCCESS;
}
/*----------------------------------------------------------------------------*/
#ifdef ASLIB
int angles(int X, /* width of image */
           int Y, /* height of image */
           double *im0, /* first image */
           double *im1, /* second image */
           double *adiff) /* normalized angle difference */
{
  int N = 2; /* for two images: im0 and im1 */
  double * ims[2] = {im0, im1}; /* the two images */
  double * image; /* holds one of the two images */
  double * angles;
  int n,x,y,k,i,j;

  /* get memory */
  angles = (double *) xmalloc( X * Y * N * sizeof(double) );


  for(n=0; n<N; n++)
    {
      image = ims[n];

      /* compute gradient angles */
      for(x=1; x<X-1; x++)
      for(y=1; y<Y-1; y++)
        {
          double dx = image[(x+1)+y*X] - image[(x-1)+y*X];
          double dy = image[x+(y+1)*X] - image[x+(y-1)*X];
          double mod = sqrt(dx*dx + dy*dy);

          if( mod <= 0.0 ) angles[x+y*X + n*X*Y] = NOTDEF;
          else             angles[x+y*X + n*X*Y] = atan2(dy,dx);
        }

    }

  /* compute gradient angle difference */
  for(i=0; i<X*Y; i++) adiff[i] = NOTDEF;
  for(x=1; x<X-1; x++)
  for(y=1; y<Y-1; y++)
    {
      double a = angles[x+y*X];
      double b = angles[x+y*X + X*Y];
      if( a != NOTDEF && b != NOTDEF ) adiff[x+y*X] = norm_angle(a,b);
      else                             adiff[x+y*X] = 1.0;
    }

  free( (void *) angles );

  return EXIT_SUCCESS;
}
#endif
/*----------------------------------------------------------------------------*/
