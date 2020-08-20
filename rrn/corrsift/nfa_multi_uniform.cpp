/*
 * Copyright (c) 2016, Axel Davy
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *     names of its contributors may be used to endorse or promote products
 *     derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED ''AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include <stdlib.h>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>

#include "nfa_multi_uniform.h"

using boost::math::chi_squared_distribution;
using boost::math::normal_distribution;

void compute_nfa_multi_uniform(float **dst, float **uniforms, int nuniforms, 
			       float * map, int npix)
{
    int i, j;
    normal_distribution<float> gaussian; /* mean 0, std 1 */
    chi_squared_distribution<float> chi2_distr(nuniforms);

    if (*dst == NULL)
        *dst = (float *)malloc(npix * sizeof(float));

    for (i = 0; i < npix; i++) 
      {
        float chi2 = 0;
        try {
           for (j = 0; j < nuniforms; j++)
	     {
	       if (map[i] > 0)
	         {
		    float v = quantile(gaussian, uniforms[j][i]);
		    chi2 += v*v;
	         }
	     }
           (*dst)[i] = cdf(complement(chi2_distr, chi2)) * npix;
           //(*dst)[i] = cdf(complement(chi2_distr, chi2)) / npix;
        } catch (boost::exception const&  ex) // overflow, etc
        {
            (*dst)[i] = 0;
        }
      }
}
