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

#ifndef _NFA_MULTI_UNIFORM
#define _NFA_MULTI_UNIFORM

/* compute_nfa_multi_uniform
 * Computes a table filled with the NFA computation.
 * Then entries are supposed to be nuniforms tables of uniforms.
 * The tables are of length N.
 * Each NFA is computed on a tuple of nuniforms uniforms.
 * The uniforms are converted to gaussians, are squared and then
 * summed. The NFA is then computed on the resulting chi2 distribution.
 * dst: pointer to the output table. If pointer is NULL, the table is allocated.
 * uniforms: table of pointers to table of uniforms.
 * nuniforms: size of uniforms
 * N: size of the table of uniforms. size of the output.
 */


void compute_nfa_multi_uniform(float **dst, float **uniforms, int nuniforms, 
			       float * map, int N); 


#endif /* _NFA_MULTI_UNIFORM */
