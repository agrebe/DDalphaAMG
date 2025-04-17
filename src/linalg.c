/*
 * Copyright (C) 2016, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern Leder, Issaku Kanamori.
 * 
 * This file is part of the DDalphaAMG solver library.
 * 
 * The DDalphaAMG solver library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * The DDalphaAMG solver library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * 
 * You should have received a copy of the GNU General Public License
 * along with the DDalphaAMG solver library. If not, see http://www.gnu.org/licenses/.
 * 
 */

#include "main.h"

void process_multi_inner_product_MP( int count, complex_double *results, vector_float *phi,
                                     vector_float psi, int start, int end, level_struct *l,
                                     struct Thread *threading ) {

  PROF_float_START( _PIP, threading );
  int i;
  for(int c=0; c<count; c++)
    results[c] = 0.0;

  for(int c=0; c<count; c++) {
    for ( i=start; i<end; i++) {
      results[c] += (complex_double) conj_float(phi[c][i])*psi[i];
    }
  }

  PROF_float_STOP( _PIP, (double)(end-start)/(double)l->inner_vector_size, threading );
}

double global_norm_MP( vector_float x, int start, int end, level_struct *l, struct Thread *threading ) {
  
  PROF_float_START( _GIP, threading );
  
  int i;
  double local_alpha = 0, global_alpha = 0;

  SYNC_CORES(threading)
  for ( i=start; i<end; i++)
    local_alpha += (complex_double) NORM_SQUARE_float(x[i]);

  if ( g.num_processes > 1 ) {
    PROF_double_START( _ALLR );
    MPI_Allreduce( &local_alpha, &global_alpha, 1, MPI_double, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_float.level_comm );
    PROF_double_STOP( _ALLR, 1 );
    PROF_float_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );
    return sqrt((double)global_alpha);
  } else {
    PROF_float_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );
    return sqrt((double)local_alpha);
  }
}
