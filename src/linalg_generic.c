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

complex_PRECISION global_inner_product_PRECISION( vector_PRECISION phi, vector_PRECISION psi, int start, int end, level_struct *l ) {
  
  complex_PRECISION local_alpha = 0, global_alpha = 0;
  local_alpha = process_inner_product_PRECISION( phi, psi, start, end, l );
  PROF_PRECISION_START( _GIP );

  if ( g.num_processes > 1 ) {
    PROF_PRECISION_START( _ALLR );
    MPI_Allreduce( &local_alpha, &global_alpha, 1, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
    PROF_PRECISION_STOP( _ALLR, 1 );
    PROF_PRECISION_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size );
    return global_alpha;
  } else {
    // all threads need the result of the norm
    PROF_PRECISION_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size );
    return local_alpha;
  }
}


complex_PRECISION process_inner_product_PRECISION( vector_PRECISION phi, vector_PRECISION psi, int start, int end, level_struct *l ) {
  
  PROF_PRECISION_START( _PIP );
  complex_PRECISION local_alpha = 0;
  
  
  for (int i = start; i < end; i ++)
    local_alpha += conj_PRECISION(phi[i])*psi[i];

  PROF_PRECISION_STOP( _PIP, (double)(end-start)/(double)l->inner_vector_size );

  return local_alpha;
}


void process_multi_inner_product_PRECISION( int count, complex_PRECISION *results, vector_PRECISION *phi, vector_PRECISION psi,
    int start, int end, level_struct *l ) {

  PROF_PRECISION_START( _PIP );
  int i;
  for(int c=0; c<count; c++)
    results[c] = 0.0;

  if ( l->depth == 0 ) {
    for(int c=0; c<count; c++)
      for ( i=start; i<end; )
        FOR12( results[c] += conj_PRECISION(phi[c][i])*psi[i]; i++; )
  } else {
    for(int c=0; c<count; c++)
      for ( i=start; i<end; )
        FOR2( results[c] += conj_PRECISION(phi[c][i])*psi[i]; i++; )
  }

  PROF_PRECISION_STOP( _PIP, (double)(end-start)/(double)l->inner_vector_size );
}


complex_PRECISION local_xy_over_xx_PRECISION( vector_PRECISION phi, vector_PRECISION psi, int start, int end, level_struct *l  ) {
  
  complex_PRECISION numerator = 0.0; PRECISION denominator = 0.0;
  
  VECTOR_FOR( int i=start, i<end, numerator += conj_PRECISION(phi[i])*psi[i]; denominator += NORM_SQUARE_PRECISION(phi[i]), i++, l );
  
  if ( abs_PRECISION(denominator) < EPS_PRECISION ) {
    return 0.0;
  }
  
  return numerator/denominator;
}

PRECISION global_norm_PRECISION( vector_PRECISION x, int start, int end, level_struct *l ) {
  return (PRECISION) sqrt((double) global_inner_product_PRECISION( x, x, start, end, l ));
}

PRECISION process_norm_PRECISION( vector_PRECISION x, int start, int end, level_struct *l ) {
  return (PRECISION) sqrt((double) process_inner_product_PRECISION( x, x, start, end, l ));
}


void vector_PRECISION_plus( vector_PRECISION z, vector_PRECISION x, vector_PRECISION y, int start, int end, level_struct *l ) {
  
  PROF_PRECISION_START( _LA2 );
  
  VECTOR_FOR( int i=start, i<end, z[i] = x[i] + y[i], i++, l );
  
  PROF_PRECISION_STOP( _LA2, (double)(end-start)/(double)l->inner_vector_size );
}


void vector_PRECISION_minus( vector_PRECISION z, vector_PRECISION x, vector_PRECISION y, int start, int end, level_struct *l ) {
  
  PROF_PRECISION_START( _LA2 );

  VECTOR_FOR( int i=start, i<end, z[i] = x[i] - y[i], i++, l );
  
  PROF_PRECISION_STOP( _LA2, (double)(end-start)/(double)l->inner_vector_size );
}

void vector_PRECISION_scale( vector_PRECISION z, vector_PRECISION x, complex_PRECISION alpha, int start, int end, level_struct *l ) {
  
  PROF_PRECISION_START( _LA6 );
  
  VECTOR_FOR( int i=start, i<end, z[i] = alpha*x[i], i++, l );
  
  PROF_PRECISION_STOP( _LA6, (double)(end-start)/(double)l->inner_vector_size );
}


void vector_PRECISION_real_scale( vector_PRECISION z, vector_PRECISION x, complex_PRECISION alpha,
                                  int start, int end, level_struct *l ) {
  
  PRECISION *r_z = (PRECISION*)z, *r_x = (PRECISION*)x, r_alpha = creal_PRECISION(alpha);
  int r_start = 2*start, r_end = 2*end;
  
  PROF_PRECISION_START( _LA2 );
  
  REAL_VECTOR_FOR( int i=r_start, i<r_end, r_z[i] = r_alpha*r_x[i], i++, l );
  
  PROF_PRECISION_STOP( _LA2, (double)(end-start)/(double)l->inner_vector_size );
}


void vector_PRECISION_copy( vector_PRECISION z, vector_PRECISION x, int start, int end, level_struct *l ) {
  
  PROF_PRECISION_START( _CPY );
  
  VECTOR_FOR( int i=start, i<end, z[i] = x[i], i++, l );
  
  PROF_PRECISION_STOP( _CPY, (double)(end-start)/(double)l->inner_vector_size );
}

void vector_PRECISION_saxpy( vector_PRECISION z, vector_PRECISION x, vector_PRECISION y, complex_PRECISION alpha, int start, int end, level_struct *l ) {
  
  PROF_PRECISION_START( _LA8 );
  
  VECTOR_FOR( int i=start, i<end, z[i] = x[i] + alpha*y[i], i++, l );
  
  PROF_PRECISION_STOP( _LA8, (double)(end-start)/(double)l->inner_vector_size );
}

void vector_PRECISION_multi_saxpy( vector_PRECISION z, vector_PRECISION *V, complex_PRECISION *alpha,
                               int sign, int count, int start, int end, level_struct *l ) {
  
  PROF_PRECISION_START( _LA8 );
  
  complex_PRECISION alpha_signed[count];
  for ( int c=0; c<count; c++ ) {
    alpha_signed[c] = sign*alpha[c];
  }
  
  for ( int c=0; c<count; c++ ) {
    for ( int i=start; i<end; ) {
      FOR12( z[i] += V[c][i]*alpha_signed[c]; i++; )
    }
  }
  
  PROF_PRECISION_STOP( _LA8, (PRECISION)(count) );
}

void vector_PRECISION_projection( vector_PRECISION z, vector_PRECISION v, int k, vector_PRECISION *W, complex_PRECISION *diag, 
                                  int orthogonal, level_struct *l ) {
  int j;

  vector_PRECISION v_tmp = NULL, *W_tmp = NULL;
  complex_PRECISION ip[k], ip_buffer[2*k];      
  
  MALLOC( v_tmp, complex_PRECISION, l->inner_vector_size );
  vector_PRECISION_define(v_tmp, 0, 0, l->inner_vector_size, l );
  
  MALLOC( W_tmp, complex_PRECISION*, k );
  W_tmp[0] = NULL; 
  MALLOC( W_tmp[0], complex_PRECISION, k*l->inner_vector_size );
  for ( j = 1; j<k; j++ )
    W_tmp[j] = W_tmp[0]+j*l->inner_vector_size;
  
  for ( j=0; j<k; j++ ) {
   vector_PRECISION_scale( W_tmp[j], W[j], diag[j], 0, l->inner_vector_size, l );
  }
  process_multi_inner_product_PRECISION( k, ip, W_tmp, v, 0, l->inner_vector_size, l );
  
  for ( j=0; j<k; j++ ) {
    ip_buffer[j] = ip[j];
  }
  MPI_Allreduce( ip_buffer, ip_buffer+k, k, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
  
  vector_PRECISION_multi_saxpy( v_tmp, W_tmp, ip_buffer+k, 1, k, 0, l->inner_vector_size, l );
   
  if (orthogonal) 
    vector_PRECISION_minus( z, v, v_tmp, 0, l->inner_vector_size, l );
  else
    vector_PRECISION_copy( z, v_tmp, 0, l->inner_vector_size, l );
  
  FREE( v_tmp, complex_PRECISION, l->inner_vector_size );
  FREE( W_tmp[0], complex_PRECISION, k*l->inner_vector_size );
  FREE( W_tmp, complex_PRECISION*, k );
}

void gram_schmidt_on_aggregates_PRECISION( vector_PRECISION *V, const int num_vect, level_struct *l ) {
  
  PROF_PRECISION_START( _GRAM_SCHMIDT_ON_AGGREGATES );
  int i, j, k, k1, k2, num_aggregates = l->s_PRECISION.num_aggregates,
      aggregate_size = l->inner_vector_size / num_aggregates, offset = l->num_lattice_site_var/2;
      
  complex_PRECISION alpha1, alpha2;
  vector_PRECISION v_pt1, v_pt2;
  PRECISION norm1, norm2;
      
  for ( j=0; j<num_aggregates; j++ ) {
    for ( k1=0; k1<num_vect; k1++ ) {
      v_pt1 = V[k1] + j*aggregate_size;
      
      for ( k2=0; k2<k1; k2++ ) {
        v_pt2 = V[k2] + j*aggregate_size;
        alpha1 = 0; alpha2 = 0;
        // V[k1] -= <V[k2],V[k1]> V[k2] | 2*j-th and 2*j+1-st aggregate
        for ( i=0; i<aggregate_size; ) {
          for ( k=0; k<offset; k++, i++ )
            alpha1 += conj_PRECISION(v_pt2[i]) * v_pt1[i];
          for ( k=0; k<offset; k++, i++ )
            alpha2 += conj_PRECISION(v_pt2[i]) * v_pt1[i];
        }
        for ( i=0; i<aggregate_size; ) {
          for ( k=0; k<offset; k++, i++ )
            v_pt1[i] -=  alpha1 * v_pt2[i];
          for ( k=0; k<offset; k++, i++ )
            v_pt1[i] -=  alpha2 * v_pt2[i];
        }
      }
      
      norm1 = 0; norm2 = 0;
      // V[k1] = V[k1]/norm(V[k1]) | 2*j-th and 2*j+1-st aggregate    
      for ( i=0; i<aggregate_size; ) {
        for ( k=0; k<offset; k++, i++ )
          norm1 += NORM_SQUARE_PRECISION(v_pt1[i]);
        for ( k=0; k<offset; k++, i++ )
          norm2 += NORM_SQUARE_PRECISION(v_pt1[i]);
      }
      norm1 = 1/sqrt(norm1); norm2 = 1/sqrt(norm2);
      for ( i=0; i<aggregate_size; ) {
        for ( k=0; k<offset; k++, i++ )
          v_pt1[i] =  norm1 * creal_PRECISION(v_pt1[i]) + I*norm1* cimag_PRECISION(v_pt1[i]);
        for ( k=0; k<offset; k++, i++ )
          v_pt1[i] =  norm2 * creal_PRECISION(v_pt1[i]) + I*norm2* cimag_PRECISION(v_pt1[i]);
      }
    }
  }
  PROF_PRECISION_STOP( _GRAM_SCHMIDT_ON_AGGREGATES, 1 );
}


void spinwise_PRECISION_skalarmultiply( vector_PRECISION eta1, vector_PRECISION eta2, vector_PRECISION phi, complex_PRECISION alpha,
                                        int start, int end, level_struct *l ) {
  
  PROF_PRECISION_START( _LA6 );  
  for ( int i=start; i<end; ) {
    FOR6( eta1[i] = alpha*phi[i]; eta2[i] = _COMPLEX_PRECISION_ZERO; i++; )
    FOR6( eta2[i] = alpha*phi[i]; eta1[i] = _COMPLEX_PRECISION_ZERO; i++; )
  }
  PROF_PRECISION_STOP( _LA6, 1 );
}


void set_boundary_PRECISION( vector_PRECISION phi, complex_PRECISION alpha, level_struct *l ) {
  
  PROF_PRECISION_START( _SET );
  
  for (int i = l->inner_vector_size; i < l->vector_size; i ++)
    phi[i] = alpha;
  
  PROF_PRECISION_STOP( _SET, (double)(l->vector_size-l->inner_vector_size)/(double)l->inner_vector_size );
}


void gram_schmidt_PRECISION( vector_PRECISION *V, complex_PRECISION *buffer, const int begin, const int n, level_struct *l ) {
  
  // NOTE: only thread safe, if "buffer" is the same buffer for all threads belonging to a common MPI process
  PROF_PRECISION_START( _LA );
  
  PRECISION beta;
  int i, j, start = 0, end = l->inner_vector_size;
  
  for ( i=begin; i<n; i++ ) {
    
    complex_PRECISION tmp[i];
    process_multi_inner_product_PRECISION( i, tmp, V, V[i], 0, l->inner_vector_size, l );
    for ( j=0; j<i; j++ ) {
      buffer[j] = tmp[j];
    }
    
    if ( i>0 ) {
      PROF_PRECISION_START( _ALLR );
      MPI_Allreduce( buffer, buffer+n, i, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
      PROF_PRECISION_STOP( _ALLR, 1 );
    }
    
    for( j=0; j<i; j++ ) {
      vector_PRECISION_saxpy( V[i], V[i], V[j], -(buffer+n)[j], start, end, l );
    }
    
      
    beta = global_norm_PRECISION( V[i], 0, l->inner_vector_size, l );
    vector_PRECISION_real_scale( V[i], V[i], creal(1.0/beta), start, end, l );
  }
  
  PROF_PRECISION_STOP( _LA, 1 );
}


void setup_gram_schmidt_PRECISION_compute_dots(
    complex_PRECISION *thread_buffer, vector_PRECISION *V, int count, int offset,
    int start, int end, level_struct *l) {

  int cache_block_size = 12*64;
  complex_PRECISION tmp[cache_block_size];

  for(int i=0; i<2*offset; i++)
    thread_buffer[i] = 0.0;

  
  for ( int i=start; i<end; i+=cache_block_size) {
    coarse_gamma5_PRECISION( tmp, V[count]+i, 0, cache_block_size, l );
    for ( int j=0; j<count; j++ ) {
      for ( int k=0; k<cache_block_size; k++) {
        thread_buffer[j]   += conj_PRECISION(V[j][i+k])*V[count][i+k];
        thread_buffer[j+offset] += conj_PRECISION(V[j][i+k])*tmp[k];
      }
    }
  }

}


void setup_gram_schmidt_PRECISION_axpys(
    complex_PRECISION *thread_buffer, vector_PRECISION *V, int count, int offset,
    int start, int end, level_struct *l) {
  
  int cache_block_size = 12*64;
  complex_PRECISION tmp[cache_block_size];

  for ( int i=start; i<end; i+=cache_block_size) {
    for ( int j=0; j<count; j++ ) {
      coarse_gamma5_PRECISION( tmp, V[j]+i, 0, cache_block_size, l );
      for ( int k=0; k<cache_block_size; k++) {
        V[count][i+k] -= thread_buffer[2*offset+j]*V[j][i+k];
        V[count][i+k] -= thread_buffer[3*offset+j]*tmp[k];
      }
    }
  }
}
