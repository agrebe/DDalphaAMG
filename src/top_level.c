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

void rhs_define( vector_double rhs, level_struct *l, struct Thread *threading ) {
  
  // no hyperthreading here
  if(threading->thread != 0)
    return;

  int start = threading->start_index[l->depth];
  int end = threading->end_index[l->depth];

  if ( g.rhs == 0 ) {
    vector_double_define( rhs, 1, start, end, l );
    START_MASTER(threading)
    if ( g.print > 0 ) printf0("rhs = ones\n");
    END_MASTER(threading)
  } else if ( g.rhs == 1 )  {
    vector_double_define( rhs, 0, start, end, l );
    if ( g.my_rank == 0 ) {
      START_LOCKED_MASTER(threading)
      rhs[0] = 1.0;
      END_LOCKED_MASTER(threading)
    }
    START_MASTER(threading)
    if ( g.print > 0 ) printf0("rhs = first unit vector\n");
    END_MASTER(threading)
  } else if ( g.rhs == 2 ) {
    // this would yield different results if we threaded it, so we don't
    START_LOCKED_MASTER(threading)
    vector_double_define_random( rhs, 0, l->inner_vector_size, l );
    END_LOCKED_MASTER(threading)
    START_MASTER(threading)
    if ( g.print > 0 ) printf0("rhs = random\n");
    END_MASTER(threading)
  } else if ( g.rhs == 3 ) {
    vector_double_define( rhs, 0, start, end, l );
  } else {
    ASSERT( g.rhs >= 0 && g.rhs <= 4 );
  }
}


int wilson_driver( vector_double solution, vector_double source, level_struct *l, struct Thread *threading ) {
  
  int iter = 0, start = threading->start_index[l->depth], end = threading->end_index[l->depth];
  
  vector_double rhs = g.mixed_precision==2?g.p_MP.dp.b:g.p.b;
  vector_double sol = g.mixed_precision==2?g.p_MP.dp.x:g.p.x;

#ifdef WILSON_BENCHMARK
  START_MASTER(threading)
  prof_init( l );
  END_MASTER(threading)
  double t = -MPI_Wtime();
  double t_min = 1000;
  for ( int i=0; i<100; i++ ) {
    double tmp_t = -MPI_Wtime();
#endif
  
  vector_double_copy( rhs, source, start, end, l );  
  if ( g.method == -1 ) {
    cgn_double( &(g.p), l, threading );
  } else if ( g.mixed_precision == 2 ) {
    iter = fgmres_MP( &(g.p_MP), l, threading );
  } else {
    iter = fgmres_double( &(g.p), l, threading );
  }
  vector_double_copy( solution, sol, start, end, l );
#ifdef WILSON_BENCHMARK
    tmp_t += MPI_Wtime();
    if ( tmp_t < t_min )
      t_min = tmp_t;
  }
  t +=MPI_Wtime();
  START_MASTER(threading)
  printf0("average over 100 solves: %lf seconds\n", t/100 );
  printf0("minimum out of 100 solves: %lf seconds\n", t_min );
  prof_print( l );
  END_MASTER(threading)
#endif
  
  return iter;
}


void solve( vector_double solution, vector_double source, level_struct *l, struct Thread *threading ) {
  
  vector_double rhs = g.mixed_precision==2?g.p_MP.dp.b:g.p.b;

  if ( g.vt.evaluation ) {
    // this would yield different results if we threaded it, so we don't
    START_LOCKED_MASTER(threading)
    vector_double_define_random( rhs, 0, l->inner_vector_size, l );
    scan_var( &(g.vt), l );
    END_LOCKED_MASTER(threading)
  } else {
    wilson_driver( solution, source, l, threading );
  }
}


void solve_driver( level_struct *l, struct Thread *threading ) {
  
  vector_double source   = g.mixed_precision==2?g.p_MP.dp.b:g.p.b;
  vector_double solution = g.mixed_precision==2?g.p_MP.dp.x:g.p.x;
  
  int Lx = g.local_lattice[0][X];
  int Ly = g.local_lattice[0][Y];
  int Lz = g.local_lattice[0][Z];
  int Lt = g.local_lattice[0][T];
  double * pion_correlator = malloc(sizeof(double) * Lt);
  START_LOCKED_MASTER(threading)
  for (int t = 0; t < Lt; t ++) pion_correlator[t] = 0;
  END_LOCKED_MASTER(threading)

  for (int j = 0; j < 12; j ++) {
    rhs_define( source, l, threading );
    if ( g.my_rank == 0 ) {
      START_LOCKED_MASTER(threading)
      source[0] = 0;
      source[j] = 1.0;
      END_LOCKED_MASTER(threading)
    }

    if ( g.method == -1 ) {
      cgn_double( &(g.p), l, threading );
    } else if ( g.mixed_precision == 2 ) {
      fgmres_MP( &(g.p_MP), l, threading );
    } else {
      fgmres_double( &(g.p), l, threading );
    }

    // compute pion correlator
    // do this on one thread
    START_LOCKED_MASTER(threading)
    printf("First entry in solution: %.14e\n", creal(solution[0]));
    double * correlator = malloc(sizeof(double) * Lt);
    for (int t = 0; t < Lt; t ++) correlator[t] = 0;
    for (int x = 0; x < Lx; x ++)
      for (int y = 0; y < Ly; y ++)
        for (int z = 0; z < Lz; z ++)
          for (int t = 0; t < Lt; t ++) {
            int location = ((t * Lz + z) * Ly + y) * Lx + x;
            for (int i = 0; i < 24; i ++) {
              double value = ((double *) solution)[location * 24 + i];
              correlator[t] += value * value;
              pion_correlator[t] += value * value;
            }
          }
    for (int t = 0; t < Lt; t ++)
      printf("C[%d] = %e\n", t, correlator[t]);
    free(correlator);
    END_LOCKED_MASTER(threading)
  }
  START_LOCKED_MASTER(threading)
  for (int t = 0; t < Lt; t ++)
    printf("C[%d] = %e\n", t, pion_correlator[t]);
  END_LOCKED_MASTER(threading)
  free(pion_correlator);
}

