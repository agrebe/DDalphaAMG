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

void rhs_define( vector_double rhs, level_struct *l ) {

  int start = 0;
  int end = l->num_inner_lattice_sites * l->num_lattice_site_var;

  if ( g.rhs == 0 ) {
    vector_double_define( rhs, 1, start, end, l );
    if ( g.print > 0 ) printf0("rhs = ones\n");
  } else if ( g.rhs == 1 )  {
    vector_double_define( rhs, 0, start, end, l );
    if ( g.my_rank == 0 ) {
      rhs[0] = 1.0;
    }
    if ( g.print > 0 ) printf0("rhs = first unit vector\n");
  } else if ( g.rhs == 2 ) {
    // this would yield different results if we threaded it, so we don't
    vector_double_define_random( rhs, 0, l->inner_vector_size, l );
    if ( g.print > 0 ) printf0("rhs = random\n");
  } else if ( g.rhs == 3 ) {
    vector_double_define( rhs, 0, start, end, l );
  } else {
    ASSERT( g.rhs >= 0 && g.rhs <= 4 );
  }
}


int wilson_driver( vector_double solution, vector_double source, level_struct *l ) {
  
  int iter = 0, start = 0, end = l->num_inner_lattice_sites * l->num_lattice_site_var;
  
  vector_double rhs = g.mixed_precision==2?g.p_MP.dp.b:g.p.b;
  vector_double sol = g.mixed_precision==2?g.p_MP.dp.x:g.p.x;

#ifdef WILSON_BENCHMARK
  prof_init( l );
  double t = -MPI_Wtime();
  double t_min = 1000;
  for ( int i=0; i<100; i++ ) {
    double tmp_t = -MPI_Wtime();
#endif
  
  vector_double_copy( rhs, source, start, end, l );  
  if ( g.method == -1 ) {
    cgn_double( &(g.p), l );
  } else if ( g.mixed_precision == 2 ) {
    iter = fgmres_MP( &(g.p_MP), l );
  } else {
    iter = fgmres_double( &(g.p), l );
  }
  vector_double_copy( solution, sol, start, end, l );
#ifdef WILSON_BENCHMARK
    tmp_t += MPI_Wtime();
    if ( tmp_t < t_min )
      t_min = tmp_t;
  }
  t +=MPI_Wtime();
  printf0("average over 100 solves: %lf seconds\n", t/100 );
  printf0("minimum out of 100 solves: %lf seconds\n", t_min );
  prof_print( l );
#endif
  
  return iter;
}


void solve( vector_double solution, vector_double source, level_struct *l ) {
  
  vector_double rhs = g.mixed_precision==2?g.p_MP.dp.b:g.p.b;
  wilson_driver( solution, source, l );
}


void solve_driver( level_struct *l ) {
  
  vector_double solution = NULL, source = NULL;
  
  MALLOC( solution, complex_double, l->inner_vector_size );
  MALLOC( source, complex_double, l->inner_vector_size );
  
  int Lx = g.local_lattice[0][X];
  int Ly = g.local_lattice[0][Y];
  int Lz = g.local_lattice[0][Z];
  int Lt = g.local_lattice[0][T];
  double * pion_correlator = malloc(sizeof(double) * Lt);
  for (int t = 0; t < Lt; t ++) pion_correlator[t] = 0;

  for (int j = 0; j < 12; j ++) {
    rhs_define( source, l );
    if ( g.my_rank == 0 ) {
      source[0] = 0;
      source[j] = 1.0;
    }

    solve( solution, source, l );

    // compute pion correlator
    // do this on one thread
//    printf("First 24 entries in solution: \n");
//    for (int i = 0; i < 24; i ++)
//      printf("%.6e %.6e\n", creal(solution[i]), cimag(solution[i]));
    printf("First entry in solution: %.5e\n", creal(solution[0]));
    double * correlator = malloc(sizeof(double) * Lt);
    for (int t = 0; t < Lt; t ++) correlator[t] = 0;
    for (int x = 0; x < Lx; x ++)
      for (int y = 0; y < Ly; y ++)
        for (int z = 0; z < Lz; z ++)
          for (int t = 0; t < Lt; t ++) {
            //int location = ((x * Lx + y) * Ly + z) * Lz + t;
            int location = ((t * Lz + z) * Ly + y) * Lx + x;
            for (int i = 0; i < 24; i ++) {
              double value = ((double *) solution)[location * 24 + i];
//              if (value != 0)
//                printf("Nonzero value at (%d, %d, %d, %d, %d): %e\n", x, y, z, t, i, value);
              correlator[t] += value * value;
              pion_correlator[t] += value * value;
            }
          }
    for (int t = 0; t < Lt; t ++)
      printf("C[%d] = %e\n", t, correlator[t]);
    free(correlator);
  }
  for (int t = 0; t < Lt; t ++)
    printf("C[%d] = %e\n", t, pion_correlator[t]);
  free(pion_correlator);


  FREE( solution, complex_double, l->inner_vector_size );
  FREE( source, complex_double, l->inner_vector_size );
}

