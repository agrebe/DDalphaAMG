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

#ifndef IO_HEADER
  #define IO_HEADER
  
  void byteswap( char *in );  
  void byteswap8( char *in );  
  void read_conf( double *input_data, char *input_name, double *conf_plaq, level_struct *l );
  void read_conf_multi( double *input_data, char *input_name, double *conf_plaq, level_struct *l );
  void vector_io( double *phi, char *filename, const int mode, level_struct *l );
  void vector_io_single_file( double *psi, double *lambda, char *filename, const int mode, int n, char *vector_type, level_struct *l );
  void d_dump( config_double D, level_struct *l );
  
  static inline int process_index( int t, int z, int y, int x, int ll[4] ) {
    
    if ( g.num_processes > 1 ) {
      int pcoords[4]; int rank;
      pcoords[T] = t/ll[T]; pcoords[Z] = z/ll[Z];
      pcoords[Y] = y/ll[Y]; pcoords[X] = x/ll[X];
      MPI_Cart_rank( g.comm_cart, pcoords, &rank );
      return rank;
    } else {
      return 0;
    }
  }
  
#endif
