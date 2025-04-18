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

/*
 * 2016 Aug. 25: Issaku Kanamori
 *   modified read_conf(), to avoid reading from non-master nodes
 *
 * 2016 Oct. 3:  Issaku Kanamori
 *   added read_conf_multi()
 */

#include "main.h"

void byteswap( char *in ) {

  char tmp[4];
  tmp[0] = in[3];
  tmp[1] = in[2];
  tmp[2] = in[1];
  tmp[3] = in[0];

  in[0] = tmp[0];
  in[1] = tmp[1];
  in[2] = tmp[2];
  in[3] = tmp[3];
}
void byteswap8(char *in ) {
  
  char tmp[8];
  tmp[0] = in[7];
  tmp[1] = in[6];
  tmp[2] = in[5];
  tmp[3] = in[4];
  tmp[4] = in[3];
  tmp[5] = in[2];
  tmp[6] = in[1];
  tmp[7] = in[0];

  in[0] = tmp[0];
  in[1] = tmp[1];
  in[2] = tmp[2];
  in[3] = tmp[3];
  in[4] = tmp[4];
  in[5] = tmp[5];
  in[6] = tmp[6];
  in[7] = tmp[7];
}

void read_conf( double *input_data, char *input_name, double *conf_plaq, level_struct *l ) {
  
/*********************************************************************************
* Reads in the configuration.
* - double *input_data: Variable where conf data is stored.
* - char *input_name: Name of the input file.
* - double *conf_plaq: Holds the plaquette of given configuration.                                                
*********************************************************************************/

  int t, z, y, x, i, j, k, mu, lsize[4], desired_rank,
      *gl=l->global_lattice, *ll=l->local_lattice, read_size = 4*18*ll[X];
  double *input_data_pt, plaq;
  FILE* fin = NULL;
  MPI_Request sreq;
  confbuffer_struct buffer[2]; // Having two buffers allows communication hiding.
  confbuffer_struct *buffer_pt = NULL;
  
  buffer[0].next = &(buffer[1]);
  buffer[1].next = &(buffer[0]);  
  buffer[0].data = NULL;
  buffer[1].data = NULL;
  
  if ( g.my_rank == 0 ) {
    MALLOC( buffer[0].data, double, read_size );
    MALLOC( buffer[1].data, double, read_size );
    buffer_pt = &buffer[0];
  
    ASSERT( (fin = fopen( input_name, "rb" )) != NULL );
  
    // Read in lattice size in each dimension
    ASSERT( fread( lsize, sizeof(int), 4, fin ) > 0 );
#ifdef BIG_ENDIAN_CNFG
    byteswap( (char *) &lsize[0] );  
    byteswap( (char *) &lsize[1] );
    byteswap( (char *) &lsize[2] );
    byteswap( (char *) &lsize[3] );
#endif
    for ( mu=0; mu<4; mu++ )
      ASSERT( lsize[mu] == gl[mu] );
  
  // Read in plaquette
    ASSERT( fread( &plaq, sizeof(double), 1, fin ) > 0 );
#ifdef BIG_ENDIAN_CNFG
    byteswap8( (char *) &plaq );
#endif
    printf0("\nDesired average plaquette: %.13lf in [0,3]\n", plaq );
    printf0("\nDesired average plaquette: %.13lf in [0,1]\n", plaq/3.0 );
    *conf_plaq = plaq; 
  }
  MPI_Bcast(conf_plaq, 1, MPI_DOUBLE, 0, g.comm_cart);

  input_data_pt = input_data;
  
  // Distribute data to according processes
  if ( g.my_rank == 0 ) 
    ASSERT( fread( buffer_pt->data, sizeof(double), read_size, fin ) > 0 );
  
  k = 0;
  for ( t=0; t<gl[T]; t++ )
    for ( z=0; z<gl[Z]; z++ )
      for ( y=0; y<gl[Y]; y++ )
        for ( x=0; x<gl[X]; x+=ll[X] ) {
          desired_rank = process_index( t, z, y, x, ll );
          
          if ( g.my_rank == 0 ) {
            MPI_Isend( buffer_pt->data, read_size, MPI_DOUBLE, desired_rank, k, g.comm_cart, &sreq );
            if ( ! ( t == gl[T]-1 && z == gl[Z]-1 && y == gl[Y]-1 && x == gl[X]-ll[X]  ) )
              ASSERT( fread( buffer_pt->next->data, sizeof(double), read_size, fin ) > 0 );
          }
          
          if ( g.my_rank == desired_rank ) {
            MPI_Recv( input_data_pt, read_size, MPI_DOUBLE, 0, k, g.comm_cart, MPI_STATUS_IGNORE );
#ifdef BIG_ENDIAN_CNFG
            for ( i=0; i<read_size; i++ ) {
              byteswap8( (char *) ( input_data_pt + i ) );
            }
#endif
            if ( t == gl[T]-1 ) {
              if ( g.anti_pbc ) {
                for ( j=0; j<read_size; j+=4*18 )
                  for ( i=0; i<18; i++ )
                    (input_data_pt+j+T*18)[i] = -(input_data_pt+j+T*18)[i];
              }
            }
            input_data_pt += read_size;
          }
          
          if ( g.my_rank == 0 ) {
            MPI_Wait( &sreq, MPI_STATUS_IGNORE );
            buffer_pt = buffer_pt->next;
          }
          
          k = (k+1)%10000;
        }
  
  if( g.my_rank == 0)
    fclose( fin );
  
  if ( g.my_rank == 0 ) {
    FREE( buffer[0].data, double, read_size );
    FREE( buffer[1].data, double, read_size );
  }

}



void read_conf_multi( double *input_data, char *input_name, double *conf_plaq, level_struct *l ) {
  
/*********************************************************************************
* Reads in the configuration from multi files.
* - double *input_data: Variable where conf data is stored.
* - char *input_name: Name of the input file.
* - double *conf_plaq: Holds the plaquette of given configuration.                                                
*********************************************************************************/

  int t, z, y, i, j,  mu, lsize[4],
      *gl=l->global_lattice, *ll=l->local_lattice, read_size = 4*18*ll[X];
  double *input_data_pt, plaq;
  FILE* fin = NULL;

  char postfix[100];
  char full_input_name[1000];

  /****************************************************************
   * 2016 Oct. 3:  Issaku Kanamori
   *   it seems MPI rank = g.my_coords[T] * Pz * Py * Px
   *                     g.my_coords[Z] * Py * Px
   *                   g.my_coords[Y] * Px
   *                    g.my_coords[X]
   *   for (Pt,Pz,Py,Px) processor lattice.
   *  ( It might depend on the implementation of MPI.  At least,
   *    I could not find any specification in the manural of MPI_Cart_* )
   *
   *  If you have a problem in opening the file due to misdistribution
   *  of configuration files to computing nodes, try 
   *  #define READ_CONF_MULTI_CHECKFILE
   *  and find out which node is trying to read which file.
   ******************************************************************/
  printf0("reading from multi file configuration\n");
  sprintf(postfix, ".pt%dpz%dpy%dpx%d", 
          g.my_coords[T], g.my_coords[Z], g.my_coords[Y], g.my_coords[X]);
  sprintf(full_input_name,"%s%s", input_name, postfix);
#ifdef READ_CONF_MULTI_CHECKFILE
  fin = fopen( full_input_name, "rb" );
  int open_flag=0;
  int open_flag_sum=0;
  if(fin == NULL){
    printf("ERROR!  node=%d,  filename=%s\n", g.my_rank, full_input_name);
    fflush(0);
    open_flag=1;
  } else {
    printf("opened: node=%d,  filename=%s\n", g.my_rank, full_input_name);
    fflush(0);
  }
  MPI_Allreduce( &open_flag, &open_flag_sum, 1, MPI_INT, MPI_SUM, g.comm_cart );
  ASSERT(open_flag_sum==0);
#else
  ASSERT( (fin = fopen( full_input_name, "rb" )) != NULL );
#endif

  // Read in lattice size in each dimension
  ASSERT( fread( lsize, sizeof(int), 4, fin ) > 0 );
#ifdef BIG_ENDIAN_CNFG
  byteswap( (char *) &lsize[0] );  
  byteswap( (char *) &lsize[1] );
  byteswap( (char *) &lsize[2] );
  byteswap( (char *) &lsize[3] );
#endif

  for ( mu=0; mu<4; mu++ )
    ASSERT( lsize[mu] == gl[mu] );

  // Read in plaquette
  ASSERT( fread( &plaq, sizeof(double), 1, fin ) > 0 );
#ifdef BIG_ENDIAN_CNFG
  byteswap8( (char *) &plaq );
#endif
  printf0("\nDesired average plaquette: %.13lf in [0,3]\n", plaq );
  printf0("\nDesired average plaquette: %.13lf in [0,1]\n", plaq/3.0 );
  *conf_plaq = plaq; 

  // check if value of the plaauette is the same 
  MPI_Bcast(conf_plaq, 1, MPI_DOUBLE, 0, g.comm_cart);
  ASSERT( plaq == *conf_plaq );

  // read the configuration
  input_data_pt = input_data;
  for ( t=0; t<ll[T]; t++ )
    for ( z=0; z<ll[Z]; z++ )
      for ( y=0; y<ll[Y]; y++ ) {
        // read ll[X] data at once  (see def. of read_size)
        ASSERT( fread( input_data_pt, sizeof(double), read_size, fin ) > 0 );
        #ifdef BIG_ENDIAN_CNFG
        for ( i=0; i<read_size; i++ ) {
          byteswap8( (char *) ( input_data_pt + i ) );
        }
        #endif
        if ( (g.my_coords[T]+1 == gl[T]/ll[T]) && (t == ll[T]-1) ) {
          if ( g.anti_pbc ) {
            for ( j=0; j<read_size; j+=4*18 )
              for ( i=0; i<18; i++ )
                (input_data_pt+j+T*18)[i] = -(input_data_pt+j+T*18)[i];
          }
        }
        input_data_pt += read_size;
      }

  fclose( fin );
}


void write_header( FILE **file, double *lambda, char* vector_type, int n, level_struct *l ) {
  
  fprintf( *file, "<header>\n" );
  fprintf( *file, "%s\n", vector_type );
  fprintf( *file, "clifford basis: %s\n", CLIFFORD_BASIS );
  fprintf( *file, "m0: %.14lf\n", l->real_shift );
  fprintf( *file, "csw: %.14lf\n", g.csw );
  fprintf( *file, "clov plaq: %.14lf\n", g.plaq_clov );
  fprintf( *file, "hopp plaq: %.14lf\n", g.plaq_hopp );
  fprintf( *file, "clov conf name: %s\n", g.in_clov );
  fprintf( *file, "hopp conf name: %s\n", g.in );
  fprintf( *file, "X: %d\n", l->global_lattice[X] );
  fprintf( *file, "Y: %d\n", l->global_lattice[Y] );
  fprintf( *file, "Z: %d\n", l->global_lattice[Z] );
  fprintf( *file, "T: %d\n", l->global_lattice[T] );
  fprintf( *file, "X local: %d\n", l->local_lattice[X] );
  fprintf( *file, "Y local: %d\n", l->local_lattice[Y] );
  fprintf( *file, "Z local: %d\n", l->local_lattice[Z] );
  fprintf( *file, "T local: %d\n", l->local_lattice[T] );
  fprintf( *file, "number of vectors: %d\n", n );
  fprintf( *file, "krylov subspace size: %d\n", 100 );
  fprintf( *file, "clifford basis: %s\n", CLIFFORD_BASIS );
  if ( lambda != NULL ) {
    fprintf( *file, "eigenvalues: " );
    for ( int i=0; i<2*n; i++ ) {
      fprintf( *file, "%.16lf ", lambda[i] );
    }
    fprintf( *file, "\n");
  }
  fprintf( *file, "</header>\n" );
}

void vector_io( double *phi, char *filename, const int mode, level_struct *l ) {
  
  int t, z, y, x, *gl=l->global_lattice, *ll=l->local_lattice, bar_size = 24*ll[X], desired_rank;
  double *phi_pt = phi, t0, t1, norm;
  FILE* file = NULL;
  MPI_Request sreq, rreq;
  confbuffer_struct buffer[2];
  confbuffer_struct *buffer_pt = NULL;
  
  buffer[0].next = &(buffer[1]);
  buffer[1].next = &(buffer[0]);  
  buffer[0].data = NULL;
  buffer[1].data = NULL;
  
  if ( g.my_rank == 0 ) {
    MALLOC( buffer[0].data, double, bar_size );
    MALLOC( buffer[1].data, double, bar_size );
    buffer_pt = &buffer[0];
  }
  
  t0 = MPI_Wtime();
  
  if ( mode == _READ ) {
    
    if ( g.my_rank == 0 ) {
      ASSERT( (file = fopen( filename, "rb" )) != NULL );
      char *cur_line = NULL;
      MALLOC( cur_line, char, STRINGLENGTH );
      fgets( cur_line, STRINGLENGTH-1, file );
      if ( strcmp( cur_line, "<header>\n" ) ) {
        fseek( file, 0L, SEEK_SET );
      } else {
        do {
          ASSERT( fgets( cur_line, STRINGLENGTH-1, file ) );
        } while ( strcmp( cur_line, "</header>\n" ) );
      }
      FREE( cur_line, char, STRINGLENGTH );
    }
    
    printf0("reading from file \"%s\" ...\n", filename );
  
    if ( g.my_rank == 0 ) {
      ASSERT( fread( buffer_pt->data, sizeof(double), bar_size, file ) > 0 );
    }
    
    for ( t=0; t<gl[T]; t++ )
      for ( z=0; z<gl[Z]; z++ )
        for ( y=0; y<gl[Y]; y++ )
          for ( x=0; x<gl[X]; x+=ll[X] ) {
            desired_rank = process_index( t, z, y, x, ll );
            
            if ( g.my_rank == 0 ) {
              MPI_Isend( buffer_pt->data, bar_size, MPI_DOUBLE, desired_rank, 0, g.comm_cart, &sreq );
              if ( ! ( t == gl[T]-1 && z == gl[Z]-1 && y == gl[Y]-1 && x == gl[X]-ll[X]  ) )
                ASSERT( fread( buffer_pt->next->data, sizeof(double), bar_size, file ) > 0 );
            }
            
            if ( g.my_rank == desired_rank ) {
              MPI_Recv( phi_pt, bar_size, MPI_DOUBLE, 0, 0, g.comm_cart, MPI_STATUS_IGNORE );
              phi_pt += bar_size;
            }
            
            if ( g.my_rank == 0 ) {
              MPI_Wait( &sreq, MPI_STATUS_IGNORE );
              buffer_pt = buffer_pt->next;
            }
          }
          
  } else if ( mode == _WRITE ) {
    if ( g.my_rank == 0 ) {
      ASSERT( (file = fopen( filename, "wb" )) != NULL );
      write_header( &file, NULL, filename, 1, l );
    }
    printf0("writing file \"%s\" ...\n", filename );
    
    for ( t=0; t<gl[T]; t++ ) {
      for ( z=0; z<gl[Z]; z++ )
        for ( y=0; y<gl[Y]; y++ )
          for ( x=0; x<gl[X]; x+=ll[X] ) {
            desired_rank = process_index( t, z, y, x, ll );
            
            if ( g.my_rank == desired_rank ) {
              MPI_Isend( phi_pt, bar_size, MPI_DOUBLE, 0, 0, g.comm_cart, &sreq );
              phi_pt += bar_size;
            }
            
            if ( g.my_rank == 0 ) {
              MPI_Irecv( buffer_pt->next->data, bar_size, MPI_DOUBLE, desired_rank, 0, g.comm_cart, &rreq );
              
              if ( ! ( t == 0 && z == 0 && y == 0 && x == 0  ) ) {
                fwrite( buffer_pt->data, sizeof(double), bar_size, file );
              }
              
              MPI_Wait( &rreq, MPI_STATUS_IGNORE );
              buffer_pt = buffer_pt->next;
            }
            
            if ( g.my_rank == desired_rank ) {
              MPI_Wait( &sreq, MPI_STATUS_IGNORE );
            }
          }
    }
          
    if ( g.my_rank == 0 ) {
      fwrite( buffer_pt->data, sizeof(double), bar_size, file );
    }
  
  } else
    ASSERT( mode == _READ || mode == _WRITE );
  
  if ( g.my_rank == 0 ){
    fclose( file );
  }

  t1 = MPI_Wtime();
  
  if ( g.my_rank == 0 ) {
    FREE( buffer[0].data, double, bar_size );
    FREE( buffer[1].data, double, bar_size );
  }
  
  norm = global_norm_double( (vector_double)phi, 0, l->inner_vector_size, l );
  printf0("norm: %e\n", norm );
  printf0("...done (%lf seconds)\n\n", t1-t0 ); 
}


void vector_io_single_file( double *psi, double *lambda, char *filename, const int mode, int n, char *vector_type, level_struct *l ) {
  
  int t, z, y, x, *gl=l->global_lattice, *ll=l->local_lattice, bar_size = 24*ll[X], desired_rank, j;
  double t0, t1;
  double *phi_pt = NULL;
  double *phi = NULL;
  FILE* file = NULL;
  MPI_Request sreq, rreq;
  confbuffer_struct buffer[2];
  confbuffer_struct *buffer_pt = NULL;
    
  buffer[0].next = &(buffer[1]);
  buffer[1].next = &(buffer[0]);  
  buffer[0].data = NULL;
  buffer[1].data = NULL;

  if ( g.my_rank == 0 ) {
    MALLOC( buffer[0].data, double, bar_size );
    MALLOC( buffer[1].data, double, bar_size );
    buffer_pt = &buffer[0];
  }
  
  t0 = MPI_Wtime();
  
  if ( mode == _READ ) {

    if ( g.my_rank == 0 ) {
      ASSERT( (file = fopen( filename, "rb" )) != NULL );
      printf0("reading from file \"%s\" ...\n", filename );

      char *cur_line = NULL;
      MALLOC( cur_line, char, STRINGLENGTH );
      do {
        ASSERT( fgets( cur_line, STRINGLENGTH-1, file ) );
      } while ( strcmp( cur_line, "</header>\n" ) );
      FREE( cur_line, char, STRINGLENGTH );
    }
    
    for ( j=0; j<n; j++ ) {
      if ( g.my_rank == 0 ) {
        ASSERT( fread( buffer_pt->data, sizeof(double), bar_size, file ) );
      }

      phi=(double *) (l->x);
      phi_pt=phi;
      for ( t=0; t<gl[T]; t++ )
        for ( z=0; z<gl[Z]; z++ )
          for ( y=0; y<gl[Y]; y++ )
            for ( x=0; x<gl[X]; x+=ll[X] ) {
              
              desired_rank = process_index( t, z, y, x, ll );
              
              if ( g.my_rank == desired_rank ) {
                MPI_Irecv( phi_pt, bar_size, MPI_DOUBLE, 0, 0, g.comm_cart, &rreq );
              }
              if ( g.my_rank == 0 ) {
                MPI_Isend( buffer_pt->data, bar_size, MPI_DOUBLE, desired_rank, 0, g.comm_cart, &sreq );
              }
              
              if ( g.my_rank == 0 ) {
                if ( ! ( t == gl[T]-1 && z == gl[Z]-1 && y == gl[Y]-1 && x == gl[X]-ll[X] ) )
                  ASSERT( fread( buffer_pt->next->data, sizeof(double), bar_size, file ) );
              }
              
              if ( g.my_rank == desired_rank ) {
                MPI_Wait( &rreq, MPI_STATUS_IGNORE );
                phi_pt += bar_size;
              }
              if ( g.my_rank == 0 ) {
                MPI_Wait( &sreq, MPI_STATUS_IGNORE );
                buffer_pt = buffer_pt->next;
              }
            }
      if ( psi == NULL ) {
        if ( g.mixed_precision )
          trans_float(l->is_float.test_vector[j], l->x, l->s_float.op.translation_table, l);
        else
          trans_double(l->is_double.test_vector[j], l->x, l->s_double.op.translation_table, l);
      } else {
        vector_double_copy( ((vector_double)psi)+j*l->inner_vector_size, l->x, 0, l->inner_vector_size, l );
      }
    }
  } else if ( mode == _WRITE ) {
    
    if ( g.my_rank == 0 ) {
      ASSERT( (file = fopen( filename, "wb" )) != NULL );
      write_header( &file, lambda, vector_type, n, l );
    }

    printf0("writing file \"%s\" ...\n", filename );
    
    for ( j=0; j<n; j++ ){
      if ( psi == NULL ) {
        if ( g.mixed_precision )
          trans_back_float( l->x, l->is_float.test_vector[j], l->s_float.op.translation_table, l );
        else
          trans_back_double( l->x, l->is_double.test_vector[j], l->s_double.op.translation_table, l );
      } else {
        vector_double_copy( l->x, ((complex_double*)psi)+j*l->inner_vector_size, 0, l->inner_vector_size, l );
      }
      phi=(double *)(l->x);
      phi_pt=phi;
      for ( t=0; t<gl[T]; t++ )
        for ( z=0; z<gl[Z]; z++ )
          for ( y=0; y<gl[Y]; y++ )
            for ( x=0; x<gl[X]; x+=ll[X] ) {
              
              desired_rank = process_index( t, z, y, x, ll );
              
              if ( g.my_rank == 0 ) {
                MPI_Irecv( buffer_pt->next->data, bar_size, MPI_DOUBLE, desired_rank, 0, g.comm_cart, &rreq );
              }
              if ( g.my_rank == desired_rank ) {
                MPI_Isend( phi_pt, bar_size, MPI_DOUBLE, 0, 0, g.comm_cart, &sreq );
              }
              
              if ( g.my_rank == 0 ) {
                if ( ! ( t == 0 && z == 0 && y == 0 && x == 0  ) ) {
                  fwrite( buffer_pt->data, sizeof(double), bar_size, file );
                }
              }
              
              if ( g.my_rank == 0 ) {
                MPI_Wait( &rreq, MPI_STATUS_IGNORE );
                buffer_pt = buffer_pt->next;
              }
              if ( g.my_rank == desired_rank ) {
                MPI_Wait( &sreq, MPI_STATUS_IGNORE );
                phi_pt += bar_size;
              }
              
            }
            
      if ( g.my_rank == 0 ) {
        fwrite( buffer_pt->data, sizeof(double), bar_size, file );
      }
    }
  } else
    ASSERT( mode == _READ || mode == _WRITE );

  if ( g.my_rank == 0 ) {
    fclose( file );
  }
  
  t1 = MPI_Wtime();
  
  if ( g.my_rank == 0 ) {
    FREE( buffer[0].data, double, bar_size );
    FREE( buffer[1].data, double, bar_size );
  }
  
  printf0("...done (%lf seconds)\n\n", t1-t0 ); 
}


void d_dump( config_double D, level_struct *l ) {
  
  if ( g.num_processes == 1 ) {
    
    int i, x, y, z, t;
    
    i=0;
    for ( t=0; t<l->local_lattice[T]; t++ )
      for ( z=0; z<l->local_lattice[Z]; z++ )
        for ( y=0; y<l->local_lattice[Y]; y++ )
          for ( x=0; x<l->local_lattice[X]; x++ ) {
            printf("site_number=%d, t=%d, z=%d, y=%d, x=%d.\n", i/36, t, z, y, x );
            printf("dir: T\n");
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+0]), cimag(D[i+0]), creal(D[i+1]), cimag(D[i+1]), creal(D[i+2]), cimag(D[i+2]));
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+3]), cimag(D[i+3]), creal(D[i+4]), cimag(D[i+4]), creal(D[i+5]), cimag(D[i+5]));
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+6]), cimag(D[i+6]), creal(D[i+7]), cimag(D[i+7]), creal(D[i+8]), cimag(D[i+8]));
            i+=9;
            printf("dir: Z\n");
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+0]), cimag(D[i+0]), creal(D[i+1]), cimag(D[i+1]), creal(D[i+2]), cimag(D[i+2]));
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+3]), cimag(D[i+3]), creal(D[i+4]), cimag(D[i+4]), creal(D[i+5]), cimag(D[i+5]));
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+6]), cimag(D[i+6]), creal(D[i+7]), cimag(D[i+7]), creal(D[i+8]), cimag(D[i+8]));
            i+=9;
            printf("dir: Y\n");
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+0]), cimag(D[i+0]), creal(D[i+1]), cimag(D[i+1]), creal(D[i+2]), cimag(D[i+2]));
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+3]), cimag(D[i+3]), creal(D[i+4]), cimag(D[i+4]), creal(D[i+5]), cimag(D[i+5]));
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+6]), cimag(D[i+6]), creal(D[i+7]), cimag(D[i+7]), creal(D[i+8]), cimag(D[i+8]));
            i+=9;
            printf("dir: X\n");
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+0]), cimag(D[i+0]), creal(D[i+1]), cimag(D[i+1]), creal(D[i+2]), cimag(D[i+2]));
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+3]), cimag(D[i+3]), creal(D[i+4]), cimag(D[i+4]), creal(D[i+5]), cimag(D[i+5]));
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+6]), cimag(D[i+6]), creal(D[i+7]), cimag(D[i+7]), creal(D[i+8]), cimag(D[i+8]));
            i+=9;
            printf("------------------------------------------------------------\n");
          }
  }
}

