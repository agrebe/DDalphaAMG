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

void inv_iter_2lvl_extension_setup_PRECISION( int setup_iter, level_struct *l );
void inv_iter_inv_fcycle_PRECISION( int setup_iter, level_struct *l );
void read_tv_from_file_PRECISION( level_struct *l );

void coarse_grid_correction_PRECISION_setup( level_struct *l ) {
  
  if ( !l->idle ) {
    
    coarse_operator_PRECISION_alloc( l );
    coarse_operator_PRECISION_setup( l->is_PRECISION.interpolation, l );
    
    if ( !l->next_level->idle ) {
      if ( l->next_level->level > 0 ) {
        schwarz_PRECISION_alloc( &(l->next_level->s_PRECISION), l->next_level );
        schwarz_layout_PRECISION_define( &(l->next_level->s_PRECISION), l->next_level );
      } else {
        operator_PRECISION_alloc( &(l->next_level->s_PRECISION.op), _ORDINARY, l->next_level );
        operator_PRECISION_define( &(l->next_level->s_PRECISION.op), l->next_level );
        interpolation_PRECISION_alloc( l->next_level );
      }
    } else {
      interpolation_PRECISION_dummy_alloc( l->next_level );
    }
    
    conf_PRECISION_gather( &(l->next_level->s_PRECISION.op), &(l->next_level->op_PRECISION), l->next_level );
    
    
    if ( !l->next_level->idle && l->next_level->level > 0 ) {
      schwarz_PRECISION_boundary_update( &(l->next_level->s_PRECISION), l->next_level );
      if ( g.method >= 4 && g.odd_even ) {
        coarse_oddeven_setup_PRECISION( &(l->next_level->s_PRECISION.op), _REORDER, l->next_level );
      }
      l->next_level->p_PRECISION.op = &(l->next_level->s_PRECISION.op);
    }
    if ( !l->next_level->idle && l->next_level->level == 0 && g.odd_even ) {
      coarse_oddeven_setup_PRECISION( &(l->next_level->s_PRECISION.op), _NO_REORDERING, l->next_level );
    } 
  }
  
  if ( l->next_level->level > 0 ) {
    next_level_setup( NULL, l->next_level );
    if ( !l->next_level->idle )
      interpolation_PRECISION_alloc( l->next_level );
    if ( !l->idle ) {
      for ( int i=0; i<MIN(l->next_level->num_eig_vect,l->num_eig_vect); i++ ) {
        restrict_PRECISION( l->next_level->is_PRECISION.test_vector[i], l->is_PRECISION.test_vector[i], l );
      }
      for ( int i=MIN(l->next_level->num_eig_vect,l->num_eig_vect); i<l->next_level->num_eig_vect; i++ ) {
        if ( !l->next_level->idle )
          vector_PRECISION_define_random( l->next_level->is_PRECISION.test_vector[i], 0,
                                          l->next_level->inner_vector_size, l->next_level );
      }
    }
    if ( !l->next_level->idle )
      interpolation_PRECISION_define( NULL, l->next_level );
    
    coarse_grid_correction_PRECISION_setup( l->next_level );
    if ( !l->next_level->idle )
      define_interpolation_PRECISION_operator( l->next_level->is_PRECISION.interpolation, l->next_level );
  }
}


void iterative_PRECISION_setup( int setup_iter, level_struct *l ) {
  if ( l->depth == 0 ) {
    switch ( g.interpolation ) {
      case 2: inv_iter_inv_fcycle_PRECISION( setup_iter, l ); break;
      case 3: inv_iter_inv_fcycle_PRECISION( setup_iter, l ); break;
      case 4: read_tv_from_file_PRECISION( l ); break;
      default: inv_iter_2lvl_extension_setup_PRECISION( setup_iter, l ); break;
    }
  }

  level_struct *lp = l;
  while( lp->level > 0 ) {
    lp = lp->next_level;
    if ( lp == NULL )
      break;
  }
}


void read_tv_from_file_PRECISION( level_struct *l ) {
  
  if ( l->depth == 0 ) {
    if ( g.tv_io_single_file ) {
      vector_io_single_file( NULL, NULL, g.tv_io_file_name, _READ, l->num_eig_vect, "test vectors", l );
      re_setup_PRECISION( l );
    } else {

      int n = l->num_eig_vect, i;
      char filename[STRINGLENGTH+1];
      vector_double tmp = NULL;
      
      MALLOC( tmp, complex_double, l->inner_vector_size );
      
      for ( i=0; i<n; i++ ) {
        sprintf( filename, "%s.%02d", g.tv_io_file_name, i );
        printf0("%s.%02d\n", g.tv_io_file_name, i );
        vector_io( (double*)tmp, filename, _READ, l );
        trans_PRECISION( l->is_PRECISION.test_vector[i], tmp, l->s_PRECISION.op.translation_table, l );
      }
      
      FREE( tmp, complex_double, l->inner_vector_size );


      re_setup_PRECISION( l );
    }
  }
}


void coarse_grid_correction_PRECISION_free( level_struct *l ) {
  
  next_level_free( l->next_level );
  
  if ( !l->idle ) {
    if ( !l->next_level->idle ) {
      if ( l->next_level->level > 0 ) {
        schwarz_PRECISION_free( &(l->next_level->s_PRECISION), l->next_level );
        if ( g.method >= 4 && g.odd_even ) {
          coarse_oddeven_free_PRECISION( l->next_level );
        }
      } else {
        operator_PRECISION_free( &(l->next_level->s_PRECISION.op), _ORDINARY, l->next_level );
        interpolation_PRECISION_free( l->next_level );
        if ( g.odd_even )
          coarse_oddeven_free_PRECISION( l->next_level );
      }
    } else {
      interpolation_PRECISION_dummy_free( l->next_level );
    }
    interpolation_PRECISION_free( l );
    coarse_operator_PRECISION_free( l );
  }  
}


void interpolation_PRECISION_define( vector_double *V, level_struct *l ) {
  
  int k, i, n = l->num_eig_vect,
      pc = 0, pi = 1, pn = n*6;
  vector_PRECISION *buffer = NULL;
  int start = 0;
  int end   = l->num_inner_lattice_sites * l->num_lattice_site_var;
    
  if ( V == NULL ) {
    
    PUBLIC_MALLOC( buffer, complex_PRECISION*, 3 );
    buffer[0] = NULL;
    PUBLIC_MALLOC( buffer[0], complex_PRECISION, l->vector_size*3 );
    
    for( i=1; i<3; i++)
      buffer[i] = buffer[0] + l->vector_size*i;
    if ( g.print > 0 ) printf0("initial definition --- depth: %d\n", l->depth );
    if ( g.print > 0 ) { printf0("\033[0;42m\033[1;37m|"); fflush(0); }
    

    for ( k=0; k<n; k++ ) {
//       if ( l->depth == 0 ) {
        vector_PRECISION_define_random( l->is_PRECISION.test_vector[k], 0, l->inner_vector_size, l );
//       }
      
      smoother_PRECISION( buffer[0], NULL, l->is_PRECISION.test_vector[k],
                          1, _NO_RES, _NO_SHIFT, l );
      vector_PRECISION_copy( l->is_PRECISION.test_vector[k], buffer[0], start, end, l );
      smoother_PRECISION( buffer[0], NULL, l->is_PRECISION.test_vector[k],
                          g.method>=4?1:2, _NO_RES, _NO_SHIFT, l );
      vector_PRECISION_copy( l->is_PRECISION.test_vector[k], buffer[0], start, end, l );
      smoother_PRECISION( buffer[0], NULL, l->is_PRECISION.test_vector[k],
                          g.method>=4?1:3, _NO_RES, _NO_SHIFT, l );
      vector_PRECISION_copy( l->is_PRECISION.test_vector[k], buffer[0], start, end, l );
        
      pc += 6;
      if ( pc >= 0.2*pi*pn ) { if ( g.print > 0 ) printf0("%4d%% |", 20*pi); if ( g.my_rank == 0 ) fflush(0); pi++; }
    }
    
    PUBLIC_FREE( buffer[0], complex_PRECISION, l->vector_size*3 );
    PUBLIC_FREE( buffer, complex_PRECISION*, 3 );
    
    for ( k=0; k<n; k++ ) {
      vector_PRECISION_real_scale( l->is_PRECISION.test_vector[k], l->is_PRECISION.test_vector[k],
                                  1.0/global_norm_PRECISION( l->is_PRECISION.test_vector[k], 0, l->inner_vector_size, l ),
                                  start, end, l );
    }
    
    if ( g.print > 0 ) printf0("\033[0m\n");
    
    } else {
    for ( i=0; i<n; i++ ) {
      trans_PRECISION( l->is_PRECISION.test_vector[i], V[i], l->s_PRECISION.op.translation_table, l );
    }
  }

  for ( k=0; k<n; k++ ) {
    vector_PRECISION_copy( l->is_PRECISION.interpolation[k], l->is_PRECISION.test_vector[k], start, end, l );
  }
  
  
    

  gram_schmidt_on_aggregates_PRECISION( l->is_PRECISION.interpolation, n, l );
  
}


void re_setup_PRECISION( level_struct *l ) {
  
  if ( l->level > 0 ) {
    if ( !l->idle ) {
      for ( int i=0; i<l->num_eig_vect; i++ ) {
        vector_PRECISION_copy( l->is_PRECISION.interpolation[i], l->is_PRECISION.test_vector[i],
            0, l->num_inner_lattice_sites * l->num_lattice_site_var, l );
      }
      gram_schmidt_on_aggregates_PRECISION( l->is_PRECISION.interpolation, l->num_eig_vect, l );
      coarse_operator_PRECISION_setup( l->is_PRECISION.interpolation, l );
      define_interpolation_PRECISION_operator( l->is_PRECISION.interpolation, l );
      conf_PRECISION_gather( &(l->next_level->s_PRECISION.op), &(l->next_level->op_PRECISION), l->next_level );
      if ( !l->next_level->idle && l->next_level->level > 0 ) {
        schwarz_PRECISION_boundary_update( &(l->next_level->s_PRECISION), l->next_level );
        if ( g.method >= 4 && g.odd_even ) {
          coarse_oddeven_setup_PRECISION_set_couplings( &(l->next_level->s_PRECISION.op), _REORDER, l->next_level );
        }
      }
      if ( !l->next_level->idle && l->next_level->level == 0 && g.odd_even ) {
        coarse_oddeven_setup_PRECISION_set_couplings( &(l->next_level->s_PRECISION.op), _NO_REORDERING, l->next_level );
      }
      re_setup_PRECISION( l->next_level );
    }
  }  
}


void inv_iter_2lvl_extension_setup_PRECISION( int setup_iter, level_struct *l ) {
  
  if ( !l->idle ) {
    vector_PRECISION buf1 = NULL;
    gmres_PRECISION_struct gmres;
    
    // TODO: bugfix - threading, etc
    
    MALLOC( buf1, complex_PRECISION, l->vector_size );
    fgmres_PRECISION_struct_init( &gmres );
    fgmres_PRECISION_struct_alloc( g.coarse_iter, g.coarse_restart, l->next_level->vector_size, g.coarse_tol, 
                                   _COARSE_GMRES, _NOTHING, NULL, apply_coarse_operator_PRECISION, &gmres, l->next_level );
    
    if ( g.odd_even && l->next_level->level == 0 )
      gmres.v_end = l->next_level->oe_op_PRECISION.num_even_sites*l->next_level->num_lattice_site_var;
    
    for ( int k=0; k<setup_iter; k++ ) {
      int pc = 0, pi = 1, pn = l->num_eig_vect*l->post_smooth_iter;
      printf0("depth: %d, 2lvl correction step number %d...\n", l->depth, k+1 ); 
      printf0("\033[0;42m\033[1;37m|"); fflush(0);
      for ( int i=0; i<l->num_eig_vect; i++ ) {
        restrict_PRECISION( gmres.b, l->is_PRECISION.test_vector[i], l );
        if ( !l->next_level->idle ) {
          if ( g.odd_even && l->next_level->level == 0 ) {
            coarse_solve_odd_even_PRECISION( &gmres, &(l->next_level->oe_op_PRECISION), l->next_level );
          } else {
            fgmres_PRECISION( &gmres, l->next_level );
          }
        }
        interpolate3_PRECISION( buf1, gmres.x, l );
        smoother_PRECISION( buf1, NULL, l->is_PRECISION.test_vector[i], l->post_smooth_iter, _RES, _NO_SHIFT, l );
        vector_PRECISION_real_scale( l->is_PRECISION.test_vector[i], buf1,
                                     1.0/global_norm_PRECISION( buf1, 0, l->inner_vector_size, l ),
                                     0, l->num_inner_lattice_sites * l->num_lattice_site_var, l );
        pc += l->post_smooth_iter;
        if ( pc >= 0.2*pi*pn ) { printf0("%4d%% |", 20*pi); fflush(0); pi++; }
      }
      printf0("\033[0m\n");
      
      for ( int i=0; i<l->num_eig_vect; i++ )
        vector_PRECISION_copy( l->is_PRECISION.interpolation[i], l->is_PRECISION.test_vector[i],
            0, l->num_inner_lattice_sites * l->num_lattice_site_var, l );
      gram_schmidt_on_aggregates_PRECISION( l->is_PRECISION.interpolation, l->num_eig_vect, l );
      coarse_operator_PRECISION_setup( l->is_PRECISION.interpolation, l );
      define_interpolation_PRECISION_operator( l->is_PRECISION.interpolation, l );
      conf_PRECISION_gather( &(l->next_level->s_PRECISION.op), &(l->next_level->op_PRECISION), l->next_level );
      if ( !l->next_level->idle && l->next_level->level > 0 ) {
        schwarz_PRECISION_boundary_update( &(l->next_level->s_PRECISION), l->next_level );
        if ( g.method >= 4 && g.odd_even ) {
          coarse_oddeven_setup_PRECISION_set_couplings( &(l->next_level->s_PRECISION.op), _REORDER, l->next_level );
        }
      }
      if ( !l->next_level->idle && l->next_level->level == 0 && g.odd_even ) {
        coarse_oddeven_setup_PRECISION_set_couplings( &(l->next_level->s_PRECISION.op), _NO_REORDERING, l->next_level );
      } 
    }
    
    if ( l->level > 1 )
      inv_iter_2lvl_extension_setup_PRECISION( setup_iter, l->next_level );

    FREE( buf1, complex_PRECISION, l->vector_size );
    fgmres_PRECISION_struct_free( &gmres, l );
  }
}


void set_kcycle_tol_PRECISION( PRECISION tol, level_struct *l ) {
  
  if ( !l->idle )
    l->p_PRECISION.tol = tol;
  
  if ( l->level > 1 )
    set_kcycle_tol_PRECISION( tol, l->next_level );
}


void test_vector_PRECISION_update( int i, level_struct *l ) {
  
  if ( l->level > 1 )
    test_vector_PRECISION_update( i, l->next_level );
  
  if ( !l->idle )
    vector_PRECISION_real_scale( l->is_PRECISION.test_vector[i], l->p_PRECISION.x,
                                 1.0/global_norm_PRECISION( l->p_PRECISION.x, 0, l->inner_vector_size, l ),
                                 0, l->num_inner_lattice_sites * l->num_lattice_site_var, l );
}


void inv_iter_inv_fcycle_PRECISION( int setup_iter, level_struct *l ) {
  
  vector_PRECISION v_buf = NULL;
  complex_PRECISION *buffer = NULL;
  
  PUBLIC_MALLOC( buffer, complex_PRECISION, 2*l->num_eig_vect );
      
  if ( l->depth == 0 )
    set_kcycle_tol_PRECISION( g.coarse_tol, l );
  
  PUBLIC_MALLOC( v_buf, complex_PRECISION, l->vector_size );
  
  if ( !l->idle ) {
    for ( int j=0; j<setup_iter; j++ ) {
      int pc = 0, pi = 1, pn = l->num_eig_vect*l->post_smooth_iter;
      
      if ( g.print > 0 ) printf0("depth: %d, bootstrap step number %d...\n", l->depth, j+1 );
      if ( g.print > 0 ) { printf0("\033[0;42m\033[1;37m|"); if ( g.my_rank == 0 ) fflush(0); }
      
      gram_schmidt_PRECISION( l->is_PRECISION.test_vector, buffer, 0, l->num_eig_vect, l );
      
      for ( int i=0; i<l->num_eig_vect; i++ ) {
        vcycle_PRECISION( l->p_PRECISION.x, NULL, l->is_PRECISION.test_vector[i], _NO_RES, l );
        
        test_vector_PRECISION_update( i, l );
        
        pc += l->post_smooth_iter;
        if ( pc >= (int)((0.2*pi)*pn) ) { if ( g.print > 0 ) { printf0("%4d%% |", 20*pi); if ( g.my_rank == 0 ) fflush(0); } pi++; }
      }

      if ( g.print > 0 ) printf0("\033[0m\n");
      
      re_setup_PRECISION( l );
      
      if ( l->depth == 0 && l->next_level->level > 0 ) {
        inv_iter_inv_fcycle_PRECISION( MAX(1,round( ((double)(j+1)*l->next_level->setup_iter)/
        ((double)setup_iter) )), l->next_level );
      }
    }
    if ( l->depth > 0 && l->next_level->level > 0 ) {
      inv_iter_inv_fcycle_PRECISION( MAX(1, round((double)(l->next_level->setup_iter*setup_iter)/
      ((double)l->setup_iter))), l->next_level );
    }
  }
  
  PUBLIC_FREE( v_buf, complex_PRECISION, l->vector_size );
  PUBLIC_FREE( buffer, complex_PRECISION, 2*l->num_eig_vect );
  
  if ( l->depth == 0 ) {
    set_kcycle_tol_PRECISION( g.kcycle_tol, l );
  }
}


