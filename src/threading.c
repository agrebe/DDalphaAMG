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

#include "main.h" // for level_struct
#include <omp.h>


void no_barrier(int id)
{
}
void no_hyperthread_barrier(void *barrier, int id)
{
}
void core_barrier(int core)
{
#pragma omp barrier
}
void hyperthread_barrier(void *barrier, int hyperthead)
{
    // no hyperthreads for now => do nothing
}


void init_common_thread_data(struct common_thread_data *common)
{
    common->barrier = &core_barrier;
    common->thread_barrier = &hyperthread_barrier;
}


void setup_threading(struct Thread *threading, struct common_thread_data *common, struct level_struct *l)
{
    // no hyperthreading for now
    setup_threading_external(threading, common, l, omp_get_num_threads(), 1, omp_get_thread_num(), 0);
}


void setup_threading_external(struct Thread *threading, struct common_thread_data *common, struct level_struct *l,
        int n_core, int n_thread, int core, int thread)
{
    threading->n_core = n_core;
    threading->n_thread = n_thread;
    threading->core = core;
    threading->thread = thread;
    threading->thread_barrier_data = 0;

    update_threading(threading, l);

    threading->barrier = common->barrier;
    threading->thread_barrier = common->thread_barrier;

}


void update_threading(struct Thread *threading, struct level_struct *l)
{
    struct level_struct *current = l;

    while(1)
    {
        *(threading->start_site+current->depth) = 0;
        *(threading->end_site+current->depth) = current->num_inner_lattice_sites;
        threading->n_site[current->depth] = threading->end_site[current->depth] - threading->start_site[current->depth];
        threading->start_index[current->depth] = threading->start_site[current->depth]*current->num_lattice_site_var;
        threading->end_index[current->depth]   =   threading->end_site[current->depth]*current->num_lattice_site_var;
        threading->n_index[current->depth]     =     threading->n_site[current->depth]*current->num_lattice_site_var;

        if(current->next_level == NULL)
            break;
        current = current->next_level;
    }
}


void setup_no_threading(struct Thread *no_threading, struct level_struct *l)
{
    no_threading->core = 0;
    no_threading->n_core = 1;
    // no hyperthreading for now
    no_threading->thread = 0;
    no_threading->n_thread = 1;

    update_threading(no_threading, l);

    no_threading->barrier = &no_barrier;
    no_threading->thread_barrier = &no_hyperthread_barrier;
    
}

void finalize_common_thread_data( struct common_thread_data *common ) {
}


void finalize_no_threading( struct Thread *no_threading ) {
}
