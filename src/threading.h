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

#ifndef THREADING_H
#define THREADING_H


#ifdef OPENMP
#include <omp.h>
#else
static inline int omp_get_thread_num( void ) {
  return 0;
}
static inline int omp_get_num_threads( void ) {
  return 1;
}
#endif

struct level_struct;

struct common_thread_data
{
    // barrier among cores
    void (*barrier)(int);
    // barrier among hyperthreads on a core
    void (*thread_barrier)(void *, int);
};

void init_common_thread_data(struct common_thread_data *common);

// holds information relevant for specific core/thread
typedef struct Thread
{
} Thread;

#endif // THREADING_H
