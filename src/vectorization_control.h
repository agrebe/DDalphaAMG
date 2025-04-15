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

#ifndef VECTORIZATION_CONTROL_H
#define VECTORIZATION_CONTROL_H

#define OPERATOR_COMPONENT_OFFSET_float  (SIMD_LENGTH_float *((l->num_eig_vect+SIMD_LENGTH_float -1)/SIMD_LENGTH_float ))
#define OPERATOR_COMPONENT_OFFSET_double (SIMD_LENGTH_double*((l->num_eig_vect+SIMD_LENGTH_double-1)/SIMD_LENGTH_double))

#define OPERATOR_TYPE_float float
#define OPERATOR_TYPE_double double

#endif // VECTORIZATION_CONTROL_H
