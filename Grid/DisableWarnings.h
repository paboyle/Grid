/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/DisableWarnings.h

Copyright (C) 2016

Author: Guido Cossu <guido.cossu@ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */

#ifndef DISABLE_WARNINGS_H
#define DISABLE_WARNINGS_H



#if defined __GNUC__ && __GNUC__>=6
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

 //disables and intel compiler specific warning (in json.hpp)
#ifdef __ICC
#pragma warning disable 488  
#endif

#ifdef __NVCC__
 //disables nvcc specific warning in json.hpp
#pragma clang diagnostic ignored "-Wdeprecated-register"

#if (__CUDACC_VER_MAJOR__ >= 11) && (__CUDACC_VER_MINOR__ >= 5)
 //disables nvcc specific warning in json.hpp
#pragma nv_diag_suppress unsigned_compare_with_zero
#pragma nv_diag_suppress cast_to_qualified_type
 //disables nvcc specific warning in many files
#pragma nv_diag_suppress esa_on_defaulted_function_ignored
#pragma nv_diag_suppress extra_semicolon
#else
 //disables nvcc specific warning in json.hpp
#pragma diag_suppress unsigned_compare_with_zero
#pragma diag_suppress cast_to_qualified_type
 //disables nvcc specific warning in many files
#pragma diag_suppress esa_on_defaulted_function_ignored
#pragma diag_suppress extra_semicolon
#endif
#endif

// Disable vectorisation in Eigen on the Power8/9 and PowerPC
#ifdef  __ALTIVEC__
#define  EIGEN_DONT_VECTORIZE
#endif
#ifdef  __VSX__
#define  EIGEN_DONT_VECTORIZE
#endif

#endif
