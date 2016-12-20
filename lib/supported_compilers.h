    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/supported_compilers.h

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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */



#ifndef COMPILER_CHECK_H
#define COMPILER_CHECK_H

// exclude unsupported compilers
#if defined(__clang__)
    #define CLANG_VERSION (__clang_major__ * 10000 + __clang_minor__ * 100 + __clang_patchlevel__)
    #if CLANG_VERSION < 30800
        #error "unsupported Clang version"
    #endif
#elif defined(__GNUC__)
    #define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
    #if GCC_VERSION < 40800
        #error "unsupported GCC version"
    #endif
#endif


#endif  // COMPILER_CHECK_H
