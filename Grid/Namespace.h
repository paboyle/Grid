/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/Namespace.h

Copyright (C) 2016

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
#pragma once

#include <type_traits>
#include <cassert>

#define NAMESPACE_BEGIN(A) namespace A {
#define NAMESPACE_END(A)   }
#define GRID_NAMESPACE_BEGIN NAMESPACE_BEGIN(Grid)
#define GRID_NAMESPACE_END   NAMESPACE_END(Grid)
#define NAMESPACE_CHECK(x) struct namespaceTEST##x {};  static_assert(std::is_same<namespaceTEST##x, ::namespaceTEST##x>::value,"Not in :: at"  ); 
