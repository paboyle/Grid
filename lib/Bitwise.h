/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/Bitwise.h

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
#ifndef GRID_BITWISE_H
#define GRID_BITWISE_H

#include <cassert>
#include <bitset>
#include <climits>
#include <Config.h>

#ifdef GRID_DEFAULT_PRECISION_SINGLE
#define GRID_REAL_BYTES 4
#endif
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
#define GRID_REAL_BYTES 8
#endif


namespace Grid {

void show_binaryrep(const unsigned char* a, size_t size);

template <typename T>
void show_binaryrep(const T& a) {
  const char* beg = reinterpret_cast<const char*>(&a);
  const char* end = beg + sizeof(a);
  unsigned int ctr = 0;
  while (beg != end) {
  	std::cout << std::bitset<CHAR_BIT>(*beg++) << ' ';
  	ctr++;
  	if (ctr % GRID_REAL_BYTES == 0) std::cout << '\n';
  }
  std::cout << '\n';
}

template <typename T>
void bitwise_xor(T& l, T& r, unsigned char* xors) {
  assert(sizeof(l) == sizeof(r));
  unsigned char* org = reinterpret_cast<unsigned char*>(&l);
  unsigned char* cur = reinterpret_cast<unsigned char*>(&r);
  int words = sizeof(l) / sizeof(*org);
  unsigned char result = 0;
  for (int w = 0; w < words; w++) xors[w] = (org[w] ^ cur[w]);
}

}; // namespace 


#endif
