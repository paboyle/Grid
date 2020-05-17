/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/MemoryStats.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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

    See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#pragma once


NAMESPACE_BEGIN(Grid);

std::string sizeString(size_t bytes);

struct MemoryStats
{
  size_t totalAllocated{0}, maxAllocated{0}, 
    currentlyAllocated{0}, totalFreed{0};
};
    
class MemoryProfiler
{
public:
  static MemoryStats *stats;
  static bool        debug;
};

#define memString(bytes) std::to_string(bytes) + " (" + sizeString(bytes) + ")"
#define profilerDebugPrint						\
  if (MemoryProfiler::stats)						\
    {									\
      auto s = MemoryProfiler::stats;					\
      std::cout << GridLogDebug << "[Memory debug] Stats " << MemoryProfiler::stats << std::endl; \
      std::cout << GridLogDebug << "[Memory debug] total  : " << memString(s->totalAllocated) \
		<< std::endl;						\
      std::cout << GridLogDebug << "[Memory debug] max    : " << memString(s->maxAllocated) \
		<< std::endl;						\
      std::cout << GridLogDebug << "[Memory debug] current: " << memString(s->currentlyAllocated) \
		<< std::endl;						\
      std::cout << GridLogDebug << "[Memory debug] freed  : " << memString(s->totalFreed) \
		<< std::endl;						\
    }

#define profilerAllocate(bytes)						\
  if (MemoryProfiler::stats)						\
    {									\
      auto s = MemoryProfiler::stats;					\
      s->totalAllocated     += (bytes);					\
      s->currentlyAllocated += (bytes);					\
      s->maxAllocated        = std::max(s->maxAllocated, s->currentlyAllocated); \
    }									\
  if (MemoryProfiler::debug)						\
    {									\
      std::cout << GridLogDebug << "[Memory debug] allocating " << memString(bytes) << std::endl; \
      profilerDebugPrint;						\
    }

#define profilerFree(bytes)						\
  if (MemoryProfiler::stats)						\
    {									\
      auto s = MemoryProfiler::stats;					\
      s->totalFreed         += (bytes);					\
      s->currentlyAllocated -= (bytes);					\
    }									\
  if (MemoryProfiler::debug)						\
    {									\
      std::cout << GridLogDebug << "[Memory debug] freeing " << memString(bytes) << std::endl; \
      profilerDebugPrint;						\
    }

void check_huge_pages(void *Buf,uint64_t BYTES);

NAMESPACE_END(Grid);

