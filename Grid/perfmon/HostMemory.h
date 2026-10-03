/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/perfmon/HostMemory.h

    Copyright (C) 2026

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

#include <sys/resource.h>
#include <unistd.h>
#include <cctype>
#include <fstream>
#ifdef __APPLE__
#include <mach/mach.h>
#endif

NAMESPACE_BEGIN(Grid);

/////////////////////////////////////////////////////////////////////////////////////////////
// Host memory as the OOM killer sees it, plus the device allocator footprint.
// Targets are Apple or Linux; on Apple only the RSS and device fields are filled.
//
//  RSS      : this process, now and its high-water mark (getrusage).
//  Node     : MemAvailable and MemTotal from /proc/meminfo (Linux).
//  Cgroup   : the nearest enclosing memory cgroup with a finite limit (the one the batch
//             system enforces), its usage now, its high-water mark and the limit (Linux,
//             v2 or v1 hierarchy). Cgroup usage includes page cache, tmpfs and kernel memory
//             charged to the job, none of which appear in any process RSS.
//  Device   : MemoryManager Lattice footprint on the device and its free-block cache. Device
//             allocations outside the MemoryManager (comms buffers) are not included.
//
// All values in GB; -1 where the source does not exist on this platform or kernel.
/////////////////////////////////////////////////////////////////////////////////////////////
struct HostMemoryStatus
{
  RealD RSS           = -1;
  RealD RSSPeak       = -1;
  RealD NodeAvailable = -1;
  RealD NodeTotal     = -1;
  RealD CgroupCurrent = -1;
  RealD CgroupPeak    = -1;
  RealD CgroupLimit   = -1;
  RealD DeviceLattice = -1;
  RealD DeviceCache   = -1;
};

#ifndef __APPLE__ // Linux

// One number from a file; -1 if the file is absent or holds no number (e.g. "max")
inline RealD HostMemoryReadBytes(const std::string &file)
{
  std::ifstream f(file);
  if ( !f.good() ) {
    return -1;
  }
  std::string word;
  f >> word;
  if ( word.empty() || !isdigit(word[0]) ) {
    return -1;
  }
  return std::stod(word);
}

// Value of a "Key: N kB" line in /proc/meminfo, in bytes; -1 if absent
inline RealD HostMemoryMeminfo(const std::string &key)
{
  std::ifstream f("/proc/meminfo");
  std::string name;
  RealD value;
  std::string unit;
  while ( f >> name >> value >> unit ) {
    if ( name == key + ":" ) {
      return value*1024.0;
    }
  }
  return -1;
}

// Fill the cgroup fields from the nearest enclosing memory cgroup with a finite limit, or
// from this process's own cgroup (limit left at -1) if no ancestor has one
inline void HostMemoryCgroup(HostMemoryStatus &s)
{
  std::ifstream f("/proc/self/cgroup");
  std::string line;
  std::string root;
  std::string path;
  std::string current;
  std::string peak;
  std::string limit;
  while ( std::getline(f,line) ) {
    size_t c1 = line.find(':');
    size_t c2 = line.find(':',c1+1);
    if ( c1 == std::string::npos || c2 == std::string::npos ) {
      continue;
    }
    std::string controllers = line.substr(c1+1,c2-c1-1);
    // v2 unified hierarchy: "0::/path"
    if ( controllers.empty() && root.empty() ) {
      root    = "/sys/fs/cgroup";
      path    = line.substr(c2+1);
      current = "memory.current";
      peak    = "memory.peak";
      limit   = "memory.max";
    }
    // v1 memory controller: "N:memory:/path", takes precedence
    if ( (","+controllers+",").find(",memory,") != std::string::npos ) {
      root    = "/sys/fs/cgroup/memory";
      path    = line.substr(c2+1);
      current = "memory.usage_in_bytes";
      peak    = "memory.max_usage_in_bytes";
      limit   = "memory.limit_in_bytes";
      break;
    }
  }
  if ( root.empty() ) {
    return;
  }

  // v1 reports "no limit" as a number near 2^63
  const RealD unlimited = 1.0e18;
  std::string dir = path;
  while ( true ) {
    std::string base = root + dir + "/";
    RealD lim = HostMemoryReadBytes(base+limit);
    if ( lim > 0 && lim < unlimited ) {
      s.CgroupCurrent = HostMemoryReadBytes(base+current)/1.0e9;
      s.CgroupPeak    = HostMemoryReadBytes(base+peak)/1.0e9;
      s.CgroupLimit   = lim/1.0e9;
      return;
    }
    if ( dir.empty() || dir == "/" ) {
      break;
    }
    dir = dir.substr(0,dir.rfind('/'));
  }
  std::string base = root + path + "/";
  s.CgroupCurrent = HostMemoryReadBytes(base+current)/1.0e9;
  s.CgroupPeak    = HostMemoryReadBytes(base+peak)/1.0e9;
}

#endif

// This rank only; no communication
inline HostMemoryStatus HostMemoryQuery(void)
{
  HostMemoryStatus s;

  struct rusage ru;
  getrusage(RUSAGE_SELF,&ru);
#ifdef __APPLE__
  s.RSSPeak = ru.ru_maxrss/1.0e9;                // bytes on macOS
  mach_task_basic_info_data_t info;
  mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
  if ( task_info(mach_task_self(),MACH_TASK_BASIC_INFO,(task_info_t)&info,&count) == KERN_SUCCESS ) {
    s.RSS = info.resident_size/1.0e9;
  }
#else // Linux
  s.RSSPeak = ru.ru_maxrss*1024.0/1.0e9;         // kilobytes on Linux
  std::ifstream statm("/proc/self/statm");
  long pages    = 0;
  long resident = 0;
  if ( statm >> pages >> resident ) {
    s.RSS = resident*(RealD)sysconf(_SC_PAGESIZE)/1.0e9;
  }
  RealD avail = HostMemoryMeminfo("MemAvailable");
  RealD total = HostMemoryMeminfo("MemTotal");
  if ( avail >= 0 ) {
    s.NodeAvailable = avail/1.0e9;
  }
  if ( total >= 0 ) {
    s.NodeTotal = total/1.0e9;
  }
  HostMemoryCgroup(s);
#endif

  s.DeviceLattice = MemoryManager::DeviceBytes/1.0e9;
  s.DeviceCache   = MemoryManager::DeviceCacheBytes()/1.0e9;
  return s;
}

/////////////////////////////////////////////////////////////////////////////////////////////
// Collective over grid. One line on log: maxima over ranks, except node MemAvailable which
// is the minimum over ranks (the tightest node).
/////////////////////////////////////////////////////////////////////////////////////////////
inline void HostMemoryReport(GridBase *grid,GridLogger &log,const std::string &label)
{
  HostMemoryStatus s = HostMemoryQuery();
  RealD minus_available = -s.NodeAvailable;
  grid->GlobalMax(s.RSS);
  grid->GlobalMax(s.RSSPeak);
  grid->GlobalMax(minus_available);
  grid->GlobalMax(s.NodeTotal);
  grid->GlobalMax(s.CgroupCurrent);
  grid->GlobalMax(s.CgroupPeak);
  grid->GlobalMax(s.CgroupLimit);
  grid->GlobalMax(s.DeviceLattice);
  grid->GlobalMax(s.DeviceCache);
  std::cout << log << "HOSTMEM " << label
            << " : rank RSS " << s.RSS << " peak " << s.RSSPeak
            << " ; node available " << -minus_available << " of " << s.NodeTotal
            << " ; cgroup " << s.CgroupCurrent << " peak " << s.CgroupPeak << " limit " << s.CgroupLimit
            << " ; device lattice " << s.DeviceLattice << " cache " << s.DeviceCache
            << " GB" << std::endl;
}

NAMESPACE_END(Grid);
