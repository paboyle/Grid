// vim: set ts=2 sw=2 expandtab:

// Copyright (c) 2016 Luchang Jin
// All rights reserved.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#ifndef INCLUDE_RNG_STATE_H
#define INCLUDE_RNG_STATE_H

#include "show.h"

#ifndef USE_OPENSSL
#include "sha256.h"
#else
#include <openssl/sha.h>
#endif

#include <stdint.h>
#include <endian.h>
#include <cstring>
#include <cmath>
#include <cassert>
#include <string>
#include <ostream>
#include <istream>
#include <vector>

#ifdef CURRENT_DEFAULT_NAMESPACE_NAME
namespace CURRENT_DEFAULT_NAMESPACE_NAME {
#endif

struct RngState;

inline void reset(RngState& rs);

inline void reset(RngState& rs, const std::string& seed);

inline void reset(RngState& rs, const long seed)
{
  reset(rs, show(seed));
}

inline void splitRngState(RngState& rs, const RngState& rs0, const std::string& sindex);

inline void splitRngState(RngState& rs, const RngState& rs0, const long sindex = 0)
{
  splitRngState(rs, rs0, show(sindex));
}

inline uint64_t randGen(RngState& rs);

inline double uRandGen(RngState& rs, const double upper = 1.0, const double lower = 0.0);

inline double gRandGen(RngState& rs, const double sigma = 1.0, const double center = 0.0);

inline void computeHashWithInput(uint32_t hash[8], const RngState& rs, const std::string& input);

struct RngState
{
  uint64_t numBytes;
  uint32_t hash[8];
  unsigned long index;
  //
  uint64_t cache[3];
  double gaussian;
  int cacheAvail;
  bool gaussianAvail;
  //
  inline void init()
  {
    reset(*this);
  }
  //
  RngState()
  {
    init();
  }
  RngState(const std::string& seed)
  {
    reset(*this, seed);
  }
  RngState(const long seed)
  {
    reset(*this, seed);
  }
  RngState(const RngState& rs0, const std::string& sindex)
  {
    splitRngState(*this, rs0, sindex);
  }
  RngState(const RngState& rs0, const long sindex)
  {
    splitRngState(*this, rs0, sindex);
  }
  //
  RngState split(const std::string& sindex)
  {
    RngState rs(*this, sindex);
    return rs;
  }
  RngState split(const long sindex)
  {
    RngState rs(*this, sindex);
    return rs;
  }
};

const size_t RNG_STATE_NUM_OF_INT32 = 2 + 8 + 2 + 3 * 2 + 2 + 1 + 1;

inline uint64_t patchTwoUint32(const uint32_t a, const uint32_t b)
{
  return (uint64_t)a << 32 | (uint64_t)b;
}

inline void splitTwoUint32(uint32_t& a, uint32_t& b, const uint64_t x)
{
  b = (uint32_t)x;
  a = (uint32_t)(x >> 32);
  assert(x == patchTwoUint32(a, b));
}

inline void exportRngState(uint32_t* v, const RngState& rs)
{
  assert(22 == RNG_STATE_NUM_OF_INT32);
  splitTwoUint32(v[0], v[1], rs.numBytes);
  for (int i = 0; i < 8; ++i) {
    v[2 + i] = rs.hash[i];
  }
  splitTwoUint32(v[10], v[11], rs.index);
  for (int i = 0; i < 3; ++i) {
    splitTwoUint32(v[12 + i * 2], v[12 + i * 2 + 1], rs.cache[i]);
  }
  union {
    double d;
    uint64_t l;
  } g;
  g.d = rs.gaussian;
  splitTwoUint32(v[18], v[19], g.l);
  v[20] = rs.cacheAvail;
  v[21] = rs.gaussianAvail;
}

inline void importRngState(RngState& rs, const uint32_t* v)
{
  assert(22 == RNG_STATE_NUM_OF_INT32);
  rs.numBytes = patchTwoUint32(v[0], v[1]);
  for (int i = 0; i < 8; ++i) {
    rs.hash[i] = v[2 + i];
  }
  rs.index = patchTwoUint32(v[10], v[11]);
  for (int i = 0; i < 3; ++i) {
    rs.cache[i] = patchTwoUint32(v[12 + i * 2], v[12 + i * 2 + 1]);
  }
  union {
    double d;
    uint64_t l;
  } g;
  g.l = patchTwoUint32(v[18], v[19]);
  rs.gaussian = g.d;
  rs.cacheAvail = v[20];
  rs.gaussianAvail = v[21];
}

inline void exportRngState(std::vector<uint32_t>& v, const RngState& rs)
{
  v.resize(RNG_STATE_NUM_OF_INT32);
  exportRngState(v.data(), rs);
}

inline void importRngState(RngState& rs, const std::vector<uint32_t>& v)
{
  assert(RNG_STATE_NUM_OF_INT32 == v.size());
  importRngState(rs, v.data());
}

inline std::ostream& operator<<(std::ostream& os, const RngState& rs)
{
  std::vector<uint32_t> v(RNG_STATE_NUM_OF_INT32);
  exportRngState(v, rs);
  for (size_t i = 0; i < v.size() - 1; ++i) {
    os << v[i] << " ";
  }
  os << v.back();
  return os;
}

inline std::istream& operator>>(std::istream& is, RngState& rs)
{
  std::vector<uint32_t> v(RNG_STATE_NUM_OF_INT32);
  for (size_t i = 0; i < v.size(); ++i) {
    is >> v[i];
  }
  importRngState(rs, v);
  return is;
}

inline std::string show(const RngState& rs)
{
  return shows(rs);
}

inline bool operator==(const RngState& rs1, const RngState& rs2)
{
  return 0 == memcmp(&rs1, &rs2, sizeof(RngState));
}

inline void reset(RngState& rs)
{
  std::memset(&rs, 0, sizeof(RngState));
  rs.numBytes = 0;
  rs.hash[0] = 0;
  rs.hash[1] = 0;
  rs.hash[2] = 0;
  rs.hash[3] = 0;
  rs.hash[4] = 0;
  rs.hash[5] = 0;
  rs.hash[6] = 0;
  rs.hash[7] = 0;
  rs.index = 0;
  rs.cache[0] = 0;
  rs.cache[1] = 0;
  rs.cache[2] = 0;
  rs.gaussian = 0.0;
  rs.cacheAvail = 0;
  rs.gaussianAvail = false;
}

inline void reset(RngState& rs, const std::string& seed)
{
  reset(rs);
  splitRngState(rs, rs, seed);
}

inline void computeHashWithInput(uint32_t hash[8], const RngState& rs, const std::string& input)
{
  std::string data(32, ' ');
  for (int i = 0; i < 8; ++i) {
    data[i*4 + 0] = (rs.hash[i] >> 24) & 0xFF;
    data[i*4 + 1] = (rs.hash[i] >> 16) & 0xFF;
    data[i*4 + 2] = (rs.hash[i] >>  8) & 0xFF;
    data[i*4 + 3] =  rs.hash[i]        & 0xFF;
  }
  data += input;
#ifndef USE_OPENSSL
  sha256::computeHash(hash, (const uint8_t*)data.c_str(), data.length());
#else
  {
    uint8_t rawHash[32];
    SHA256((unsigned char*)data.c_str(), data.length(), rawHash);
    for (int i = 0; i < 8; ++i) {
      hash[i] = (((uint32_t)rawHash[i*4 + 0]) << 24)
              + (((uint32_t)rawHash[i*4 + 1]) << 16)
              + (((uint32_t)rawHash[i*4 + 2]) <<  8)
              + ( (uint32_t)rawHash[i*4 + 3]);
    }
  }
#endif
}

inline void splitRngState(RngState& rs, const RngState& rs0, const std::string& sindex)
  // produce a new rng ``rs'' uniquely identified by ``rs0'' and ``sindex''
  // will not affect old rng ``rs0''
  // the function should behave correctly even if ``rs'' is actually ``rs0''
{
  std::string input = ssprintf("[%lu] {%s}", rs0.index, sindex.c_str());
  rs.numBytes = rs0.numBytes + 64 * ((32 + input.length() + 1 + 8 - 1) / 64 + 1);
  computeHashWithInput(rs.hash, rs0, input);
  rs.index = 0;
  rs.cache[0] = 0;
  rs.cache[1] = 0;
  rs.cache[2] = 0;
  rs.gaussian = 0.0;
  rs.cacheAvail = 0;
  rs.gaussianAvail = false;
}

inline uint64_t randGen(RngState& rs)
{
  assert(0 <= rs.cacheAvail && rs.cacheAvail <= 3);
  rs.index += 1;
  if (rs.cacheAvail > 0) {
    rs.cacheAvail -= 1;
    uint64_t r = rs.cache[rs.cacheAvail];
    rs.cache[rs.cacheAvail] = 0;
    return r;
  } else {
    uint32_t hash[8];
    computeHashWithInput(hash, rs, ssprintf("[%lu]", rs.index));
    rs.cache[0] = patchTwoUint32(hash[0], hash[1]);
    rs.cache[1] = patchTwoUint32(hash[2], hash[3]);
    rs.cache[2] = patchTwoUint32(hash[4], hash[5]);
    rs.cacheAvail = 3;
    return patchTwoUint32(hash[6], hash[7]);
  }
}

inline double uRandGen(RngState& rs, const double upper, const double lower)
{
  uint64_t u = randGen(rs);
  const double fac = 1.0 / (256.0 * 256.0 * 256.0 * 256.0) / (256.0 * 256.0 * 256.0 * 256.0);
  return u * fac * (upper - lower) + lower;
}

inline double gRandGen(RngState& rs, const double sigma, const double center)
{
  rs.index += 1;
  if (rs.gaussianAvail) {
    rs.gaussianAvail = false;
    return rs.gaussian * sigma + center;
  } else {
    // pick 2 uniform numbers in the square extending from
    // -1 to 1 in each direction, see if they are in the
    // unit circle, and try again if they are not.
    int num_try = 1;
    double v1, v2, rsq;
    do {
      v1 = uRandGen(rs, 1.0, -1.0);
      v2 = uRandGen(rs, 1.0, -1.0);
      if ((num_try % 1000)==0) {
        printf("gRandGen : WARNING num_try=%d v1=%e v2=%e\n",num_try,v1,v2);
      }
      rsq = v1*v1 + v2*v2;
      num_try++;
    } while ((num_try < 10000) && (rsq >= 1.0 || rsq == 0));
    if (num_try > 9999) {
      printf("gRandGen : WARNING failed after 10000 tries (corrupted RNG?), returning ridiculous numbers (1e+10)\n");
      return 1e+10;
    }
    double fac = std::sqrt(-2.0 * std::log(rsq)/rsq);
    rs.gaussian = v1 * fac;
    rs.gaussianAvail = true;
    return v2 * fac * sigma + center;
  }
}

#ifdef CURRENT_DEFAULT_NAMESPACE_NAME
}
#endif

#endif
