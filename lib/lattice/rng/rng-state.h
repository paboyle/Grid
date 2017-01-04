// vim: set ts=2 sw=2 expandtab:

// Copyright (c) 2017 Luchang Jin
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

// Code within namespace sha256 are originally from Stephan Brumme.
// see http://create.stephan-brumme.com/disclaimer.html

#pragma once

#ifndef INCLUDE_RNG_STATE_H
#define INCLUDE_RNG_STATE_H

#ifndef USE_OPENSSL
#include "sha256.h"
#else
#include <openssl/sha.h>
#endif

#include <stdint.h>
#include <cassert>
#include <string>
#include <cstring>
#include <cmath>
#include <cstdio>
#include <sstream>

#ifdef CURRENT_DEFAULT_NAMESPACE_NAME
namespace CURRENT_DEFAULT_NAMESPACE_NAME {
#endif

const int RngStateSizeUint32 = 6 * 2;

struct RngState
{
  uint32_t state[RngStateSizeUint32];
  //
  uint64_t cache[3];
  double gaussian;
  int cacheAvail;
  bool gaussianAvail;
  //
  RngState()
  {
    init();
  }
  //
  RngState(const uint64_t seed, const uint64_t shift = 0)
  {
    init(seed, shift);
  }
  //
  void init()
  {
    std::memset(state, 0, sizeof(state));
    cache[0] = 0;
    cache[1] = 0;
    cache[2] = 0;
    gaussian = 0.0;
    cacheAvail = 0;
    gaussianAvail= false;
  }
  //
  void init(const uint64_t seed, const uint64_t shift = 0)
  {
    init(seed, 0, 0, 0, shift);
  }
  //
  void init(const uint64_t seed, const uint64_t type, const uint64_t traj, const uint64_t index = 0, const uint64_t shift = 0);
};

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

inline void RngState::init(const uint64_t seed, const uint64_t type, const uint64_t traj, const uint64_t index, const uint64_t shift)
{
  init();
  uint64_t ss[RngStateSizeUint32/2];
  std::memset(ss, 0, sizeof(ss));
  ss[0] = seed;
  ss[1] = type;
  ss[2] = traj;
  ss[3] = index;
  ss[RngStateSizeUint32/2-1] = shift;
  // std::memcpy(state, ss, sizeof(ss));
  for (int i = 0; i < RngStateSizeUint32/2; ++i) {
    splitTwoUint32(state[i*2], state[i*2+1], ss[i]);
  }
}

const size_t RNG_STATE_NUM_OF_INT32 = RngStateSizeUint32 + 2 * 3 + 2 + 1 + 1;

inline void exportRngState(uint32_t* v, const RngState& rs)
{
  assert(22 == RNG_STATE_NUM_OF_INT32);
  assert(12 == RngStateSizeUint32);
  for (int i = 0; i < 12; ++i) {
    v[0 + i] = rs.state[i];
  }
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
  assert(12 == RngStateSizeUint32);
  for (int i = 0; i < 12; ++i) {
    rs.state[i] = v[0 + i];
  }
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
  std::ostringstream out;
  out << rs;
  return out.str();
}

inline bool operator==(const RngState& rs1, const RngState& rs2)
{
  return 0 == std::memcmp(&rs1, &rs2, sizeof(RngState));
}

inline bool operator!=(const RngState& rs1, const RngState& rs2)
{
  return !(rs1 == rs2);
}

inline void shiftRngState(RngState& rs, const uint64_t shift = 1)
{
  uint64_t pre = shift;
  int i = RngStateSizeUint32 - 1;
  do {
    pre += rs.state[i];
    rs.state[i] = (uint32_t)pre;
    pre = pre >> 32;
    i -= 1;
  } while (pre > 0 && i >= 0);
}

inline void computeHash(uint32_t hash[8], const RngState& rs)
{
  std::string data(RngStateSizeUint32 * sizeof(uint32_t), ' ');
  for (int i = 0; i < RngStateSizeUint32; ++i) {
    data[i*4 + 0] = (rs.state[i] >> 24) & 0xFF;
    data[i*4 + 1] = (rs.state[i] >> 16) & 0xFF;
    data[i*4 + 2] = (rs.state[i] >>  8) & 0xFF;
    data[i*4 + 3] =  rs.state[i]        & 0xFF;
  }
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

inline uint64_t randGen(RngState& rs)
{
  assert(0 <= rs.cacheAvail && rs.cacheAvail <= 3);
  if (rs.cacheAvail > 0) {
    rs.cacheAvail -= 1;
    uint64_t r = rs.cache[rs.cacheAvail];
    rs.cache[rs.cacheAvail] = 0;
    return r;
  } else {
    uint32_t hash[8];
    computeHash(hash, rs);
    shiftRngState(rs);
    rs.cache[0] = patchTwoUint32(hash[0], hash[1]);
    rs.cache[1] = patchTwoUint32(hash[2], hash[3]);
    rs.cache[2] = patchTwoUint32(hash[4], hash[5]);
    rs.cacheAvail = 3;
    return patchTwoUint32(hash[6], hash[7]);
  }
}

inline double uRandGen(RngState& rs, const double upper = 1.0, const double lower = 0.0)
{
  uint64_t u = randGen(rs);
  const double fac = 1.0 / (256.0 * 256.0 * 256.0 * 256.0) / (256.0 * 256.0 * 256.0 * 256.0);
  return u * fac * (upper - lower) + lower;
}

inline double gRandGen(RngState& rs, const double center = 0.0, const double sigma = 1.0)
{
  if (rs.gaussianAvail) {
    rs.gaussianAvail = false;
    return rs.gaussian * sigma + center;
  } else {
    const char* cname = "";
    const char* fname = "gRandGen()";
    int num_try = 1;
    double v1, v2, rsq;
    do {
      v1 = uRandGen(rs, 1.0, -1.0);
      v2 = uRandGen(rs, 1.0, -1.0);
      if ((num_try % 1000)==0) {
        std::printf("num_try=%d v1=%e v2=%e\n",num_try,v1,v2);
      }
      rsq = v1*v1 + v2*v2;
      num_try++;
    } while ((num_try < 10000) && (rsq >= 1.0 || rsq == 0));
    if (num_try > 9999) {
      std::printf("failed after 10000 tries (corrupted RNG?), returning ridiculous numbers (1e+10)\n");
      return 1e+10;
    }
    // pick 2 uniform numbers in the square extending from
    // -1 to 1 in each direction, see if they are in the
    // unit circle, and try again if they are not.
    //
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
