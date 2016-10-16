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

// Code within namespace sha256 are from Stephan Brumme.
// see http://create.stephan-brumme.com/disclaimer.html

#ifndef RNG_STATE_RNG_STATE_H
#define RNG_STATE_RNG_STATE_H

#include "show.h"

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
  const uint64_t* p = (const uint64_t*)&rs.gaussian;
  splitTwoUint32(v[18], v[19], *p);
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
  uint64_t g = patchTwoUint32(v[18], v[19]);
  rs.gaussian = reinterpret_cast<double&>(g);
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

namespace sha256 {

  const size_t BlockSize = 512 / 8;

  const size_t HashBytes = 32;

  const size_t HashValues = HashBytes / 4;

  inline uint32_t rotate(uint32_t a, uint32_t c)
  {
    return (a >> c) | (a << (32 - c));
  }

  inline uint32_t swap(uint32_t x)
  {
    return (x >> 24) |
      ((x >>  8) & 0x0000FF00) |
      ((x <<  8) & 0x00FF0000) |
      (x << 24);
  }

  inline uint32_t f1(uint32_t e, uint32_t f, uint32_t g)
  // mix functions for processBlock()
  {
    uint32_t term1 = rotate(e, 6) ^ rotate(e, 11) ^ rotate(e, 25);
    uint32_t term2 = (e & f) ^ (~e & g); //(g ^ (e & (f ^ g)))
    return term1 + term2;
  }

  inline uint32_t f2(uint32_t a, uint32_t b, uint32_t c)
  // mix functions for processBlock()
  {
    uint32_t term1 = rotate(a, 2) ^ rotate(a, 13) ^ rotate(a, 22);
    uint32_t term2 = ((a | b) & c) | (a & b); //(a & (b ^ c)) ^ (b & c);
    return term1 + term2;
  }

  inline void processBlock(uint32_t newHash[8], const uint32_t oldHash[8], const uint8_t data[64])
    // process 64 bytes of data
    // newHash and oldHash and be the same
  {
    // get last hash
    uint32_t a = oldHash[0];
    uint32_t b = oldHash[1];
    uint32_t c = oldHash[2];
    uint32_t d = oldHash[3];
    uint32_t e = oldHash[4];
    uint32_t f = oldHash[5];
    uint32_t g = oldHash[6];
    uint32_t h = oldHash[7];
    // data represented as 16x 32-bit words
    const uint32_t* input = (uint32_t*) data;
    // convert to big endian
    uint32_t words[64];
    int i;
    for (i = 0; i < 16; i++) {
#if defined(__BYTE_ORDER) && (__BYTE_ORDER != 0) && (__BYTE_ORDER == __BIG_ENDIAN)
      words[i] =      input[i];
#else
      words[i] = swap(input[i]);
#endif
    }
    uint32_t x,y; // temporaries
    // first round
    x = h + f1(e,f,g) + 0x428a2f98 + words[ 0]; y = f2(a,b,c); d += x; h = x + y;
    x = g + f1(d,e,f) + 0x71374491 + words[ 1]; y = f2(h,a,b); c += x; g = x + y;
    x = f + f1(c,d,e) + 0xb5c0fbcf + words[ 2]; y = f2(g,h,a); b += x; f = x + y;
    x = e + f1(b,c,d) + 0xe9b5dba5 + words[ 3]; y = f2(f,g,h); a += x; e = x + y;
    x = d + f1(a,b,c) + 0x3956c25b + words[ 4]; y = f2(e,f,g); h += x; d = x + y;
    x = c + f1(h,a,b) + 0x59f111f1 + words[ 5]; y = f2(d,e,f); g += x; c = x + y;
    x = b + f1(g,h,a) + 0x923f82a4 + words[ 6]; y = f2(c,d,e); f += x; b = x + y;
    x = a + f1(f,g,h) + 0xab1c5ed5 + words[ 7]; y = f2(b,c,d); e += x; a = x + y;
    // secound round
    x = h + f1(e,f,g) + 0xd807aa98 + words[ 8]; y = f2(a,b,c); d += x; h = x + y;
    x = g + f1(d,e,f) + 0x12835b01 + words[ 9]; y = f2(h,a,b); c += x; g = x + y;
    x = f + f1(c,d,e) + 0x243185be + words[10]; y = f2(g,h,a); b += x; f = x + y;
    x = e + f1(b,c,d) + 0x550c7dc3 + words[11]; y = f2(f,g,h); a += x; e = x + y;
    x = d + f1(a,b,c) + 0x72be5d74 + words[12]; y = f2(e,f,g); h += x; d = x + y;
    x = c + f1(h,a,b) + 0x80deb1fe + words[13]; y = f2(d,e,f); g += x; c = x + y;
    x = b + f1(g,h,a) + 0x9bdc06a7 + words[14]; y = f2(c,d,e); f += x; b = x + y;
    x = a + f1(f,g,h) + 0xc19bf174 + words[15]; y = f2(b,c,d); e += x; a = x + y;
    // extend to 24 words
    for (; i < 24; i++)
      words[i] = words[i-16] +
        (rotate(words[i-15],  7) ^ rotate(words[i-15], 18) ^ (words[i-15] >>  3)) +
        words[i-7] +
        (rotate(words[i- 2], 17) ^ rotate(words[i- 2], 19) ^ (words[i- 2] >> 10));
    // third round
    x = h + f1(e,f,g) + 0xe49b69c1 + words[16]; y = f2(a,b,c); d += x; h = x + y;
    x = g + f1(d,e,f) + 0xefbe4786 + words[17]; y = f2(h,a,b); c += x; g = x + y;
    x = f + f1(c,d,e) + 0x0fc19dc6 + words[18]; y = f2(g,h,a); b += x; f = x + y;
    x = e + f1(b,c,d) + 0x240ca1cc + words[19]; y = f2(f,g,h); a += x; e = x + y;
    x = d + f1(a,b,c) + 0x2de92c6f + words[20]; y = f2(e,f,g); h += x; d = x + y;
    x = c + f1(h,a,b) + 0x4a7484aa + words[21]; y = f2(d,e,f); g += x; c = x + y;
    x = b + f1(g,h,a) + 0x5cb0a9dc + words[22]; y = f2(c,d,e); f += x; b = x + y;
    x = a + f1(f,g,h) + 0x76f988da + words[23]; y = f2(b,c,d); e += x; a = x + y;
    // extend to 32 words
    for (; i < 32; i++)
      words[i] = words[i-16] +
        (rotate(words[i-15],  7) ^ rotate(words[i-15], 18) ^ (words[i-15] >>  3)) +
        words[i-7] +
        (rotate(words[i- 2], 17) ^ rotate(words[i- 2], 19) ^ (words[i- 2] >> 10));
    // fourth round
    x = h + f1(e,f,g) + 0x983e5152 + words[24]; y = f2(a,b,c); d += x; h = x + y;
    x = g + f1(d,e,f) + 0xa831c66d + words[25]; y = f2(h,a,b); c += x; g = x + y;
    x = f + f1(c,d,e) + 0xb00327c8 + words[26]; y = f2(g,h,a); b += x; f = x + y;
    x = e + f1(b,c,d) + 0xbf597fc7 + words[27]; y = f2(f,g,h); a += x; e = x + y;
    x = d + f1(a,b,c) + 0xc6e00bf3 + words[28]; y = f2(e,f,g); h += x; d = x + y;
    x = c + f1(h,a,b) + 0xd5a79147 + words[29]; y = f2(d,e,f); g += x; c = x + y;
    x = b + f1(g,h,a) + 0x06ca6351 + words[30]; y = f2(c,d,e); f += x; b = x + y;
    x = a + f1(f,g,h) + 0x14292967 + words[31]; y = f2(b,c,d); e += x; a = x + y;
    // extend to 40 words
    for (; i < 40; i++)
      words[i] = words[i-16] +
        (rotate(words[i-15],  7) ^ rotate(words[i-15], 18) ^ (words[i-15] >>  3)) +
        words[i-7] +
        (rotate(words[i- 2], 17) ^ rotate(words[i- 2], 19) ^ (words[i- 2] >> 10));
    // fifth round
    x = h + f1(e,f,g) + 0x27b70a85 + words[32]; y = f2(a,b,c); d += x; h = x + y;
    x = g + f1(d,e,f) + 0x2e1b2138 + words[33]; y = f2(h,a,b); c += x; g = x + y;
    x = f + f1(c,d,e) + 0x4d2c6dfc + words[34]; y = f2(g,h,a); b += x; f = x + y;
    x = e + f1(b,c,d) + 0x53380d13 + words[35]; y = f2(f,g,h); a += x; e = x + y;
    x = d + f1(a,b,c) + 0x650a7354 + words[36]; y = f2(e,f,g); h += x; d = x + y;
    x = c + f1(h,a,b) + 0x766a0abb + words[37]; y = f2(d,e,f); g += x; c = x + y;
    x = b + f1(g,h,a) + 0x81c2c92e + words[38]; y = f2(c,d,e); f += x; b = x + y;
    x = a + f1(f,g,h) + 0x92722c85 + words[39]; y = f2(b,c,d); e += x; a = x + y;
    // extend to 48 words
    for (; i < 48; i++)
      words[i] = words[i-16] +
        (rotate(words[i-15],  7) ^ rotate(words[i-15], 18) ^ (words[i-15] >>  3)) +
        words[i-7] +
        (rotate(words[i- 2], 17) ^ rotate(words[i- 2], 19) ^ (words[i- 2] >> 10));
    // sixth round
    x = h + f1(e,f,g) + 0xa2bfe8a1 + words[40]; y = f2(a,b,c); d += x; h = x + y;
    x = g + f1(d,e,f) + 0xa81a664b + words[41]; y = f2(h,a,b); c += x; g = x + y;
    x = f + f1(c,d,e) + 0xc24b8b70 + words[42]; y = f2(g,h,a); b += x; f = x + y;
    x = e + f1(b,c,d) + 0xc76c51a3 + words[43]; y = f2(f,g,h); a += x; e = x + y;
    x = d + f1(a,b,c) + 0xd192e819 + words[44]; y = f2(e,f,g); h += x; d = x + y;
    x = c + f1(h,a,b) + 0xd6990624 + words[45]; y = f2(d,e,f); g += x; c = x + y;
    x = b + f1(g,h,a) + 0xf40e3585 + words[46]; y = f2(c,d,e); f += x; b = x + y;
    x = a + f1(f,g,h) + 0x106aa070 + words[47]; y = f2(b,c,d); e += x; a = x + y;
    // extend to 56 words
    for (; i < 56; i++)
      words[i] = words[i-16] +
        (rotate(words[i-15],  7) ^ rotate(words[i-15], 18) ^ (words[i-15] >>  3)) +
        words[i-7] +
        (rotate(words[i- 2], 17) ^ rotate(words[i- 2], 19) ^ (words[i- 2] >> 10));
    // seventh round
    x = h + f1(e,f,g) + 0x19a4c116 + words[48]; y = f2(a,b,c); d += x; h = x + y;
    x = g + f1(d,e,f) + 0x1e376c08 + words[49]; y = f2(h,a,b); c += x; g = x + y;
    x = f + f1(c,d,e) + 0x2748774c + words[50]; y = f2(g,h,a); b += x; f = x + y;
    x = e + f1(b,c,d) + 0x34b0bcb5 + words[51]; y = f2(f,g,h); a += x; e = x + y;
    x = d + f1(a,b,c) + 0x391c0cb3 + words[52]; y = f2(e,f,g); h += x; d = x + y;
    x = c + f1(h,a,b) + 0x4ed8aa4a + words[53]; y = f2(d,e,f); g += x; c = x + y;
    x = b + f1(g,h,a) + 0x5b9cca4f + words[54]; y = f2(c,d,e); f += x; b = x + y;
    x = a + f1(f,g,h) + 0x682e6ff3 + words[55]; y = f2(b,c,d); e += x; a = x + y;
    // extend to 64 words
    for (; i < 64; i++)
      words[i] = words[i-16] +
        (rotate(words[i-15],  7) ^ rotate(words[i-15], 18) ^ (words[i-15] >>  3)) +
        words[i-7] +
        (rotate(words[i- 2], 17) ^ rotate(words[i- 2], 19) ^ (words[i- 2] >> 10));
    // eigth round
    x = h + f1(e,f,g) + 0x748f82ee + words[56]; y = f2(a,b,c); d += x; h = x + y;
    x = g + f1(d,e,f) + 0x78a5636f + words[57]; y = f2(h,a,b); c += x; g = x + y;
    x = f + f1(c,d,e) + 0x84c87814 + words[58]; y = f2(g,h,a); b += x; f = x + y;
    x = e + f1(b,c,d) + 0x8cc70208 + words[59]; y = f2(f,g,h); a += x; e = x + y;
    x = d + f1(a,b,c) + 0x90befffa + words[60]; y = f2(e,f,g); h += x; d = x + y;
    x = c + f1(h,a,b) + 0xa4506ceb + words[61]; y = f2(d,e,f); g += x; c = x + y;
    x = b + f1(g,h,a) + 0xbef9a3f7 + words[62]; y = f2(c,d,e); f += x; b = x + y;
    x = a + f1(f,g,h) + 0xc67178f2 + words[63]; y = f2(b,c,d); e += x; a = x + y;
    // update hash
    newHash[0] = a + oldHash[0];
    newHash[1] = b + oldHash[1];
    newHash[2] = c + oldHash[2];
    newHash[3] = d + oldHash[3];
    newHash[4] = e + oldHash[4];
    newHash[5] = f + oldHash[5];
    newHash[6] = g + oldHash[6];
    newHash[7] = h + oldHash[7];
  }

  inline void processInput(
      uint32_t hash[8],
      const uint32_t oldHash[8], const uint64_t numBytes,
      const uint8_t* input, const size_t inputSize)
    // process final block, less than 64 bytes
    // newHash and oldHash and be the same
  {
    // the input bytes are considered as bits strings, where the first bit is the most significant bit of the byte
    // - append "1" bit to message
    // - append "0" bits until message length in bit mod 512 is 448
    // - append length as 64 bit integer
    // process initial parts of input
    std::memmove(hash, oldHash, 32);
    const int nBlocks = inputSize / 64;
    for (int i = 0; i < nBlocks; ++i) {
      processBlock(hash, hash, input + i * 64);
    }
    // initialize buffer from input
    const size_t bufferSize = inputSize - nBlocks * 64;
    unsigned char buffer[BlockSize];
    std::memcpy(buffer, input + nBlocks * 64, bufferSize);
    // number of bits
    size_t paddedLength = bufferSize * 8;
    // plus one bit set to 1 (always appended)
    paddedLength++;
    // number of bits must be (numBits % 512) = 448
    size_t lower11Bits = paddedLength & 511;
    if (lower11Bits <= 448) {
      paddedLength +=       448 - lower11Bits;
    } else {
      paddedLength += 512 + 448 - lower11Bits;
    }
    // convert from bits to bytes
    paddedLength /= 8;
    // only needed if additional data flows over into a second block
    unsigned char extra[BlockSize];
    // append a "1" bit, 128 => binary 10000000
    if (bufferSize < BlockSize) {
      buffer[bufferSize] = 128;
    } else {
      extra[0] = 128;
    }
    size_t i;
    for (i = bufferSize + 1; i < BlockSize; i++) {
      buffer[i] = 0;
    }
    for (; i < paddedLength; i++) {
      extra[i - BlockSize] = 0;
    }
    // add message length in bits as 64 bit number
    uint64_t msgBits = 8 * (numBytes + inputSize);
    // find right position
    unsigned char* addLength;
    if (paddedLength < BlockSize) {
      addLength = buffer + paddedLength;
    } else {
      addLength = extra + paddedLength - BlockSize;
    }
    // must be big endian
    *addLength++ = (unsigned char)((msgBits >> 56) & 0xFF);
    *addLength++ = (unsigned char)((msgBits >> 48) & 0xFF);
    *addLength++ = (unsigned char)((msgBits >> 40) & 0xFF);
    *addLength++ = (unsigned char)((msgBits >> 32) & 0xFF);
    *addLength++ = (unsigned char)((msgBits >> 24) & 0xFF);
    *addLength++ = (unsigned char)((msgBits >> 16) & 0xFF);
    *addLength++ = (unsigned char)((msgBits >>  8) & 0xFF);
    *addLength   = (unsigned char)( msgBits        & 0xFF);
    // process blocks
    processBlock(hash, hash, buffer);
    // flowed over into a second block ?
    if (paddedLength > BlockSize) {
      processBlock(hash, hash, extra);
    }
  }

}

inline void reset(RngState& rs)
{
  std::memset(&rs, 0, sizeof(RngState));
  rs.numBytes = 0;
  rs.hash[0] = 0x6a09e667;
  rs.hash[1] = 0xbb67ae85;
  rs.hash[2] = 0x3c6ef372;
  rs.hash[3] = 0xa54ff53a;
  rs.hash[4] = 0x510e527f;
  rs.hash[5] = 0x9b05688c;
  rs.hash[6] = 0x1f83d9ab;
  rs.hash[7] = 0x5be0cd19;
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

inline void splitRngState(RngState& rs, const RngState& rs0, const std::string& sindex)
  // produce a new rng ``rs'' uniquely identified by ``rs0'' and ``sindex''
  // will not affect old rng ``rs0''
  // the function should behave correctly even if ``rs'' is actually ``rs0''
{
  std::string data = ssprintf("[%lu] {%s}", rs0.index, sindex.c_str());
  const int nBlocks = (data.length() - 1) / 64 + 1;
  data.resize(nBlocks * 64, ' ');
  sha256::processBlock(rs.hash, rs0.hash, (const uint8_t*)data.c_str());
  for (int i = 1; i < nBlocks; ++i) {
    sha256::processBlock(rs.hash, rs.hash, (const uint8_t*)data.c_str() + i * 64);
  }
  rs.numBytes = rs0.numBytes + nBlocks * 64;
  rs.index = 0;
  rs.cache[0] = 0;
  rs.cache[1] = 0;
  rs.cache[2] = 0;
  rs.gaussian = 0.0;
  rs.cacheAvail = 0;
  rs.gaussianAvail = false;
}

inline void computeHashWithInput(uint32_t hash[8], const RngState& rs, const std::string& input)
{
  sha256::processInput(hash, rs.hash, rs.numBytes, (const uint8_t*)input.c_str(), input.length());
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
