// vim: set ts=2 sw=2 expandtab:

// Copyright (c) 2016 Luchang Jin
// All rights reserved.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
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

#include "rng-state.h"

#include <array>
#include <cstring>
#include <ostream>
#include <istream>

#ifdef CURRENT_DEFAULT_NAMESPACE_NAME
namespace CURRENT_DEFAULT_NAMESPACE_NAME {
#endif

struct SprngSha256
{
  RngState rs;
  //
  using result_type = uint64_t;
  //
  static constexpr result_type default_seed = 0;
  //
  explicit SprngSha256(result_type val = default_seed)
  {
    seed(val);
  }
  template<typename Sseq, typename = typename
    std::enable_if<!std::is_same<Sseq, SprngSha256>::value>
    ::type>
  explicit SprngSha256(Sseq& q)
  {
    seed(q);
  }
  //
  static constexpr result_type min()
  {
    return 0;
  }
  //
  static constexpr result_type max()
  {
    return UINT64_MAX;
  }
  //
  void seed(result_type val = default_seed)
  {
    reset(rs, (long)val);
  }
  template <class Sseq>
  typename std::enable_if<std::is_class<Sseq>::value>::type
  seed(Sseq& q)
  {
    std::array<uint32_t, 8> seq;
    q.generate(seq.begin(), seq.end());
    reset(rs);
    for (size_t i = 0; i < seq.size(); ++i) {
      splitRngState(rs, rs, seq[i]);
    }
  }
  //
  result_type operator()()
  {
    return randGen(rs);
  }
  //
  void discard(unsigned long long z)
  {
    for (unsigned long long i = 0; i < z; ++i) {
      randGen(rs);
    }
  }
};

inline std::ostream& operator<<(std::ostream& os, const SprngSha256& ss)
{
  os << ss.rs;
  return os;
}

inline std::istream& operator>>(std::istream& is, SprngSha256& ss)
{
  is >> ss.rs;
  return is;
}

inline bool operator==(const SprngSha256& ss1, const SprngSha256& ss2)
{
  return ss1.rs == ss2.rs;
}

#ifdef CURRENT_DEFAULT_NAMESPACE_NAME
}
#endif
