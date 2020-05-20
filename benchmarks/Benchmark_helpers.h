/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./benchmarks/Benchmark_helpers.h

    Copyright (C) 2015 - 2020

Author: Daniel Richtmann <daniel.richtmann@ur.de>

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

struct PerfNumbers {
  double seconds_, calls_, intensity_, perf_, traffic_;
  PerfNumbers(double seconds, double calls, double intensity, double perf, double traffic)
    : seconds_(seconds), calls_(calls), intensity_(intensity), perf_(perf), traffic_(traffic) {}
};

inline std::ostream& operator<<(std::ostream& stream, const PerfNumbers& pn) {
  // clang-format off
  return stream << std::scientific
                <<        pn.seconds_   << " s"
                << " " << pn.calls_     << " x"
                << " " << pn.intensity_ << " F/B"
                << " " << pn.perf_      << " F/s"
                << " " << pn.traffic_   << " B/s";
  // clang-format on
}

#define BenchmarkFunction(function, flop, byte, nIterMin, nSecMin, ...) \
  do { \
    uint64_t nIterFinal = 1; \
    if(nIterMin != 1) { \
      uint64_t nIterInitial = 10; \
      double t0 = usecond(); \
      for(uint64_t i = 0; i < nIterInitial; ++i) { \
        __SSC_START; \
        function(__VA_ARGS__); \
        __SSC_STOP; \
      } \
      double td  = (usecond()-t0)/1e6; \
      nIterFinal = std::min(nIterMin, (uint64_t)(nSecMin*nIterInitial/td)+1); \
    } \
    double t0 = usecond(); \
    for(uint64_t i = 0; i < nIterFinal; ++i) { \
      __SSC_START; \
      function(__VA_ARGS__); \
      __SSC_STOP; \
    } \
    double td = (usecond()-t0)/1e6; \
    PerfNumbers pn(td, nIterFinal, (byte == 0. ? 0. : flop/byte), flop*nIterFinal/td, byte*nIterFinal/td); \
    std::cout << GridLogPerformance << "Kernel : " #function << " : " << pn << std::endl; \
  } while(0)

#define BenchmarkFunctionMRHS(function, flop, byte, nIterMin, nSecMin, nRHS, ...) \
  do { \
    uint64_t nIterFinal = 1; \
    if(nIterMin != 1) { \
      uint64_t nIterInitial = 10; \
      double t0 = usecond(); \
      for(uint64_t i = 0; i < nIterInitial; ++i) { \
        for(uint64_t rhs = 0; rhs < nRHS; ++rhs) { \
          __SSC_START; \
          function(__VA_ARGS__); \
          __SSC_STOP; \
        } \
      } \
      double td  = (usecond()-t0)/1e6; \
      nIterFinal = std::min(nIterMin, (uint64_t)(nSecMin*nIterInitial/td)+1); \
    } \
    double t0 = usecond(); \
    for(uint64_t i = 0; i < nIterFinal; ++i) { \
      for(uint64_t rhs = 0; rhs < nRHS; ++rhs) { \
        __SSC_START; \
        function(__VA_ARGS__); \
        __SSC_STOP; \
      } \
    } \
    double td = (usecond()-t0)/1e6; \
    PerfNumbers pn(td, nIterFinal, (byte == 0. ? 0. : flop/byte), flop*nIterFinal/td, byte*nIterFinal/td); \
    std::cout << GridLogPerformance << "Kernel : " #function << " : " << pn << std::endl; \
  } while(0)

template<class Field>
void assertResultMatchesReference(RealD tolerance, const Field &reference, const Field &result) {
  conformable(reference.Grid(), result.Grid());
  Field diff(reference.Grid());

  diff        = reference - result;
  auto absDev = norm2(diff);
  auto relDev = absDev / norm2(reference);

  std::cout << GridLogMessage
            << "ref = " << norm2(reference)
            << " res = " << norm2(result)
            << " absolute deviation = " << absDev
            << " relative deviation = " << relDev;

  if(relDev <= tolerance) {
    std::cout << " < " << tolerance << " -> check passed" << std::endl;
  } else {
    std::cout << " > " << tolerance << " -> check failed" << std::endl;
  }
  assert(relDev <= tolerance);
}

NAMESPACE_END(Grid);
