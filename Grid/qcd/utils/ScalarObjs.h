/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/utils/WilsonLoops.h

    Copyright (C) 2015

Author: neo <cossu@post.kek.jp>

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
#ifndef SCALAR_OBJS_H
#define SCALAR_OBJS_H
namespace Grid {

  // FIXME drop the QCD namespace in Nd
  

// Scalar field obs
template <class Impl>
class ScalarObs {
 public:
  //////////////////////////////////////////////////
  // squared field
  //////////////////////////////////////////////////

  static void phisquared(typename Impl::Field &fsq,
                         const typename Impl::Field &f) {
    fsq = f * f;
  }
  //////////////////////////////////////////////////
  // phi^4 interaction term
  //////////////////////////////////////////////////

  static void phifourth(typename Impl::Field &fsq,
                        const typename Impl::Field &f) {
    fsq = f * f * f * f;
  }

  //////////////////////////////////////////////////
  // phi(x)phi(x+mu)
  //////////////////////////////////////////////////

  static void phider(typename Impl::Field &fsq,
                     const typename Impl::Field &f) {
    fsq = Cshift(f, 0, -1) * f;
    for (int mu = 1; mu < QCD::Nd; mu++) fsq += Cshift(f, mu, -1) * f;
  }

  //////////////////////////////////////////////////
  // Vol sum of the previous obs.
  //////////////////////////////////////////////////

  static RealD sumphider(const typename Impl::Field &f) {
    typename Impl::Field tmp(f._grid);
    tmp = Cshift(f, 0, -1) * f;
    for (int mu = 1; mu < QCD::Nd; mu++) {
      tmp += Cshift(f, mu, -1) * f;
    }
    return -sum(trace(tmp));
  }

  static RealD sumphisquared(const typename Impl::Field &f) {
    typename Impl::Field tmp(f._grid);
    tmp = f * f;
    return sum(trace(tmp));
  }

  static RealD sumphifourth(const typename Impl::Field &f) {
    typename Impl::Field tmp(f._grid);
    phifourth(tmp, f);
    return sum(trace(tmp));
  }
};


}

#endif
