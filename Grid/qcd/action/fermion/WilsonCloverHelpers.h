/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/action/fermion/WilsonCloverHelpers.h

    Copyright (C) 2021 - 2022

    Author: Daniel Richtmann <daniel.richtmann@gmail.com>

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

// Helper routines that implement common clover functionality

NAMESPACE_BEGIN(Grid);

template<class Impl> class WilsonCloverHelpers {
public:

  INHERIT_IMPL_TYPES(Impl);
  INHERIT_CLOVER_TYPES(Impl);

  // Computing C_{\mu \nu}(x) as in Eq.(B.39) in Zbigniew Sroczynski's PhD thesis
  static GaugeLinkField Cmunu(std::vector<GaugeLinkField> &U, GaugeLinkField &lambda, int mu, int nu)
  {
    conformable(lambda.Grid(), U[0].Grid());
    GaugeLinkField out(lambda.Grid()), tmp(lambda.Grid());
    // insertion in upper staple
    // please check redundancy of shift operations

    // C1+
    tmp = lambda * U[nu];
    out = Impl::ShiftStaple(Impl::CovShiftForward(tmp, nu, Impl::CovShiftBackward(U[mu], mu, Impl::CovShiftIdentityBackward(U[nu], nu))), mu);

    // C2+
    tmp = U[mu] * Impl::ShiftStaple(adj(lambda), mu);
    out += Impl::ShiftStaple(Impl::CovShiftForward(U[nu], nu, Impl::CovShiftBackward(tmp, mu, Impl::CovShiftIdentityBackward(U[nu], nu))), mu);

    // C3+
    tmp = U[nu] * Impl::ShiftStaple(adj(lambda), nu);
    out += Impl::ShiftStaple(Impl::CovShiftForward(U[nu], nu, Impl::CovShiftBackward(U[mu], mu, Impl::CovShiftIdentityBackward(tmp, nu))), mu);

    // C4+
    out += Impl::ShiftStaple(Impl::CovShiftForward(U[nu], nu, Impl::CovShiftBackward(U[mu], mu, Impl::CovShiftIdentityBackward(U[nu], nu))), mu) * lambda;

    // insertion in lower staple
    // C1-
    out -= Impl::ShiftStaple(lambda, mu) * Impl::ShiftStaple(Impl::CovShiftBackward(U[nu], nu, Impl::CovShiftBackward(U[mu], mu, U[nu])), mu);

    // C2-
    tmp = adj(lambda) * U[nu];
    out -= Impl::ShiftStaple(Impl::CovShiftBackward(tmp, nu, Impl::CovShiftBackward(U[mu], mu, U[nu])), mu);

    // C3-
    tmp = lambda * U[nu];
    out -= Impl::ShiftStaple(Impl::CovShiftBackward(U[nu], nu, Impl::CovShiftBackward(U[mu], mu, tmp)), mu);

    // C4-
    out -= Impl::ShiftStaple(Impl::CovShiftBackward(U[nu], nu, Impl::CovShiftBackward(U[mu], mu, U[nu])), mu) * lambda;

    return out;
  }

  static CloverField fillCloverYZ(const GaugeLinkField &F)
  {
    CloverField T(F.Grid());
    T = Zero();
    autoView(T_v,T,AcceleratorWrite);
    autoView(F_v,F,AcceleratorRead);
    accelerator_for(i, T.Grid()->oSites(),1,
    {
      T_v[i]()(0, 1) = timesMinusI(F_v[i]()());
      T_v[i]()(1, 0) = timesMinusI(F_v[i]()());
      T_v[i]()(2, 3) = timesMinusI(F_v[i]()());
      T_v[i]()(3, 2) = timesMinusI(F_v[i]()());
    });

    return T;
  }

  static CloverField fillCloverXZ(const GaugeLinkField &F)
  {
    CloverField T(F.Grid());
    T = Zero();
    
    autoView(T_v, T,AcceleratorWrite);
    autoView(F_v, F,AcceleratorRead);
    accelerator_for(i, T.Grid()->oSites(),1,
    {
      T_v[i]()(0, 1) = -F_v[i]()();
      T_v[i]()(1, 0) = F_v[i]()();
      T_v[i]()(2, 3) = -F_v[i]()();
      T_v[i]()(3, 2) = F_v[i]()();
    });

    return T;
  }

  static CloverField fillCloverXY(const GaugeLinkField &F)
  {
    CloverField T(F.Grid());
    T = Zero();

    autoView(T_v,T,AcceleratorWrite);
    autoView(F_v,F,AcceleratorRead);
    accelerator_for(i, T.Grid()->oSites(),1,
    {
      T_v[i]()(0, 0) = timesMinusI(F_v[i]()());
      T_v[i]()(1, 1) = timesI(F_v[i]()());
      T_v[i]()(2, 2) = timesMinusI(F_v[i]()());
      T_v[i]()(3, 3) = timesI(F_v[i]()());
    });

    return T;
  }

  static CloverField fillCloverXT(const GaugeLinkField &F)
  {
    CloverField T(F.Grid());
    T = Zero();

    autoView( T_v , T, AcceleratorWrite);
    autoView( F_v , F, AcceleratorRead);
    accelerator_for(i, T.Grid()->oSites(),1,
    {
      T_v[i]()(0, 1) = timesI(F_v[i]()());
      T_v[i]()(1, 0) = timesI(F_v[i]()());
      T_v[i]()(2, 3) = timesMinusI(F_v[i]()());
      T_v[i]()(3, 2) = timesMinusI(F_v[i]()());
    });

    return T;
  }

  static CloverField fillCloverYT(const GaugeLinkField &F)
  {
    CloverField T(F.Grid());
    T = Zero();
    
    autoView( T_v ,T,AcceleratorWrite);
    autoView( F_v ,F,AcceleratorRead);
    accelerator_for(i, T.Grid()->oSites(),1,
    {
      T_v[i]()(0, 1) = -(F_v[i]()());
      T_v[i]()(1, 0) = (F_v[i]()());
      T_v[i]()(2, 3) = (F_v[i]()());
      T_v[i]()(3, 2) = -(F_v[i]()());
    });

    return T;
  }

  static CloverField fillCloverZT(const GaugeLinkField &F)
  {
    CloverField T(F.Grid());

    T = Zero();

    autoView( T_v , T,AcceleratorWrite);
    autoView( F_v , F,AcceleratorRead);
    accelerator_for(i, T.Grid()->oSites(),1,
    {
      T_v[i]()(0, 0) = timesI(F_v[i]()());
      T_v[i]()(1, 1) = timesMinusI(F_v[i]()());
      T_v[i]()(2, 2) = timesMinusI(F_v[i]()());
      T_v[i]()(3, 3) = timesI(F_v[i]()());
    });

    return T;
  }
};

NAMESPACE_END(Grid);
