/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/action/fermion/WilsonCloverFermionImplementation.h

    Copyright (C) 2017 - 2022

    Author: paboyle <paboyle@ph.ed.ac.uk>
    Author: Daniel Richtmann <daniel.richtmann@gmail.com>
    Author: Mattia Bruno <mattia.bruno@cern.ch>

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

#include <Grid/Grid.h>
#include <Grid/qcd/spin/Dirac.h>
#include <Grid/qcd/action/fermion/WilsonCloverHelpers.h>

////////////////////////////////////////////
// Standard Clover
//   (4+m0) + csw * clover_term
// Exp Clover
//   (4+m0) * exp(csw/(4+m0) clover_term)
//   = (4+m0) + csw * clover_term + ...
////////////////////////////////////////////

NAMESPACE_BEGIN(Grid);


//////////////////////////////////
// Generic Standard Clover
//////////////////////////////////

template<class Impl>
class CloverHelpers: public WilsonCloverHelpers<Impl> {
public:

  INHERIT_IMPL_TYPES(Impl);
  INHERIT_CLOVER_TYPES(Impl);

  typedef WilsonCloverHelpers<Impl> Helpers;

  static void Instantiate(CloverField& CloverTerm, CloverField& CloverTermInv, RealD csw_t, RealD diag_mass) {
    GridBase *grid = CloverTerm.Grid();
    CloverTerm += diag_mass;

    int lvol = grid->lSites();
    int DimRep = Impl::Dimension;
    {
      autoView(CTv,CloverTerm,CpuRead);
      autoView(CTIv,CloverTermInv,CpuWrite);
      thread_for(site, lvol, {
        Coordinate lcoor;
        grid->LocalIndexToLocalCoor(site, lcoor);
        Eigen::MatrixXcd EigenCloverOp = Eigen::MatrixXcd::Zero(Ns * DimRep, Ns * DimRep);
        Eigen::MatrixXcd EigenInvCloverOp = Eigen::MatrixXcd::Zero(Ns * DimRep, Ns * DimRep);
        typename SiteClover::scalar_object Qx = Zero(), Qxinv = Zero();
        peekLocalSite(Qx, CTv, lcoor);

        for (int j = 0; j < Ns; j++)
          for (int k = 0; k < Ns; k++)
            for (int a = 0; a < DimRep; a++)
              for (int b = 0; b < DimRep; b++){
                auto zz =  Qx()(j, k)(a, b);
                EigenCloverOp(a + j * DimRep, b + k * DimRep) = std::complex<double>(zz);
              }

        EigenInvCloverOp = EigenCloverOp.inverse();
        for (int j = 0; j < Ns; j++)
          for (int k = 0; k < Ns; k++)
            for (int a = 0; a < DimRep; a++)
              for (int b = 0; b < DimRep; b++)
                Qxinv()(j, k)(a, b) = EigenInvCloverOp(a + j * DimRep, b + k * DimRep);
               pokeLocalSite(Qxinv, CTIv, lcoor);
      });
    }
  }

  static GaugeLinkField Cmunu(std::vector<GaugeLinkField> &U, GaugeLinkField &lambda, int mu, int nu) {
    return Helpers::Cmunu(U, lambda, mu, nu);
  }

};


//////////////////////////////////
// Generic Exp Clover
//////////////////////////////////

template<class Impl>
class ExpCloverHelpers: public WilsonCloverHelpers<Impl> {
public:

  INHERIT_IMPL_TYPES(Impl);
  INHERIT_CLOVER_TYPES(Impl);

  template <typename vtype> using iImplClover = iScalar<iMatrix<iMatrix<vtype, Impl::Dimension>, Ns>>;
  typedef WilsonCloverHelpers<Impl> Helpers;

  // Can this be avoided?
  static void IdentityTimesC(const CloverField& in, RealD c) {
    int DimRep = Impl::Dimension;

    autoView(in_v, in, AcceleratorWrite);

    accelerator_for(ss, in.Grid()->oSites(), 1, {
      for (int sa=0; sa<Ns; sa++)
        for (int ca=0; ca<DimRep; ca++)
          in_v[ss]()(sa,sa)(ca,ca) = c;
    });
  }

  static int getNMAX(RealD prec, RealD R) {
    /* compute stop condition for exponential */
    int NMAX=1;
    RealD cond=R*R/2.;

    while (cond*std::exp(R)>prec) {
      NMAX++;
      cond*=R/(double)(NMAX+1);
    }
    return NMAX;
  }

  static int getNMAX(Lattice<iImplClover<vComplexD>> &t, RealD R) {return getNMAX(1e-12,R);}
  static int getNMAX(Lattice<iImplClover<vComplexF>> &t, RealD R) {return getNMAX(1e-6,R);}

  static void Instantiate(CloverField& Clover, CloverField& CloverInv, RealD csw_t, RealD diag_mass) {
    GridBase* grid = Clover.Grid();
    CloverField ExpClover(grid);

    int NMAX = getNMAX(Clover, 3.*csw_t/diag_mass);

    Clover *= (1.0/diag_mass);

    // Taylor expansion, slow but generic
    // Horner scheme: a0 + a1 x + a2 x^2 + .. = a0 + x (a1 + x(...))
    // qN = cN
    // qn = cn + qn+1 X
    std::vector<RealD> cn(NMAX+1);
    cn[0] = 1.0;
    for (int i=1; i<=NMAX; i++)
      cn[i] = cn[i-1] / RealD(i);

    ExpClover = Zero();
    IdentityTimesC(ExpClover, cn[NMAX]);
    for (int i=NMAX-1; i>=0; i--)
      ExpClover = ExpClover * Clover + cn[i];

    // prepare inverse
    CloverInv = (-1.0)*Clover;

    Clover = ExpClover * diag_mass;

    ExpClover = Zero();
    IdentityTimesC(ExpClover, cn[NMAX]);
    for (int i=NMAX-1; i>=0; i--)
      ExpClover = ExpClover * CloverInv + cn[i];

    CloverInv = ExpClover * (1.0/diag_mass);

  }

  static GaugeLinkField Cmunu(std::vector<GaugeLinkField> &U, GaugeLinkField &lambda, int mu, int nu) {
    assert(0);
    return lambda;
  }

};


//////////////////////////////////
// Compact Standard Clover
//////////////////////////////////


template<class Impl>
class CompactCloverHelpers: public CompactWilsonCloverHelpers<Impl>,
                            public WilsonCloverHelpers<Impl> {
public:

  INHERIT_IMPL_TYPES(Impl);
  INHERIT_CLOVER_TYPES(Impl);
  INHERIT_COMPACT_CLOVER_TYPES(Impl);

  typedef WilsonCloverHelpers<Impl> Helpers;
  typedef CompactWilsonCloverHelpers<Impl> CompactHelpers;

  static void InstantiateClover(CloverField& Clover, CloverField& CloverInv, RealD csw_t, RealD diag_mass) {
    Clover += diag_mass;
  }

  static void InvertClover(CloverField& InvClover,
                            const CloverDiagonalField& diagonal,
                            const CloverTriangleField& triangle,
                            CloverDiagonalField&       diagonalInv,
                            CloverTriangleField&       triangleInv,
                            bool fixedBoundaries) {

    CompactHelpers::Invert(diagonal, triangle, diagonalInv, triangleInv);
  }

  // TODO: implement Cmunu for better performances with compact layout, but don't do it
  // here, but rather in WilsonCloverHelpers.h -> CompactWilsonCloverHelpers
  static GaugeLinkField Cmunu(std::vector<GaugeLinkField> &U, GaugeLinkField &lambda, int mu, int nu) {
    return Helpers::Cmunu(U, lambda, mu, nu);
  }
};

//////////////////////////////////
// Compact Exp Clover
//////////////////////////////////

template<class Impl>
class CompactExpCloverHelpers: public CompactWilsonCloverHelpers<Impl> {
public:

  INHERIT_IMPL_TYPES(Impl);
  INHERIT_CLOVER_TYPES(Impl);
  INHERIT_COMPACT_CLOVER_TYPES(Impl);

  template <typename vtype> using iImplClover = iScalar<iMatrix<iMatrix<vtype, Impl::Dimension>, Ns>>;
  typedef CompactWilsonCloverHelpers<Impl> CompactHelpers;

  // Can this be avoided?
  static void IdentityTimesC(const CloverField& in, RealD c) {
    int DimRep = Impl::Dimension;

    autoView(in_v, in, AcceleratorWrite);

    accelerator_for(ss, in.Grid()->oSites(), 1, {
      for (int sa=0; sa<Ns; sa++)
        for (int ca=0; ca<DimRep; ca++)
          in_v[ss]()(sa,sa)(ca,ca) = c;
    });
  }

  static int getNMAX(RealD prec, RealD R) {
    /* compute stop condition for exponential */
    int NMAX=1;
    RealD cond=R*R/2.;

    while (cond*std::exp(R)>prec) {
      NMAX++;
      cond*=R/(double)(NMAX+1);
    }
    return NMAX;
  }

  static int getNMAX(Lattice<iImplClover<vComplexD>> &t, RealD R) {return getNMAX(1e-12,R);}
  static int getNMAX(Lattice<iImplClover<vComplexF>> &t, RealD R) {return getNMAX(1e-6,R);}

  static void InstantiateClover(CloverField& Clover, CloverField& CloverInv, RealD csw_t, RealD diag_mass) {

    GridBase* grid = Clover.Grid();
    CloverField ExpClover(grid);

    int NMAX = getNMAX(Clover, 3.*csw_t/diag_mass);

    Clover *= (1.0/diag_mass);

    // Taylor expansion, slow but generic
    // Horner scheme: a0 + a1 x + a2 x^2 + .. = a0 + x (a1 + x(...))
    // qN = cN
    // qn = cn + qn+1 X
    std::vector<RealD> cn(NMAX+1);
    cn[0] = 1.0;
    for (int i=1; i<=NMAX; i++)
      cn[i] = cn[i-1] / RealD(i);

    ExpClover = Zero();
    IdentityTimesC(ExpClover, cn[NMAX]);
    for (int i=NMAX-1; i>=0; i--)
      ExpClover = ExpClover * Clover + cn[i];

    // prepare inverse
    CloverInv = (-1.0)*Clover;

    Clover = ExpClover * diag_mass;

    ExpClover = Zero();
    IdentityTimesC(ExpClover, cn[NMAX]);
    for (int i=NMAX-1; i>=0; i--)
      ExpClover = ExpClover * CloverInv + cn[i];

    CloverInv = ExpClover * (1.0/diag_mass);

  }

  static void InvertClover(CloverField& InvClover,
                            const CloverDiagonalField& diagonal,
                            const CloverTriangleField& triangle,
                            CloverDiagonalField&       diagonalInv,
                            CloverTriangleField&       triangleInv,
                            bool fixedBoundaries) {

    if (fixedBoundaries)
    {
      CompactHelpers::Invert(diagonal, triangle, diagonalInv, triangleInv);
    }
    else
    {
      CompactHelpers::ConvertLayout(InvClover, diagonalInv, triangleInv);
    }
  }

  static GaugeLinkField Cmunu(std::vector<GaugeLinkField> &U, GaugeLinkField &lambda, int mu, int nu) {
    assert(0);
    return lambda;
  }

};


NAMESPACE_END(Grid);
