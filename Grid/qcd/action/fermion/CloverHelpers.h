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

  static int getNMAX(Lattice<iImplClover<vComplexD2>> &t, RealD R) {return getNMAX(1e-12,R);}
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

  static void MassTerm(CloverField& Clover, RealD diag_mass) {
    Clover += diag_mass;
  }

  static void Exponentiate_Clover(CloverDiagonalField& Diagonal,
                          CloverTriangleField& Triangle,
                          RealD csw_t, RealD diag_mass) {

    // Do nothing
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

  static void MassTerm(CloverField& Clover, RealD diag_mass) {
    // do nothing!
    // mass term is multiplied to exp(Clover) below
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

  static int getNMAX(Lattice<iImplCloverDiagonal<vComplexD>> &t, RealD R) {return getNMAX(1e-12,R);}
  static int getNMAX(Lattice<iImplCloverDiagonal<vComplexF>> &t, RealD R) {return getNMAX(1e-6,R);}

  static void ExponentiateHermitean6by6(const iMatrix<ComplexD,6> &arg, const RealD& alpha, const std::vector<RealD>& cN, const int Niter, iMatrix<ComplexD,6>& dest){

  	  typedef iMatrix<ComplexD,6> mat;

  	  RealD qn[6];
  	  RealD qnold[6];
  	  RealD p[5];
  	  RealD trA2, trA3, trA4;

  	  mat A2, A3, A4, A5;
  	  A2 = alpha * alpha * arg * arg;
  	  A3 = alpha * arg * A2;
  	  A4 = A2 * A2;
  	  A5 = A2 * A3;

  	  trA2 = toReal( trace(A2) );
  	  trA3 = toReal( trace(A3) );
  	  trA4 = toReal( trace(A4));

  	  p[0] = toReal( trace(A3 * A3)) / 6.0 - 0.125 * trA4 * trA2 - trA3 * trA3 / 18.0 + trA2 * trA2 * trA2/ 48.0;
  	  p[1] = toReal( trace(A5)) / 5.0 - trA3 * trA2 / 6.0;
  	  p[2] = toReal( trace(A4)) / 4.0 - 0.125 * trA2 * trA2;
  	  p[3] = trA3 / 3.0;
  	  p[4] = 0.5 * trA2;

  	  qnold[0] = cN[Niter];
  	  qnold[1] = 0.0;
  	  qnold[2] = 0.0;
  	  qnold[3] = 0.0;
  	  qnold[4] = 0.0;
  	  qnold[5] = 0.0;

  	  for(int i = Niter-1; i >= 0; i--)
  	  {
  	   qn[0] = p[0] * qnold[5] + cN[i];
  	   qn[1] = p[1] * qnold[5] + qnold[0];
  	   qn[2] = p[2] * qnold[5] + qnold[1];
  	   qn[3] = p[3] * qnold[5] + qnold[2];
  	   qn[4] = p[4] * qnold[5] + qnold[3];
  	   qn[5] = qnold[4];

  	   qnold[0] = qn[0];
  	   qnold[1] = qn[1];
  	   qnold[2] = qn[2];
  	   qnold[3] = qn[3];
  	   qnold[4] = qn[4];
  	   qnold[5] = qn[5];
  	  }

  	  mat unit(1.0);

  	  dest = (qn[0] * unit + qn[1] * alpha * arg + qn[2] * A2 + qn[3] * A3 + qn[4] * A4 + qn[5] * A5);

    }

  static void Exponentiate_Clover(CloverDiagonalField& Diagonal, CloverTriangleField& Triangle, RealD csw_t, RealD diag_mass) {

    GridBase* grid = Diagonal.Grid();
    int NMAX = getNMAX(Diagonal, 3.*csw_t/diag_mass);

    //
    // Implementation completely in Daniel's layout
    //

    // Taylor expansion with Cayley-Hamilton recursion
    // underlying Horner scheme as above
    std::vector<RealD> cn(NMAX+1);
    cn[0] = 1.0;
    for (int i=1; i<=NMAX; i++){
      cn[i] = cn[i-1] / RealD(i);
    }

      // Taken over from Daniel's implementation
      conformable(Diagonal, Triangle);

      long lsites = grid->lSites();
    {
      typedef typename SiteCloverDiagonal::scalar_object scalar_object_diagonal;
      typedef typename SiteCloverTriangle::scalar_object scalar_object_triangle;
      typedef iMatrix<ComplexD,6> mat;

      autoView(diagonal_v,  Diagonal,  CpuRead);
      autoView(triangle_v,  Triangle,  CpuRead);
      autoView(diagonalExp_v, Diagonal, CpuWrite);
      autoView(triangleExp_v, Triangle, CpuWrite);

      thread_for(site, lsites, { // NOTE: Not on GPU because of (peek/poke)LocalSite

    	  mat srcCloverOpUL(0.0); // upper left block
    	  mat srcCloverOpLR(0.0); // lower right block
    	  mat ExpCloverOp;

        scalar_object_diagonal diagonal_tmp     = Zero();
        scalar_object_diagonal diagonal_exp_tmp = Zero();
        scalar_object_triangle triangle_tmp     = Zero();
        scalar_object_triangle triangle_exp_tmp = Zero();

        Coordinate lcoor;
        grid->LocalIndexToLocalCoor(site, lcoor);

        peekLocalSite(diagonal_tmp, diagonal_v, lcoor);
        peekLocalSite(triangle_tmp, triangle_v, lcoor);

        int block;
        block = 0;
        for(int i = 0; i < 6; i++){
        	for(int j = 0; j < 6; j++){
        		if (i == j){
        			srcCloverOpUL(i,j) = static_cast<ComplexD>(TensorRemove(diagonal_tmp()(block)(i)));
        		}
        		else{
        			srcCloverOpUL(i,j) = static_cast<ComplexD>(TensorRemove(CompactHelpers::triangle_elem(triangle_tmp, block, i, j)));
        		}
        	}
        }
        block = 1;
        for(int i = 0; i < 6; i++){
          	for(int j = 0; j < 6; j++){
           		if (i == j){
           			srcCloverOpLR(i,j) = static_cast<ComplexD>(TensorRemove(diagonal_tmp()(block)(i)));
           		}
           		else{
           			srcCloverOpLR(i,j) = static_cast<ComplexD>(TensorRemove(CompactHelpers::triangle_elem(triangle_tmp, block, i, j)));
           		}
            }
        }

        // exp(Clover)

        ExponentiateHermitean6by6(srcCloverOpUL,1.0/diag_mass,cn,NMAX,ExpCloverOp);

        block = 0;
        for(int i = 0; i < 6; i++){
        	for(int j = 0; j < 6; j++){
            	if (i == j){
            		diagonal_exp_tmp()(block)(i) = ExpCloverOp(i,j);
            	}
            	else if(i < j){
            		triangle_exp_tmp()(block)(CompactHelpers::triangle_index(i, j)) = ExpCloverOp(i,j);
            	}
           	}
        }

        ExponentiateHermitean6by6(srcCloverOpLR,1.0/diag_mass,cn,NMAX,ExpCloverOp);

        block = 1;
        for(int i = 0; i < 6; i++){
        	for(int j = 0; j < 6; j++){
              	if (i == j){
              		diagonal_exp_tmp()(block)(i) = ExpCloverOp(i,j);
               	}
               	else if(i < j){
               		triangle_exp_tmp()(block)(CompactHelpers::triangle_index(i, j)) = ExpCloverOp(i,j);
               	}
            }
        }

        pokeLocalSite(diagonal_exp_tmp, diagonalExp_v, lcoor);
        pokeLocalSite(triangle_exp_tmp, triangleExp_v, lcoor);
      });
    }

    Diagonal *= diag_mass;
    Triangle *= diag_mass;
  }


  static GaugeLinkField Cmunu(std::vector<GaugeLinkField> &U, GaugeLinkField &lambda, int mu, int nu) {
    assert(0);
    return lambda;
  }

};


NAMESPACE_END(Grid);
