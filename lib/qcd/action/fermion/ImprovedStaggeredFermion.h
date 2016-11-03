/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/ImprovedStaggered.h

Copyright (C) 2015

Author: Azusa Yamaguchi, Peter Boyle

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
#ifndef GRID_QCD_IMPR_STAG_FERMION_H
#define GRID_QCD_IMPR_STAG_FERMION_H

namespace Grid {

namespace QCD {

class ImprovedStaggeredFermionStatic {
 public:
  static const std::vector<int> directions;
  static const std::vector<int> displacements;
  static const int npoint = 16;
};

template <class Impl>
class ImprovedStaggeredFermion : public StaggeredKernels<Impl>, public ImprovedStaggeredFermionStatic {
 public:
  INHERIT_IMPL_TYPES(Impl);
  typedef StaggeredKernels<Impl> Kernels;

  ///////////////////////////////////////////////////////////////
  // Implement the abstract base
  ///////////////////////////////////////////////////////////////
  GridBase *GaugeGrid(void) { return _grid; }
  GridBase *GaugeRedBlackGrid(void) { return _cbgrid; }
  GridBase *FermionGrid(void) { return _grid; }
  GridBase *FermionRedBlackGrid(void) { return _cbgrid; }

  //////////////////////////////////////////////////////////////////
  // override multiply; cut number routines if pass dagger argument
  // and also make interface more uniformly consistent
  //////////////////////////////////////////////////////////////////
  RealD M(const FermionField &in, FermionField &out);
  RealD Mdag(const FermionField &in, FermionField &out);

  /////////////////////////////////////////////////////////
  // half checkerboard operations
  /////////////////////////////////////////////////////////
  void Meooe(const FermionField &in, FermionField &out);
  void MeooeDag(const FermionField &in, FermionField &out);

  // allow override for twisted mass and clover
  virtual void Mooee(const FermionField &in, FermionField &out);
  virtual void MooeeDag(const FermionField &in, FermionField &out);
  virtual void MooeeInv(const FermionField &in, FermionField &out);
  virtual void MooeeInvDag(const FermionField &in, FermionField &out);

  ////////////////////////
  // Derivative interface
  ////////////////////////
  // Interface calls an internal routine
  void DhopDeriv  (GaugeField &mat, const FermionField &U, const FermionField &V, int dag);
  void DhopDerivOE(GaugeField &mat, const FermionField &U, const FermionField &V, int dag);
  void DhopDerivEO(GaugeField &mat, const FermionField &U, const FermionField &V, int dag);

  ///////////////////////////////////////////////////////////////
  // non-hermitian hopping term; half cb or both
  ///////////////////////////////////////////////////////////////
  void Dhop  (const FermionField &in, FermionField &out, int dag);
  void DhopOE(const FermionField &in, FermionField &out, int dag);
  void DhopEO(const FermionField &in, FermionField &out, int dag);

  ///////////////////////////////////////////////////////////////
  // Multigrid assistance; force term uses too
  ///////////////////////////////////////////////////////////////
  void Mdir(const FermionField &in, FermionField &out, int dir, int disp);
  void DhopDir(const FermionField &in, FermionField &out, int dir, int disp);

  ///////////////////////////////////////////////////////////////
  // Extra methods added by derived
  ///////////////////////////////////////////////////////////////
  void DerivInternal(StencilImpl &st, 
		     DoubledGaugeField &U,DoubledGaugeField &UUU,
		     GaugeField &mat, 
		     const FermionField &A, const FermionField &B, int dag);

  void DhopInternal(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U,DoubledGaugeField &UUU,
                    const FermionField &in, FermionField &out, int dag);

  // Constructor
  ImprovedStaggeredFermion(GaugeField &_Umu, GridCartesian &Fgrid,
			   GridRedBlackCartesian &Hgrid, RealD _mass,
			   RealD _c1=9.0/8.0, RealD _c2=-1.0/24.0,RealD _u0,
			   const ImplParams &p = ImplParams());

  // DoubleStore impl dependent
  void ImportGauge(const GaugeField &_Umu);

  ///////////////////////////////////////////////////////////////
  // Data members require to support the functionality
  ///////////////////////////////////////////////////////////////

  //    protected:
 public:
  // any other parameters of action ???

  RealD mass;
  RealD u0;
  RealD c1;
  RealD c2;

  GridBase *_grid;
  GridBase *_cbgrid;

  // Defines the stencils for even and odd
  StencilImpl Stencil;
  StencilImpl StencilEven;
  StencilImpl StencilOdd;

  // Copy of the gauge field , with even and odd subsets
  DoubledGaugeField Umu;
  DoubledGaugeField UmuEven;
  DoubledGaugeField UmuOdd;

  DoubledGaugeField UUUmu;
  DoubledGaugeField UUUmuEven;
  DoubledGaugeField UUUmuOdd;

  LebesgueOrder Lebesgue;
  LebesgueOrder LebesgueEvenOdd;
};

typedef ImprovedStaggeredFermion<StaggeredImplF> ImprovedStaggeredFermionF;
typedef ImprovedStaggeredFermion<StaggeredImplD> ImprovedStaggeredFermionD;

}
}
#endif
