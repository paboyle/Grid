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

  FermionField _tmp;
  FermionField &tmp(void) { return _tmp; }

  ////////////////////////////////////////
  // Performance monitoring
  ////////////////////////////////////////
  void Report(void);
  void ZeroCounters(void);
  double DhopTotalTime;
  double DhopCalls;
  double DhopCommTime;
  double DhopComputeTime;
  double DhopComputeTime2;
  double DhopFaceTime;

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
  void Mooee(const FermionField &in, FermionField &out);
  void MooeeDag(const FermionField &in, FermionField &out);
  void MooeeInv(const FermionField &in, FermionField &out);
  void MooeeInvDag(const FermionField &in, FermionField &out);

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
  void DhopInternalSerialComms(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U,DoubledGaugeField &UUU,
                    const FermionField &in, FermionField &out, int dag);
  void DhopInternalOverlappedComms(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U,DoubledGaugeField &UUU,
                    const FermionField &in, FermionField &out, int dag);

  //////////////////////////////////////////////////////////////////////////
  // Grid own interface Constructor
  //////////////////////////////////////////////////////////////////////////
  ImprovedStaggeredFermion(GaugeField &_Uthin, GaugeField &_Ufat, GridCartesian &Fgrid,
			   GridRedBlackCartesian &Hgrid, RealD _mass,
			   RealD _c1, RealD _c2,RealD _u0,
			   const ImplParams &p = ImplParams());

  //////////////////////////////////////////////////////////////////////////
  // MILC constructor no gauge fields
  //////////////////////////////////////////////////////////////////////////
  ImprovedStaggeredFermion(GridCartesian &Fgrid, GridRedBlackCartesian &Hgrid, RealD _mass,
			   RealD _c1=1.0, RealD _c2=1.0,RealD _u0=1.0,
			   const ImplParams &p = ImplParams());

  // DoubleStore impl dependent
  void ImportGauge      (const GaugeField &_Uthin ) { assert(0); }
  void ImportGauge      (const GaugeField &_Uthin  ,const GaugeField &_Ufat);
  void ImportGaugeSimple(const GaugeField &_UUU    ,const GaugeField &_U);
  void ImportGaugeSimple(const DoubledGaugeField &_UUU,const DoubledGaugeField &_U);
  DoubledGaugeField &GetU(void)   { return Umu ; } ;
  DoubledGaugeField &GetUUU(void) { return UUUmu; };
  void CopyGaugeCheckerboards(void);

  ///////////////////////////////////////////////////////////////
  // Data members require to support the functionality
  ///////////////////////////////////////////////////////////////

  //    protected:
 public:
  // any other parameters of action ???
  virtual int   isTrivialEE(void) { return 1; };
  virtual RealD Mass(void) { return mass; }
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
  
  ///////////////////////////////////////////////////////////////
  // Conserved current utilities
  ///////////////////////////////////////////////////////////////
  void ContractConservedCurrent(PropagatorField &q_in_1,
                                PropagatorField &q_in_2,
                                PropagatorField &q_out,
                                Current curr_type,
                                unsigned int mu);
  void SeqConservedCurrent(PropagatorField &q_in, 
                           PropagatorField &q_out,
                           Current curr_type, 
                           unsigned int mu,
                           unsigned int tmin, 
                           unsigned int tmax,
			   ComplexField &lattice_cmplx);
};

typedef ImprovedStaggeredFermion<StaggeredImplF> ImprovedStaggeredFermionF;
typedef ImprovedStaggeredFermion<StaggeredImplD> ImprovedStaggeredFermionD;

}
}
#endif
