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
#ifndef GRID_QCD_NAIVE_STAG_FERMION_H
#define GRID_QCD_NAIVE_STAG_FERMION_H

NAMESPACE_BEGIN(Grid);

class NaiveStaggeredFermionStatic {
public:
  static const std::vector<int> directions;
  static const std::vector<int> displacements;
  static const int npoint = 8;
};

template <class Impl>
class NaiveStaggeredFermion : public StaggeredKernels<Impl>, public NaiveStaggeredFermionStatic {
public:
  INHERIT_IMPL_TYPES(Impl);
  typedef StaggeredKernels<Impl> Kernels;

  FermionField _tmp;
  FermionField &tmp(void) { return _tmp; }

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
  void M(const FermionField &in, FermionField &out);
  void Mdag(const FermionField &in, FermionField &out);

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
  void MdirAll(const FermionField &in, std::vector<FermionField> &out);
  void DhopDir(const FermionField &in, FermionField &out, int dir, int disp);

  ///////////////////////////////////////////////////////////////
  // Extra methods added by derived
  ///////////////////////////////////////////////////////////////
  void DerivInternal(StencilImpl &st, 
		     DoubledGaugeField &U,
		     GaugeField &mat, 
		     const FermionField &A, const FermionField &B, int dag);

  void DhopInternal(StencilImpl &st, DoubledGaugeField &U,
                    const FermionField &in, FermionField &out, int dag);
  void DhopInternalSerialComms(StencilImpl &st, DoubledGaugeField &U,
			       const FermionField &in, FermionField &out, int dag);
  void DhopInternalOverlappedComms(StencilImpl &st, DoubledGaugeField &U,
				   const FermionField &in, FermionField &out, int dag);

  //////////////////////////////////////////////////////////////////////////
  // Grid own interface Constructor
  //////////////////////////////////////////////////////////////////////////
  NaiveStaggeredFermion(GaugeField &_U, GridCartesian &Fgrid,
			GridRedBlackCartesian &Hgrid, RealD _mass,
			RealD _c1, RealD _u0,
			const ImplParams &p = ImplParams());
  NaiveStaggeredFermion(GridCartesian &Fgrid,
			GridRedBlackCartesian &Hgrid, RealD _mass,
			RealD _c1, RealD _u0,
			const ImplParams &p = ImplParams());

  // DoubleStore impl dependent
  void ImportGauge      (const GaugeField &_U );
  DoubledGaugeField &GetU(void)   { return Umu ; } ;
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

  ///////////////////////////////////////////////////////////////
  // Conserved current utilities
  ///////////////////////////////////////////////////////////////
  void ContractConservedCurrent(PropagatorField &q_in_1,
                                PropagatorField &q_in_2,
                                PropagatorField &q_out,
                                PropagatorField &src,
                                Current curr_type,
                                unsigned int mu);
  void SeqConservedCurrent(PropagatorField &q_in,
                           PropagatorField &q_out,
                           PropagatorField &srct,
                           Current curr_type,
                           unsigned int mu, 
                           unsigned int tmin,
                           unsigned int tmax,
			   ComplexField &lattice_cmplx);
};

typedef NaiveStaggeredFermion<StaggeredImplF> NaiveStaggeredFermionF;
typedef NaiveStaggeredFermion<StaggeredImplD> NaiveStaggeredFermionD;

NAMESPACE_END(Grid);

#endif
