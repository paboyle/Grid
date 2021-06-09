/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/DirichletFermionOperator.h

    Copyright (C) 2021

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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

////////////////////////////////////////////////////////////////
// Wrap a fermion operator in Dirichlet BC's at node boundary
////////////////////////////////////////////////////////////////
    
template<class Impl>
class DirichletFermionOperator : public FermionOperator<Impl>
{
public:

  INHERIT_IMPL_TYPES(Impl);

  // Data members
  int CommsMode;
  Coordinate Block;
  DirichletFilter<GaugeField> Filter;
  FermionOperator<Impl> & FermOp;
  
  // Constructor / bespoke
  DirichletFermionOperator(FermionOperator<Impl> & _FermOp, Coordinate &_Block)
    : FermOp(_FermOp), Block(_Block), Filter(Block)
  {
    // Save what the comms mode should be under normal BCs
    CommsMode = WilsonKernelsStatic::Comms;
    assert((CommsMode == WilsonKernelsStatic::CommsAndCompute)
         ||(CommsMode == WilsonKernelsStatic::CommsThenCompute));

    // Check the block size divides local lattice
    GridBase *grid = FermOp.GaugeGrid();

    int blocks_per_rank = 1;
    Coordinate LocalDims = grid->LocalDimensions();
    Coordinate GlobalDims= grid->GlobalDimensions();
    assert(Block.size()==LocalDims.size());

    for(int d=0;d<LocalDims.size();d++){
      if (Block[d]&&(Block[d]<=GlobalDims[d])){
	int r = LocalDims[d] % Block[d];
	assert(r == 0);
	blocks_per_rank *= (LocalDims[d] / Block[d]);
      }
    }
    // Even blocks per node required // could be relaxed but inefficient use of hardware as idle nodes in boundary operator R
    assert( blocks_per_rank != 0);

    // Possible checks that SIMD lanes are used with full occupancy???
  };
  virtual ~DirichletFermionOperator(void) = default;

  void DirichletOn(void)   {
    assert(WilsonKernelsStatic::Comms!= WilsonKernelsStatic::CommsDirichlet);
    //    WilsonKernelsStatic::Comms = WilsonKernelsStatic::CommsDirichlet;
  }
  void DirichletOff(void)  {
    //    assert(WilsonKernelsStatic::Comms== WilsonKernelsStatic::CommsDirichlet);
    //    WilsonKernelsStatic::Comms = CommsMode;
  }

  // Implement the full interface
  virtual FermionField &tmp(void) { return FermOp.tmp(); };

  virtual GridBase *FermionGrid(void)         { return FermOp.FermionGrid(); }
  virtual GridBase *FermionRedBlackGrid(void) { return FermOp.FermionRedBlackGrid(); }
  virtual GridBase *GaugeGrid(void)           { return FermOp.GaugeGrid(); }
  virtual GridBase *GaugeRedBlackGrid(void)   { return FermOp.GaugeRedBlackGrid(); }
  
  // override multiply
  virtual void  M    (const FermionField &in, FermionField &out) { DirichletOn(); FermOp.M(in,out);    DirichletOff();  };
  virtual void  Mdag (const FermionField &in, FermionField &out) { DirichletOn(); FermOp.Mdag(in,out); DirichletOff();  };

  // half checkerboard operaions
  virtual void   Meooe       (const FermionField &in, FermionField &out) { DirichletOn(); FermOp.Meooe(in,out);    DirichletOff(); };  
  virtual void   MeooeDag    (const FermionField &in, FermionField &out) { DirichletOn(); FermOp.MeooeDag(in,out); DirichletOff(); };
  virtual void   Mooee       (const FermionField &in, FermionField &out) { DirichletOn(); FermOp.Mooee(in,out);    DirichletOff(); };
  virtual void   MooeeDag    (const FermionField &in, FermionField &out) { DirichletOn(); FermOp.MooeeDag(in,out); DirichletOff(); };
  virtual void   MooeeInv    (const FermionField &in, FermionField &out) { DirichletOn(); FermOp.MooeeInv(in,out); DirichletOff(); };
  virtual void   MooeeInvDag (const FermionField &in, FermionField &out) { DirichletOn(); FermOp.MooeeInvDag(in,out); DirichletOff(); };

  // non-hermitian hopping term; half cb or both
  virtual void Dhop  (const FermionField &in, FermionField &out,int dag) { DirichletOn(); FermOp.Dhop(in,out,dag);    DirichletOff(); };
  virtual void DhopOE(const FermionField &in, FermionField &out,int dag) { DirichletOn(); FermOp.DhopOE(in,out,dag);  DirichletOff(); };
  virtual void DhopEO(const FermionField &in, FermionField &out,int dag) { DirichletOn(); FermOp.DhopEO(in,out,dag);  DirichletOff(); };
  virtual void DhopDir(const FermionField &in, FermionField &out,int dir,int disp) { DirichletOn(); FermOp.DhopDir(in,out,dir,disp);  DirichletOff(); };

  // force terms; five routines; default to Dhop on diagonal
  virtual void MDeriv  (GaugeField &mat,const FermionField &U,const FermionField &V,int dag){FermOp.MDeriv(mat,U,V,dag);};
  virtual void MoeDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag){FermOp.MoeDeriv(mat,U,V,dag);};
  virtual void MeoDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag){FermOp.MeoDeriv(mat,U,V,dag);};
  virtual void MooDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag){FermOp.MooDeriv(mat,U,V,dag);};
  virtual void MeeDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag){FermOp.MeeDeriv(mat,U,V,dag);};

  virtual void DhopDeriv  (GaugeField &mat,const FermionField &U,const FermionField &V,int dag){FermOp.DhopDeriv(mat,U,V,dag);};
  virtual void DhopDerivEO(GaugeField &mat,const FermionField &U,const FermionField &V,int dag){FermOp.DhopDerivEO(mat,U,V,dag);};
  virtual void DhopDerivOE(GaugeField &mat,const FermionField &U,const FermionField &V,int dag){FermOp.DhopDerivOE(mat,U,V,dag);};

  virtual void  Mdiag  (const FermionField &in, FermionField &out) { Mooee(in,out);};
  virtual void  Mdir   (const FermionField &in, FermionField &out,int dir,int disp){FermOp.Mdir(in,out,dir,disp);};
  virtual void  MdirAll(const FermionField &in, std::vector<FermionField> &out)    {FermOp.MdirAll(in,out);};

  ///////////////////////////////////////////////
  // Updates gauge field during HMC
  ///////////////////////////////////////////////
  DoubledGaugeField &GetDoubledGaugeField(void){ return FermOp.GetDoubledGaugeField(); };
  DoubledGaugeField &GetDoubledGaugeFieldE(void){ return FermOp.GetDoubledGaugeFieldE(); };
  DoubledGaugeField &GetDoubledGaugeFieldO(void){ return FermOp.GetDoubledGaugeFieldO(); };
  virtual void ImportGauge(const GaugeField & _U)
  {
    GaugeField U = _U;
    // Filter gauge field to apply Dirichlet
    Filter.applyFilter(U);
    FermOp.ImportGauge(U);
  }
  ///////////////////////////////////////////////
  // Physical field import/export
  ///////////////////////////////////////////////
  virtual void Dminus(const FermionField &psi, FermionField &chi)    { FermOp.Dminus(psi,chi); }
  virtual void DminusDag(const FermionField &psi, FermionField &chi) { FermOp.DminusDag(psi,chi); }
  virtual void ImportFourDimPseudoFermion(const FermionField &input,FermionField &imported)   { FermOp.ImportFourDimPseudoFermion(input,imported);}
  virtual void ExportFourDimPseudoFermion(const FermionField &solution,FermionField &exported){ FermOp.ExportFourDimPseudoFermion(solution,exported);}
  virtual void ImportPhysicalFermionSource(const FermionField &input,FermionField &imported)  { FermOp.ImportPhysicalFermionSource(input,imported);}
  virtual void ImportUnphysicalFermion(const FermionField &input,FermionField &imported)      { FermOp.ImportUnphysicalFermion(input,imported);}
  virtual void ExportPhysicalFermionSolution(const FermionField &solution,FermionField &exported) {FermOp.ExportPhysicalFermionSolution(solution,exported);}
  virtual void ExportPhysicalFermionSource(const FermionField &solution,FermionField &exported)   {FermOp.ExportPhysicalFermionSource(solution,exported);}
  //////////////////////////////////////////////////////////////////////
  // Should never be used
  //////////////////////////////////////////////////////////////////////
  virtual void MomentumSpacePropagator(FermionField &out,const FermionField &in,RealD _m,std::vector<double> twist) { assert(0);};
  virtual void FreePropagator(const FermionField &in,FermionField &out,RealD mass,std::vector<Complex> boundary,std::vector<double> twist) {assert(0);}
  virtual void FreePropagator(const FermionField &in,FermionField &out,RealD mass) { assert(0);}
  virtual void ContractConservedCurrent(PropagatorField &q_in_1,
					PropagatorField &q_in_2,
					PropagatorField &q_out,
					PropagatorField &phys_src,
					Current curr_type,
					unsigned int mu)
  {assert(0);};
  virtual void SeqConservedCurrent(PropagatorField &q_in, 
				   PropagatorField &q_out,
				   PropagatorField &phys_src,
				   Current curr_type,
				   unsigned int mu,
				   unsigned int tmin, 
				   unsigned int tmax,
				   ComplexField &lattice_cmplx)
  {assert(0);};
      // Only reimplemented in Wilson5D 
      // Default to just a zero correlation function
  virtual void ContractJ5q(FermionField &q_in   ,ComplexField &J5q) { J5q=Zero(); };
  virtual void ContractJ5q(PropagatorField &q_in,ComplexField &J5q) { J5q=Zero(); };
  
};


NAMESPACE_END(Grid);

