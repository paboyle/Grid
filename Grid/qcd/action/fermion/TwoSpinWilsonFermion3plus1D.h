/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/TwoSpinWilsonFermion3plus1D.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#pragma one 

NAMESPACE_BEGIN(Grid);

class TwoSpinWilsonFermion3plus1DStatic { 
public:
  // S-direction is INNERMOST and takes no part in the parity.
  static const std::vector<int> directions;
  static const std::vector<int> displacements;
  static constexpr int npoint = 6;
  static std::vector<int> MakeDirections(void);
  static std::vector<int> MakeDisplacements(void);
};

template<class Impl>
class TwoSpinWilsonFermion3plus1D : public TwoSpinWilsonKernels<Impl>, public TwoSpinWilsonFermion3plus1DStatic
{
public:
  INHERIT_IMPL_TYPES(Impl);
  typedef TwoSpinWilsonKernels<Impl> Kernels;

  FermionField _tmp;
  FermionField &tmp(void) { return _tmp; }

  int Dirichlet;
  Coordinate Block; 

  ///////////////////////////////////////////////////////////////
  // Implement the abstract base
  ///////////////////////////////////////////////////////////////
  GridBase *GaugeGrid(void)              { return _ThreeDimGrid ;}
  GridBase *GaugeRedBlackGrid(void)      { return _ThreeDimRedBlackGrid ;}
  GridBase *FermionGrid(void)            { return _FourDimGrid;}
  GridBase *FermionRedBlackGrid(void)    { return _FourDimRedBlackGrid;}

  // full checkerboard operations; leave unimplemented as abstract for now
  virtual void   M    (const FermionField &in, FermionField &out){assert(0);};
  virtual void   Mdag (const FermionField &in, FermionField &out){assert(0);};

  // half checkerboard operations; leave unimplemented as abstract for now
  virtual void   Meooe       (const FermionField &in, FermionField &out);
  virtual void   Mooee       (const FermionField &in, FermionField &out);
  virtual void   MooeeInv    (const FermionField &in, FermionField &out);

  virtual void   MeooeDag    (const FermionField &in, FermionField &out);
  virtual void   MooeeDag    (const FermionField &in, FermionField &out);
  virtual void   MooeeInvDag (const FermionField &in, FermionField &out);
  virtual void   Mdir   (const FermionField &in, FermionField &out,int dir,int disp){assert(0);};   // case by case Wilson, Clover, Cayley, ContFrac, PartFrac
  virtual void   MdirAll(const FermionField &in, std::vector<FermionField> &out){assert(0);};   // case by case Wilson, Clover, Cayley, ContFrac, PartFrac

  // These can be overridden by fancy 5d chiral action
  virtual void DhopDeriv  (GaugeField &mat,const FermionField &U,const FermionField &V,int dag);
  virtual void DhopDerivEO(GaugeField &mat,const FermionField &U,const FermionField &V,int dag);
  virtual void DhopDerivOE(GaugeField &mat,const FermionField &U,const FermionField &V,int dag);

  //  void MomentumSpacePropagatorHt_5d(FermionField &out,const FermionField &in,RealD mass,std::vector<double> twist) ;
  void MomentumSpacePropagatorHt(FermionField &out,const FermionField &in,RealD mass,std::vector<double> twist) ;
  void MomentumSpacePropagatorHw(FermionField &out,const FermionField &in,RealD mass,std::vector<double> twist) ;
  
  // Implement hopping term non-hermitian hopping term; half cb or both
  // Implement s-diagonal DW
  void DW    (const FermionField &in, FermionField &out,int dag);
  void Dhop  (const FermionField &in, FermionField &out,int dag);
  void DhopOE(const FermionField &in, FermionField &out,int dag);
  void DhopEO(const FermionField &in, FermionField &out,int dag);

  void DhopComms  (const FermionField &in, FermionField &out);
  void DhopCalc   (const FermionField &in, FermionField &out,uint64_t *ids);
  
  // add a DhopComm
  // -- suboptimal interface will presently trigger multiple comms.
  void DhopDir(const FermionField &in, FermionField &out,int dir,int disp);
  void DhopDirAll(const FermionField &in,std::vector<FermionField> &out);
  void DhopDirComms(const FermionField &in);
  void DhopDirCalc(const FermionField &in, FermionField &out,int point);
    
  ///////////////////////////////////////////////////////////////
  // New methods added 
  ///////////////////////////////////////////////////////////////
  void DerivInternal(StencilImpl & st,
		     DoubledGaugeField & U,
		     GaugeField &mat,
		     const FermionField &A,
		     const FermionField &B,
		     int dag);
    
  void DhopInternal(StencilImpl & st,
		    DoubledGaugeField &U,
		    const FermionField &in, 
		    FermionField &out,
		    int dag);

  void DhopInternalOverlappedComms(StencilImpl & st,
				   DoubledGaugeField &U,
				   const FermionField &in, 
				   FermionField &out,
				   int dag);

  void DhopInternalSerialComms(StencilImpl & st,
			       DoubledGaugeField &U,
			       const FermionField &in, 
			       FermionField &out,
			       int dag);
    
  // Constructors
  TwoSpinWilsonFermion3plus1D(GaugeField &_Umu,
		  GridCartesian         &FourDimGrid,
		  GridRedBlackCartesian &FourDimRedBlackGrid,
		  GridCartesian         &ThreeDimGrid,
		  GridRedBlackCartesian &ThreeDimRedBlackGrid,
		  double _M5,const ImplParams &p= ImplParams());

  virtual void DirichletBlock(const Coordinate & block)
  {
  }
    
  // DoubleStore
  void ImportGauge(const GaugeField &_Umu);
    
  ///////////////////////////////////////////////////////////////
  // Data members require to support the functionality
  ///////////////////////////////////////////////////////////////
public:
    
  // Add these to the support from Wilson
  GridBase *_ThreeDimGrid;
  GridBase *_ThreeDimRedBlackGrid;
  GridBase *_FourDimGrid;
  GridBase *_FourDimRedBlackGrid;
    
  double                        M5;
  int Ls;
    
  //Defines the stencils for even and odd
  StencilImpl Stencil; 
  StencilImpl StencilEven; 
  StencilImpl StencilOdd; 
    
  // Copy of the gauge field , with even and odd subsets
  DoubledGaugeField Umu;
  DoubledGaugeField UmuEven;
  DoubledGaugeField UmuOdd;

};

NAMESPACE_END(Grid);

