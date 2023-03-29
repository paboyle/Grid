/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/action/fermion/CayleyFermion5D.h

    Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
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

#include <Grid/qcd/action/fermion/WilsonFermion5D.h>

NAMESPACE_BEGIN(Grid);

template<class Impl>
class CayleyFermion5D : public WilsonFermion5D<Impl>
{
public:
  INHERIT_IMPL_TYPES(Impl);
public:

  // override multiply
  virtual void   M    (const FermionField &in, FermionField &out);
  virtual void   Mdag (const FermionField &in, FermionField &out);

  // half checkerboard operations
  virtual void   Meooe       (const FermionField &in, FermionField &out);
  virtual void   MeooeDag    (const FermionField &in, FermionField &out);
  virtual void   Mooee       (const FermionField &in, FermionField &out);
  virtual void   MooeeDag    (const FermionField &in, FermionField &out);
  virtual void   MooeeInv    (const FermionField &in, FermionField &out);
  virtual void   MooeeInvDag (const FermionField &in, FermionField &out);
  virtual void   Meo5D (const FermionField &psi, FermionField &chi);

  virtual void   M5D   (const FermionField &psi, FermionField &chi);
  virtual void   M5Ddag(const FermionField &psi, FermionField &chi);

  ///////////////////////////////////////////////////////////////
  // Physical surface field utilities
  ///////////////////////////////////////////////////////////////
  virtual void Dminus(const FermionField &psi, FermionField &chi);
  virtual void DminusDag(const FermionField &psi, FermionField &chi);
  virtual void ExportPhysicalFermionSolution(const FermionField &solution5d,FermionField &exported4d);
  virtual void ExportPhysicalFermionSource(const FermionField &solution5d, FermionField &exported4d);
  virtual void ImportPhysicalFermionSource(const FermionField &input4d,FermionField &imported5d);
  virtual void ImportUnphysicalFermion(const FermionField &solution5d, FermionField &exported4d);

  ///////////////////////////////////////////////////////////////
  // Support for MADWF tricks
  ///////////////////////////////////////////////////////////////
  RealD Mass(void) { return (mass_plus + mass_minus) / 2.0; };
  RealD MassPlus(void) { return mass_plus; };
  RealD MassMinus(void) { return mass_minus; };

  void  SetMass(RealD _mass) { 
    mass_plus=mass_minus=_mass; 
    SetCoefficientsInternal(_zolo_hi,_gamma,_b,_c);  // Reset coeffs
  } ;
  void  SetMass(RealD _mass_plus, RealD _mass_minus) { 
    mass_plus=_mass_plus;
    mass_minus=_mass_minus;
    SetCoefficientsInternal(_zolo_hi,_gamma,_b,_c);  // Reset coeffs
  } ;
  void  P(const FermionField &psi, FermionField &chi);
  void  Pdag(const FermionField &psi, FermionField &chi);
  
  /////////////////////////////////////////////////////
  // Instantiate different versions depending on Impl
  /////////////////////////////////////////////////////
  void M5D(const FermionField &psi,
	   const FermionField &phi,
	   FermionField &chi,
	   Vector<Coeff_t> &lower,
	   Vector<Coeff_t> &diag,
	   Vector<Coeff_t> &upper);

  void M5Ddag(const FermionField &psi,
	      const FermionField &phi,
	      FermionField &chi,
	      Vector<Coeff_t> &lower,
	      Vector<Coeff_t> &diag,
	      Vector<Coeff_t> &upper);

  virtual void   Instantiatable(void)=0;

  // force terms; five routines; default to Dhop on diagonal
  virtual void MDeriv  (GaugeField &mat,const FermionField &U,const FermionField &V,int dag);
  virtual void MoeDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag);
  virtual void MeoDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag);

  // Efficient support for multigrid coarsening
  virtual void  Mdir   (const FermionField &in, FermionField &out,int dir,int disp);
  virtual void  MdirAll(const FermionField &in, std::vector<FermionField> &out);

  void   Meooe5D       (const FermionField &in, FermionField &out);
  void   MeooeDag5D    (const FermionField &in, FermionField &out);

  //    protected:
  RealD mass_plus, mass_minus;

  // Save arguments to SetCoefficientsInternal
  Vector<Coeff_t> _gamma;
  RealD                _zolo_hi;
  RealD                _b;
  RealD                _c;

  // Cayley form Moebius (tanh and zolotarev)
  Vector<Coeff_t> omega;
  Vector<Coeff_t> bs;    // S dependent coeffs
  Vector<Coeff_t> cs;
  Vector<Coeff_t> as;
  // For preconditioning Cayley form
  Vector<Coeff_t> bee;
  Vector<Coeff_t> cee;
  Vector<Coeff_t> aee;
  Vector<Coeff_t> beo;
  Vector<Coeff_t> ceo;
  Vector<Coeff_t> aeo;
  // LDU factorisation of the eeoo matrix
  Vector<Coeff_t> lee;
  Vector<Coeff_t> leem;
  Vector<Coeff_t> uee;
  Vector<Coeff_t> ueem;
  Vector<Coeff_t> dee;

  // Matrices of 5d ee inverse params
  Vector<iSinglet<Simd> >  MatpInv;
  Vector<iSinglet<Simd> >  MatmInv;
  Vector<iSinglet<Simd> >  MatpInvDag;
  Vector<iSinglet<Simd> >  MatmInvDag;

  ///////////////////////////////////////////////////////////////
  // Conserved current utilities
  ///////////////////////////////////////////////////////////////

  // Virtual can't template
  void ContractConservedCurrent(PropagatorField &q_in_1,
				PropagatorField &q_in_2,
				PropagatorField &q_out,
				PropagatorField &phys_src,
				Current curr_type, 
				unsigned int mu);

   void SeqConservedCurrent(PropagatorField &q_in,
			   PropagatorField &q_out,
			   PropagatorField &phys_src,
			   Current curr_type,
			   unsigned int mu,
			   unsigned int tmin,
			   unsigned int tmax,
			   ComplexField &lattice_cmplx);

  void ContractJ5q(PropagatorField &q_in,ComplexField &J5q);
  void ContractJ5q(FermionField &q_in,ComplexField &J5q);

  ///////////////////////////////////////////////////////////////
  // Constructors
  ///////////////////////////////////////////////////////////////
  CayleyFermion5D(GaugeField &_Umu,
		  GridCartesian         &FiveDimGrid,
		  GridRedBlackCartesian &FiveDimRedBlackGrid,
		  GridCartesian         &FourDimGrid,
		  GridRedBlackCartesian &FourDimRedBlackGrid,
		  RealD _mass,RealD _M5,const ImplParams &p= ImplParams());


protected:
  virtual void SetCoefficientsZolotarev(RealD zolohi,Approx::zolotarev_data *zdata,RealD b,RealD c);
  virtual void SetCoefficientsTanh(Approx::zolotarev_data *zdata,RealD b,RealD c);
  virtual void SetCoefficientsInternal(RealD zolo_hi,Vector<Coeff_t> & gamma,RealD b,RealD c);
};

NAMESPACE_END(Grid);

