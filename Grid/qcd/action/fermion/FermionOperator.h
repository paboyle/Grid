    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/FermionOperator.h

    Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: Vera Guelpers <V.M.Guelpers@soton.ac.uk>

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
#ifndef  GRID_QCD_FERMION_OPERATOR_H
#define  GRID_QCD_FERMION_OPERATOR_H

namespace Grid {

  namespace QCD {

    ////////////////////////////////////////////////////////////////
    // Allow to select  between gauge representation rank bc's, flavours etc.
    // and single/double precision.
    ////////////////////////////////////////////////////////////////
    
    template<class Impl>
    class FermionOperator : public CheckerBoardedSparseMatrixBase<typename Impl::FermionField>, public Impl
    {
    public:

      INHERIT_IMPL_TYPES(Impl);

      FermionOperator(const ImplParams &p= ImplParams()) : Impl(p) {};
      virtual ~FermionOperator(void) = default;

      virtual FermionField &tmp(void) = 0;

      GridBase * Grid(void)   { return FermionGrid(); };   // this is all the linalg routines need to know
      GridBase * RedBlackGrid(void) { return FermionRedBlackGrid(); };

      virtual GridBase *FermionGrid(void)         =0;
      virtual GridBase *FermionRedBlackGrid(void) =0;
      virtual GridBase *GaugeGrid(void)           =0;
      virtual GridBase *GaugeRedBlackGrid(void)   =0;

      // override multiply
      virtual RealD  M    (const FermionField &in, FermionField &out)=0;
      virtual RealD  Mdag (const FermionField &in, FermionField &out)=0;

      // half checkerboard operaions
      virtual void   Meooe       (const FermionField &in, FermionField &out)=0;
      virtual void   MeooeDag    (const FermionField &in, FermionField &out)=0;
      virtual void   Mooee       (const FermionField &in, FermionField &out)=0;
      virtual void   MooeeDag    (const FermionField &in, FermionField &out)=0;
      virtual void   MooeeInv    (const FermionField &in, FermionField &out)=0;
      virtual void   MooeeInvDag (const FermionField &in, FermionField &out)=0;

      // non-hermitian hopping term; half cb or both
      virtual void Dhop  (const FermionField &in, FermionField &out,int dag)=0;
      virtual void DhopOE(const FermionField &in, FermionField &out,int dag)=0;
      virtual void DhopEO(const FermionField &in, FermionField &out,int dag)=0;
      virtual void DhopDir(const FermionField &in, FermionField &out,int dir,int disp)=0; // implemented by WilsonFermion and WilsonFermion5D

      // force terms; five routines; default to Dhop on diagonal
      virtual void MDeriv  (GaugeField &mat,const FermionField &U,const FermionField &V,int dag){DhopDeriv(mat,U,V,dag);};
      virtual void MoeDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag){DhopDerivOE(mat,U,V,dag);};
      virtual void MeoDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag){DhopDerivEO(mat,U,V,dag);};
      virtual void MooDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag){mat=zero;}; // Clover can override these
      virtual void MeeDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag){mat=zero;};

      virtual void DhopDeriv  (GaugeField &mat,const FermionField &U,const FermionField &V,int dag)=0;
      virtual void DhopDerivEO(GaugeField &mat,const FermionField &U,const FermionField &V,int dag)=0;
      virtual void DhopDerivOE(GaugeField &mat,const FermionField &U,const FermionField &V,int dag)=0;


      virtual void  Mdiag  (const FermionField &in, FermionField &out) { Mooee(in,out);};   // Same as Mooee applied to both CB's
      virtual void  Mdir   (const FermionField &in, FermionField &out,int dir,int disp)=0;   // case by case Wilson, Clover, Cayley, ContFrac, PartFrac


      virtual void  MomentumSpacePropagator(FermionField &out,const FermionField &in,RealD _m,std::vector<double> twist) { assert(0);};

      virtual void  FreePropagator(const FermionField &in,FermionField &out,RealD mass,std::vector<double> twist) {
	FFT theFFT((GridCartesian *) in._grid);

	FermionField in_k(in._grid);
	FermionField prop_k(in._grid);

	//phase for boundary condition
	ComplexField coor(in._grid);
	ComplexField ph(in._grid);  ph = zero;
	FermionField in_buf(in._grid); in_buf = zero;
	Scalar ci(0.0,1.0);
	assert(twist.size() == Nd);//check that twist is Nd
	for(unsigned int nu = 0; nu < Nd; nu++)
	{
          LatticeCoordinate(coor, nu);
	  ph = ph + twist[nu]*coor*((1./(in._grid->_fdimensions[nu])));
	}
	in_buf = exp(Scalar(-2.0*M_PI)*ci*ph)*in;

	theFFT.FFT_all_dim(in_k,in_buf,FFT::forward);
        this->MomentumSpacePropagator(prop_k,in_k,mass,twist);
	theFFT.FFT_all_dim(out,prop_k,FFT::backward);

	//phase for boundary condition
	out = out * exp(Scalar(2.0*M_PI)*ci*ph);

      };
      virtual void FreePropagator(const FermionField &in,FermionField &out,RealD mass) {
		std::vector<double> twist(Nd,0.0); //default: periodic boundarys in all directions
	        FreePropagator(in,out,mass,twist);
      };

      ///////////////////////////////////////////////
      // Updates gauge field during HMC
      ///////////////////////////////////////////////
      virtual void ImportGauge(const GaugeField & _U)=0;

      //////////////////////////////////////////////////////////////////////
      // Conserved currents, either contract at sink or insert sequentially.
      //////////////////////////////////////////////////////////////////////

      virtual void ContractConservedCurrent(PropagatorField &q_in_1,
                                            PropagatorField &q_in_2,
                                            PropagatorField &q_out,
                                            Current curr_type,
                                            unsigned int mu)=0;
      virtual void SeqConservedCurrent(PropagatorField &q_in, 
                                       PropagatorField &q_out,
                                       Current curr_type,
                                       unsigned int mu,
                                       unsigned int tmin, 
                                       unsigned int tmax,
                                       ComplexField &lattice_cmplx)=0;

      // Only reimplemented in Wilson5D 
      // Default to just a zero correlation function
      virtual void ContractJ5q(FermionField &q_in   ,ComplexField &J5q) { J5q=zero; };
      virtual void ContractJ5q(PropagatorField &q_in,ComplexField &J5q) { J5q=zero; };

      ///////////////////////////////////////////////
      // Physical field import/export
      ///////////////////////////////////////////////
      virtual void Dminus(const FermionField &psi, FermionField &chi)    { chi=psi; }
      virtual void DminusDag(const FermionField &psi, FermionField &chi) { chi=psi; }
      virtual void ImportPhysicalFermionSource(const FermionField &input,FermionField &imported)
      {
	imported = input;
      };
      virtual void ImportUnphysicalFermion(const FermionField &input,FermionField &imported)
      {
	imported=input;
      };
      virtual void ExportPhysicalFermionSolution(const FermionField &solution,FermionField &exported)
      {
	exported=solution;
      };
      virtual void ExportPhysicalFermionSource(const FermionField &solution,FermionField &exported)
      {
	exported=solution;
      };
    };

  }
}

#endif
