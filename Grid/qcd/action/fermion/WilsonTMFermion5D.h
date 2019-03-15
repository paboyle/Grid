    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/WilsonTMFermion5D.h

    Copyright (C) 2015

Author: paboyle <paboyle@ph.ed.ac.uk> ; NB Christoph did similar in GPT

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

#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/WilsonFermion.h>


namespace Grid {

  namespace QCD {
    
    template<class Impl>
      class WilsonTMFermion5D : public WilsonFermion5D<Impl>
      {
      public:
	INHERIT_IMPL_TYPES(Impl);
      public:

	virtual void   Instantiatable(void) {};

	// Constructors
        WilsonTMFermion5D(GaugeField &_Umu,
			  GridCartesian         &Fgrid,
			  GridRedBlackCartesian &Frbgrid, 
			  GridCartesian         &Ugrid,
			  GridRedBlackCartesian &Urbgrid, 
			  const std::vector<RealD> _mass,
			  const std::vector<RealD> _mu,
			  const ImplParams &p= ImplParams()
			  ) :
	WilsonFermion5D<Impl>(_Umu,
			      Fgrid,
			      Frbgrid,
			      Ugrid,
			      Urbgrid,
			      4.0,p)
	
	  {
	    update(_mass,_mu);
	  }

	virtual void Meooe(const FermionField &in, FermionField &out) {
	  if (in.checkerboard == Odd) {
	    this->DhopEO(in, out, DaggerNo);
	  } else {
	    this->DhopOE(in, out, DaggerNo);
	  }
	}

	virtual void MeooeDag(const FermionField &in, FermionField &out) {
	  if (in.checkerboard == Odd) {
	    this->DhopEO(in, out, DaggerYes);
	  } else {
	    this->DhopOE(in, out, DaggerYes);
	  }
	}	
	
	// allow override for twisted mass and clover
	virtual void Mooee(const FermionField &in, FermionField &out) {
	  out.checkerboard = in.checkerboard;
	  //axpibg5x(out,in,a,b); // out = a*in + b*i*G5*in
	  for (int s=0;s<(int)this->mass.size();s++) {
	    ComplexD a = 4.0+this->mass[s];
	    ComplexD b(0.0,this->mu[s]);
	    axpbg5y_ssp(out,a,in,b,in,s,s);
	  }
	}

	virtual void MooeeDag(const FermionField &in, FermionField &out) {
	  out.checkerboard = in.checkerboard;
	  for (int s=0;s<(int)this->mass.size();s++) {
	    ComplexD a = 4.0+this->mass[s];
	    ComplexD b(0.0,-this->mu[s]);
	    axpbg5y_ssp(out,a,in,b,in,s,s);
	  }
	}
	virtual void MooeeInv(const FermionField &in, FermionField &out) {
	  for (int s=0;s<(int)this->mass.size();s++) {
	    RealD m    = this->mass[s];
	    RealD tm   = this->mu[s];
	    RealD mtil = 4.0+this->mass[s];
	    RealD sq   = mtil*mtil+tm*tm;
	    ComplexD a    = mtil/sq;
	    ComplexD b(0.0, -tm /sq);
	    axpbg5y_ssp(out,a,in,b,in,s,s);
	  }
	}
	virtual void MooeeInvDag(const FermionField &in, FermionField &out) {
	  for (int s=0;s<(int)this->mass.size();s++) {
	    RealD m    = this->mass[s];
	    RealD tm   = this->mu[s];
	    RealD mtil = 4.0+this->mass[s];
	    RealD sq   = mtil*mtil+tm*tm;
	    ComplexD a    = mtil/sq;
	    ComplexD b(0.0,tm /sq);
	    axpbg5y_ssp(out,a,in,b,in,s,s);
	  }
	}

	virtual RealD M(const FermionField &in, FermionField &out) {
	  out.checkerboard = in.checkerboard;
	  this->Dhop(in, out, DaggerNo);
	  FermionField tmp(out._grid);
	  for (int s=0;s<(int)this->mass.size();s++) {
	    ComplexD a = 4.0+this->mass[s];
	    ComplexD b(0.0,this->mu[s]);
	    axpbg5y_ssp(tmp,a,in,b,in,s,s);
	  }
	  return axpy_norm(out, 1.0, tmp, out);
	}
	
	// needed for fast PV
	void update(const std::vector<RealD>& _mass, const std::vector<RealD>& _mu) {
	  assert(_mass.size() == _mu.size());
	  assert(_mass.size() == this->FermionGrid()->_fdimensions[0]);
	  this->mass = _mass;
	  this->mu = _mu;
	}
	
      private:
	std::vector<RealD> mu;
	std::vector<RealD> mass;
	
      };
   
    typedef WilsonTMFermion5D<WilsonImplF> WilsonTMFermion5DF; 
    typedef WilsonTMFermion5D<WilsonImplD> WilsonTMFermion5DD; 

}}
