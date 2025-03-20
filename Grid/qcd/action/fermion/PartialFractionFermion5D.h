/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/PartialFractionFermion5D.h

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
#ifndef  GRID_QCD_PARTIAL_FRACTION_H
#define  GRID_QCD_PARTIAL_FRACTION_H

#include <Grid/qcd/action/fermion/WilsonFermion5D.h>

NAMESPACE_BEGIN(Grid);

template<class Impl>
class PartialFractionFermion5D : public WilsonFermion5D<Impl>
{
public:
  INHERIT_IMPL_TYPES(Impl);

  const int part_frac_chroma_convention=0;

  void   Meooe_internal(const FermionField &in, FermionField &out,int dag);
  void   Mooee_internal(const FermionField &in, FermionField &out,int dag);
  void   MooeeInv_internal(const FermionField &in, FermionField &out,int dag);
  void   M_internal(const FermionField &in, FermionField &out,int dag);

  // override multiply
  virtual void   M    (const FermionField &in, FermionField &out);
  virtual void   Mdag (const FermionField &in, FermionField &out);

  // half checkerboard operaions
  virtual void   Meooe       (const FermionField &in, FermionField &out);
  virtual void   MeooeDag    (const FermionField &in, FermionField &out);
  virtual void   Mooee       (const FermionField &in, FermionField &out);
  virtual void   MooeeDag    (const FermionField &in, FermionField &out);
  virtual void   MooeeInv    (const FermionField &in, FermionField &out);
  virtual void   MooeeInvDag (const FermionField &in, FermionField &out);

  // force terms; five routines; default to Dhop on diagonal
  virtual void MDeriv  (GaugeField &mat,const FermionField &U,const FermionField &V,int dag);
  virtual void MoeDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag);
  virtual void MeoDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag);

  virtual void   Instantiatable(void) =0; // ensure no make-eee

  // Efficient support for multigrid coarsening
  virtual void  Mdir (const FermionField &in, FermionField &out,int dir,int disp);
  virtual void  MdirAll(const FermionField &in, std::vector<FermionField> &out);

  ///////////////////////////////////////////////////////////////
  // Physical surface field utilities
  ///////////////////////////////////////////////////////////////
  virtual void ExportPhysicalFermionSolution(const FermionField &solution5d,FermionField &exported4d);
  virtual void ImportPhysicalFermionSource  (const FermionField &input4d,FermionField &imported5d);

  // Constructors
  PartialFractionFermion5D(GaugeField &_Umu,
			   GridCartesian         &FiveDimGrid,
			   GridRedBlackCartesian &FiveDimRedBlackGrid,
			   GridCartesian         &FourDimGrid,
			   GridRedBlackCartesian &FourDimRedBlackGrid,
			   RealD _mass,RealD M5,const ImplParams &p= ImplParams());

  PartialFractionFermion5D(GaugeField &_Umu,
			   GridCartesian         &FiveDimGrid,
			   GridRedBlackCartesian &FiveDimRedBlackGrid,
			   GridCartesian         &FourDimGrid,
			   GridRedBlackCartesian &FourDimRedBlackGrid,
			   RealD _mass,RealD M5,std::vector<RealD> &_qmu,const ImplParams &p= ImplParams());

  void FreePropagator(const FermionField &in,FermionField &out,RealD mass,std::vector<Complex> boundary, std::vector<double> twist)
  {
    std::cout << "Free Propagator for PartialFraction"<<std::endl;
    FermionField in_k(in.Grid());
    FermionField prop_k(in.Grid());
    
    FFT theFFT((GridCartesian *) in.Grid());

    //phase for boundary condition
    ComplexField coor(in.Grid());
    ComplexField ph(in.Grid());  ph = Zero();
    FermionField in_buf(in.Grid()); in_buf = Zero();
    typedef typename Simd::scalar_type Scalar;
    Scalar ci(0.0,1.0);
    assert(twist.size() == Nd);//check that twist is Nd
    assert(boundary.size() == Nd);//check that boundary conditions is Nd
    int shift = 0;
    for(unsigned int nu = 0; nu < Nd; nu++)
      {
	// Shift coordinate lattice index by 1 to account for 5th dimension.
	LatticeCoordinate(coor, nu + shift);
	double boundary_phase = ::acos(real(boundary[nu]));
	ph = ph + boundary_phase*coor*((1./(in.Grid()->_fdimensions[nu+shift])));
	//momenta for propagator shifted by twist+boundary
	twist[nu] = twist[nu] + boundary_phase/((2.0*M_PI));
      }
    in_buf = exp(ci*ph*(-1.0))*in;

    theFFT.FFT_all_dim(in_k,in,FFT::forward);
    if ( this->qmu.size() ){
      this->MomentumSpacePropagatorHwQ(prop_k,in_k,mass,twist,this->qmu);
    } else {
      this->MomentumSpacePropagatorHw(prop_k,in_k,mass,twist);
    }
    theFFT.FFT_all_dim(out,prop_k,FFT::backward);
    
    //phase for boundary condition
    out = out * exp(ci*ph);
  };

  virtual void FreePropagator(const FermionField &in,FermionField &out,RealD mass) {
    std::vector<double> twist(Nd,0.0); //default: periodic boundarys in all directions
    std::vector<Complex> boundary;
    for(int i=0;i<Nd;i++) boundary.push_back(1);//default: periodic boundary conditions
    FreePropagator(in,out,mass,boundary,twist);
  };

  void set_qmu(std::vector<RealD> _qmu) { qmu=_qmu; assert(qmu.size()==Nd);};
  void addQmu(const FermionField &in, FermionField &out, int dag);

protected:

  virtual void SetCoefficientsTanh(Approx::zolotarev_data *zdata,RealD scale);
  virtual void SetCoefficientsZolotarev(RealD zolo_hi,Approx::zolotarev_data *zdata);

  std::vector<RealD> qmu;

  // Part frac
  RealD mass;
  RealD dw_diag;
  RealD R;
  RealD amax;
  RealD scale;
  std::vector<double> p; 
  std::vector<double> q;

};

NAMESPACE_END(Grid);

#endif
