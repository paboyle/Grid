    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_iwasaki_action_newstaple.cc

    Copyright (C) 2015

Author: Christopher Kelly <ckelly@bnl.gov>
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
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

////////////////////////////////////////////////////////////////////////
// PlaqPlusRectangleActoin
////////////////////////////////////////////////////////////////////////
template<class Gimpl>
class PlaqPlusRectangleActionOrig : public Action<typename Gimpl::GaugeField> {
public:

  INHERIT_GIMPL_TYPES(Gimpl);

private:
  RealD c_plaq;
  RealD c_rect;

public:
  PlaqPlusRectangleActionOrig(RealD b,RealD c): c_plaq(b),c_rect(c){};

  virtual std::string action_name(){return "PlaqPlusRectangleActionOrig";}
      
  virtual void refresh(const GaugeField &U, GridSerialRNG &sRNG, GridParallelRNG& pRNG) {}; // noop as no pseudoferms
      
  virtual std::string LogParameters(){
    std::stringstream sstream;
    sstream << GridLogMessage << "["<<action_name() <<"] c_plaq: " << c_plaq << std::endl;
    sstream << GridLogMessage << "["<<action_name() <<"] c_rect: " << c_rect << std::endl;
    return sstream.str();
  }


  virtual RealD S(const GaugeField &U) {
    RealD vol = U.Grid()->gSites();

    RealD plaq = WilsonLoops<Gimpl>::avgPlaquette(U);
    RealD rect = WilsonLoops<Gimpl>::avgRectangle(U);

    RealD action=c_plaq*(1.0 -plaq)*(Nd*(Nd-1.0))*vol*0.5
      +c_rect*(1.0 -rect)*(Nd*(Nd-1.0))*vol;

    return action;
  };

  virtual void deriv(const GaugeField &Umu,GaugeField & dSdU) {
    //extend Ta to include Lorentz indexes
    RealD factor_p = c_plaq/RealD(Nc)*0.5;
    RealD factor_r = c_rect/RealD(Nc)*0.5;

    GridBase *grid = Umu.Grid();

    std::vector<GaugeLinkField> U (Nd,grid);
    std::vector<GaugeLinkField> U2(Nd,grid);

    for(int mu=0;mu<Nd;mu++){
      U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
      WilsonLoops<Gimpl>::RectStapleDouble(U2[mu],U[mu],mu);
    }

    GaugeLinkField dSdU_mu(grid);
    GaugeLinkField staple(grid);

    for (int mu=0; mu < Nd; mu++){

      // Staple in direction mu

      WilsonLoops<Gimpl>::Staple(staple,Umu,mu);

      dSdU_mu = Ta(U[mu]*staple)*factor_p;

      WilsonLoops<Gimpl>::RectStaple(Umu,staple,U2,U,mu);

      dSdU_mu = dSdU_mu + Ta(U[mu]*staple)*factor_r;
	  
      PokeIndex<LorentzIndex>(dSdU, dSdU_mu, mu);
    }

  };

};

// Convenience for common physically defined cases.
//
// RBC c1 parameterisation is not really RBC but don't have good
// reference and we are happy to change name if prior use of this plaq coeff
// parameterisation is made known to us. 
template<class Gimpl>
class RBCGaugeActionOrig : public PlaqPlusRectangleActionOrig<Gimpl> {
public:
  INHERIT_GIMPL_TYPES(Gimpl);
  RBCGaugeActionOrig(RealD beta,RealD c1) : PlaqPlusRectangleActionOrig<Gimpl>(beta*(1.0-8.0*c1), beta*c1) {};
  virtual std::string action_name(){return "RBCGaugeActionOrig";}
};

template<class Gimpl>
class IwasakiGaugeActionOrig : public RBCGaugeActionOrig<Gimpl> {
public:
  INHERIT_GIMPL_TYPES(Gimpl);
  IwasakiGaugeActionOrig(RealD beta) : RBCGaugeActionOrig<Gimpl>(beta,-0.331) {};
  virtual std::string action_name(){return "IwasakiGaugeActionOrig";}
};


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Coordinate latt_size  = GridDefaultLatt();
  Coordinate simd_layout= GridDefaultSimd(Nd,vComplexD::Nsimd());
  Coordinate mpi_layout = GridDefaultMpi();
  std::cout << " mpi "<<mpi_layout<<std::endl;
  std::cout << " simd "<<simd_layout<<std::endl;
  std::cout << " latt "<<latt_size<<std::endl;
  GridCartesian GRID(latt_size,simd_layout,mpi_layout);

  GridParallelRNG   pRNG(&GRID);
  pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
  LatticeGaugeField U(&GRID);

  SU<Nc>::HotConfiguration(pRNG,U);

  //#define PRD
#ifdef PRD
  typedef PeriodicGimplD Gimpl;
#else
  typedef ConjugateGimplD Gimpl;
  std::vector<int> conj_dirs(Nd,0); conj_dirs[0]=1; conj_dirs[3]=1;
  Gimpl::setDirections(conj_dirs);
#endif

  typedef typename WilsonLoops<Gimpl>::GaugeMat GaugeMat;
  typedef typename WilsonLoops<Gimpl>::GaugeLorentz GaugeLorentz;

  GaugeLorentz derivOrig(&GRID), derivNew(&GRID);
  double beta = 2.13;
  IwasakiGaugeActionOrig<Gimpl> action_orig(beta);
  IwasakiGaugeAction<Gimpl> action_new(beta);

  double torig=0, tnew=0;
  int ntest = 10;
  for(int i=0;i<ntest;i++){
    double t0 = usecond();
    action_orig.deriv(U, derivOrig);
    double t1 = usecond();
    action_new.deriv(U, derivNew);
    double t2 = usecond();

    GaugeLorentz diff = derivOrig - derivNew;
    double n = norm2(diff);
    std::cout << GridLogMessage << "Difference " << n << " (expect 0)" << std::endl;
    assert(n<1e-10);

    std::cout << GridLogMessage << "Timings orig: " << (t1-t0)/1000 << "ms,  new: " << (t2-t1)/1000 << "ms" << std::endl;
    torig += (t1-t0)/1000; tnew += (t2-t1)/1000;
  }
  std::cout << GridLogMessage << "Avg timings " << ntest << " iterations: orig:" << torig/ntest << "ms,   new:" << tnew/ntest << "ms" << std::endl;
  
  Grid_finalize();
}
