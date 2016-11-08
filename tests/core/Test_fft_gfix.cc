    /*************************************************************************************

    grid` physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_cshift.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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

using namespace Grid;
using namespace Grid::QCD;

template <class Gimpl> 
class FourierAcceleratedGaugeFixer  : public Gimpl {
  public:
  INHERIT_GIMPL_TYPES(Gimpl);

  typedef typename Gimpl::GaugeLinkField GaugeMat;
  typedef typename Gimpl::GaugeField GaugeLorentz;

  static void GaugeLinkToLieAlgebraField(const std::vector<GaugeMat> &U,std::vector<GaugeMat> &A) {
    for(int mu=0;mu<Nd;mu++){
//      ImplComplex cmi(0.0,-1.0);
      ComplexD cmi(0.0,-1.0);
      A[mu] = Ta(U[mu]) * cmi;
    }
  }
  static void DmuAmu(const std::vector<GaugeMat> &A,GaugeMat &dmuAmu) {
    dmuAmu=zero;
    for(int mu=0;mu<Nd;mu++){
      dmuAmu = dmuAmu + A[mu] - Cshift(A[mu],mu,-1);
    }
  }  
  static void SteepestDescentGaugeFix(GaugeLorentz &Umu,RealD & alpha,int maxiter,RealD Omega_tol, RealD Phi_tol) {
    GridBase *grid = Umu._grid;

    RealD org_plaq      =WilsonLoops<Gimpl>::avgPlaquette(Umu);
    RealD org_link_trace=WilsonLoops<Gimpl>::linkTrace(Umu); 
    RealD old_trace = org_link_trace;
    RealD trG;

    std::vector<GaugeMat> U(Nd,grid);
                 GaugeMat dmuAmu(grid);

    for(int i=0;i<maxiter;i++){
      for(int mu=0;mu<Nd;mu++) U[mu]= PeekIndex<LorentzIndex>(Umu,mu);
      //trG = SteepestDescentStep(U,alpha,dmuAmu);
      trG = FourierAccelSteepestDescentStep(U,alpha,dmuAmu);
      for(int mu=0;mu<Nd;mu++) PokeIndex<LorentzIndex>(Umu,U[mu],mu);
      // Monitor progress and convergence test 
      // infrequently to minimise cost overhead
      if ( i %20 == 0 ) { 
	RealD plaq      =WilsonLoops<Gimpl>::avgPlaquette(Umu);
	RealD link_trace=WilsonLoops<Gimpl>::linkTrace(Umu); 

	std::cout << GridLogMessage << " Iteration "<<i<< " plaq= "<<plaq<< " dmuAmu " << norm2(dmuAmu)<< std::endl;
	
	RealD Phi  = 1.0 - old_trace / link_trace ;
	RealD Omega= 1.0 - trG;


	std::cout << GridLogMessage << " Iteration "<<i<< " Phi= "<<Phi<< " Omega= " << Omega<< " trG " << trG <<std::endl;
	if ( (Omega < Omega_tol) && ( ::fabs(Phi) < Phi_tol) ) {
	  std::cout << GridLogMessage << "Converged ! "<<std::endl;
	  return;
	}

	old_trace = link_trace;

      }
    }
  };
  static RealD SteepestDescentStep(std::vector<GaugeMat> &U,RealD & alpha, GaugeMat & dmuAmu) {
    GridBase *grid = U[0]._grid;

    std::vector<GaugeMat> A(Nd,grid);
    GaugeMat g(grid);

    GaugeLinkToLieAlgebraField(U,A);
    ExpiAlphaDmuAmu(A,g,alpha,dmuAmu);


    RealD vol = grid->gSites();
    RealD trG = TensorRemove(sum(trace(g))).real()/vol/Nc;

    SU<Nc>::GaugeTransform(U,g);

    return trG;
  }

  static RealD FourierAccelSteepestDescentStep(std::vector<GaugeMat> &U,RealD & alpha, GaugeMat & dmuAmu) {

    GridBase *grid = U[0]._grid;

    RealD vol = grid->gSites();

    FFT theFFT((GridCartesian *)grid);

    LatticeComplex  Fp(grid);
    LatticeComplex  psq(grid); psq=zero;
    LatticeComplex  pmu(grid); 
    LatticeComplex   one(grid); one = Complex(1.0,0.0);

    GaugeMat g(grid);
    GaugeMat dmuAmu_p(grid);
    std::vector<GaugeMat> A(Nd,grid);

    GaugeLinkToLieAlgebraField(U,A);

    DmuAmu(A,dmuAmu);

    theFFT.FFT_all_dim(dmuAmu_p,dmuAmu,FFT::forward);

    //////////////////////////////////
    // Work out Fp = psq_max/ psq...
    //////////////////////////////////
    std::vector<int> latt_size = grid->GlobalDimensions();
    std::vector<int> coor(grid->_ndimension,0);
    for(int mu=0;mu<Nd;mu++) {

      RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];
      LatticeCoordinate(pmu,mu);
      pmu = TwoPiL * pmu ;
      psq = psq + 4.0*sin(pmu*0.5)*sin(pmu*0.5); 
    }

    ComplexD psqMax(16.0);
    Fp =  psqMax*one/psq;

    static int once;
    if ( once == 0 ) { 
      std::cout << " Fp " << Fp <<std::endl;
      once ++;
    }
    pokeSite(TComplex(1.0),Fp,coor);

    dmuAmu_p  = dmuAmu_p * Fp; 

    theFFT.FFT_all_dim(dmuAmu,dmuAmu_p,FFT::backward);

    GaugeMat ciadmam(grid);
    ComplexD cialpha(0.0,-alpha);
    ciadmam = dmuAmu*cialpha;
    SU<Nc>::taExp(ciadmam,g);

    RealD trG = TensorRemove(sum(trace(g))).real()/vol/Nc;

    SU<Nc>::GaugeTransform(U,g);

    return trG;
  }

  static void ExpiAlphaDmuAmu(const std::vector<GaugeMat> &A,GaugeMat &g,RealD & alpha, GaugeMat &dmuAmu) {
    GridBase *grid = g._grid;
    ComplexD cialpha(0.0,-alpha);
    GaugeMat ciadmam(grid);
    DmuAmu(A,dmuAmu);
    ciadmam = dmuAmu*cialpha;
    SU<Nc>::taExp(ciadmam,g);
  }  
/*
  ////////////////////////////////////////////////////////////////
  // NB The FT for fields living on links has an extra phase in it
  // Could add these to the FFT class as a later task since this code
  // might be reused elsewhere ????
  ////////////////////////////////////////////////////////////////
  static void InverseFourierTransformAmu(FFT &theFFT,const std::vector<GaugeMat> &Ap,std::vector<GaugeMat> &Ax) {
    GridBase * grid = theFFT.Grid();
    std::vector<int> latt_size = grid->GlobalDimensions();

    ComplexField  pmu(grid);
    ComplexField  pha(grid);
    GaugeMat      Apha(grid);

    ComplexD ci(0.0,1.0);

    for(int mu=0;mu<Nd;mu++){

      RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];
      LatticeCoordinate(pmu,mu);
      pmu = TwoPiL * pmu ;
      pha = exp(pmu *  (0.5 *ci)); // e(ipmu/2) since Amu(x+mu/2)

      Apha = Ap[mu] * pha;

      theFFT.FFT_all_dim(Apha,Ax[mu],FFT::backward);
    }
  }
  static void FourierTransformAmu(FFT & theFFT,const std::vector<GaugeMat> &Ax,std::vector<GaugeMat> &Ap) {
    GridBase * grid = theFFT.Grid();
    std::vector<int> latt_size = grid->GlobalDimensions();

    ComplexField  pmu(grid);
    ComplexField  pha(grid);
    ComplexD ci(0.0,1.0);
    
    // Sign convention for FFTW calls:
    // A(x)= Sum_p e^ipx A(p) / V
    // A(p)= Sum_p e^-ipx A(x)

    for(int mu=0;mu<Nd;mu++){
      RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];
      LatticeCoordinate(pmu,mu);
      pmu = TwoPiL * pmu ;
      pha = exp(-pmu *  (0.5 *ci)); // e(+ipmu/2) since Amu(x+mu/2)

      theFFT.FFT_all_dim(Ax[mu],Ap[mu],FFT::backward);
      Ap[mu] = Ap[mu] * pha;
    }
  }
*/
};

int main (int argc, char ** argv)
{
  std::vector<int> seeds({1,2,3,4});

  Grid_init(&argc,&argv);

  int threads = GridThread::GetThreads();

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout( { vComplexD::Nsimd(),1,1,1});
  std::vector<int> mpi_layout  = GridDefaultMpi();

  int vol = 1;
  for(int d=0;d<latt_size.size();d++){
    vol = vol * latt_size[d];
  }

  GridCartesian         GRID(latt_size,simd_layout,mpi_layout);
  GridSerialRNG          sRNG;  sRNG.SeedFixedIntegers(seeds); // naughty seeding
  GridParallelRNG          pRNG(&GRID);   pRNG.SeedFixedIntegers(seeds);

  FFT theFFT(&GRID);

  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  std::cout<< "*****************************************************************" <<std::endl;
  std::cout<< "* Testing we can gauge fix steep descent a RGT of Unit gauge    *" <<std::endl;
  std::cout<< "*****************************************************************" <<std::endl;

  LatticeGaugeField   Umu(&GRID);
  LatticeGaugeField   Uorg(&GRID);
  LatticeColourMatrix   g(&GRID); // Gauge xform

  
  SU3::ColdConfiguration(pRNG,Umu); // Unit gauge
  Uorg=Umu;

  SU3::RandomGaugeTransform(pRNG,Umu,g); // Unit gauge
  RealD plaq=WilsonLoops<PeriodicGimplF>::avgPlaquette(Umu);
  std::cout << " Initial plaquette "<<plaq << std::endl;



  RealD alpha=0.1;
  FourierAcceleratedGaugeFixer<PeriodicGimplF>::SteepestDescentGaugeFix(Umu,alpha,10000,1.0e-10, 1.0e-10);


  plaq=WilsonLoops<PeriodicGimplF>::avgPlaquette(Umu);
  std::cout << " Final plaquette "<<plaq << std::endl;

  Uorg = Uorg - Umu;
  std::cout << " Norm Difference "<< norm2(Uorg) << std::endl;


  //  std::cout<< "*****************************************************************" <<std::endl;
  //  std::cout<< "* Testing Fourier accelerated fixing                            *" <<std::endl;
  //  std::cout<< "*****************************************************************" <<std::endl;

  //  std::cout<< "*****************************************************************" <<std::endl;
  //  std::cout<< "* Testing non-unit configuration                                *" <<std::endl;
  //  std::cout<< "*****************************************************************" <<std::endl;



  Grid_finalize();
}
