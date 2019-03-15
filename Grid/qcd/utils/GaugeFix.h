    /*************************************************************************************

    grid` physics library, www.github.com/paboyle/Grid 

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
//#include <Grid/Grid.h>

#ifndef GRID_QCD_GAUGE_FIX_H
#define GRID_QCD_GAUGE_FIX_H
namespace Grid {
namespace QCD {

template <class Gimpl> 
class FourierAcceleratedGaugeFixer  : public Gimpl {
 public:
  INHERIT_GIMPL_TYPES(Gimpl);

  typedef typename Gimpl::GaugeLinkField GaugeMat;
  typedef typename Gimpl::GaugeField GaugeLorentz;

  static void GaugeLinkToLieAlgebraField(const std::vector<GaugeMat> &U,std::vector<GaugeMat> &A) {
    for(int mu=0;mu<Nd;mu++){
      Complex cmi(0.0,-1.0);
      A[mu] = Ta(U[mu]) * cmi;
    }
  }
  static void DmuAmu(const std::vector<GaugeMat> &A,GaugeMat &dmuAmu) {
    dmuAmu=zero;
    for(int mu=0;mu<Nd;mu++){
      dmuAmu = dmuAmu + A[mu] - Cshift(A[mu],mu,-1);
    }
  }  
  static void SteepestDescentGaugeFix(GaugeLorentz &Umu,Real & alpha,int maxiter,Real Omega_tol, Real Phi_tol,bool Fourier=false) {
    GridBase *grid = Umu._grid;
    GaugeMat xform(grid);
    SteepestDescentGaugeFix(Umu,xform,alpha,maxiter,Omega_tol,Phi_tol,Fourier);
  }
  static void SteepestDescentGaugeFix(GaugeLorentz &Umu,GaugeMat &xform,Real & alpha,int maxiter,Real Omega_tol, Real Phi_tol,bool Fourier=false) {

    GridBase *grid = Umu._grid;

    Real org_plaq      =WilsonLoops<Gimpl>::avgPlaquette(Umu);
    Real org_link_trace=WilsonLoops<Gimpl>::linkTrace(Umu); 
    Real old_trace = org_link_trace;
    Real trG;
    
    xform=1.0;

    std::vector<GaugeMat> U(Nd,grid);

    GaugeMat dmuAmu(grid);

    for(int i=0;i<maxiter;i++){

      for(int mu=0;mu<Nd;mu++) U[mu]= PeekIndex<LorentzIndex>(Umu,mu);

      if ( Fourier==false ) { 
	trG = SteepestDescentStep(U,xform,alpha,dmuAmu);
      } else { 
	trG = FourierAccelSteepestDescentStep(U,xform,alpha,dmuAmu);
      }
      for(int mu=0;mu<Nd;mu++) PokeIndex<LorentzIndex>(Umu,U[mu],mu);
      // Monitor progress and convergence test 
      // infrequently to minimise cost overhead
      if ( i %20 == 0 ) { 
	Real plaq      =WilsonLoops<Gimpl>::avgPlaquette(Umu);
	Real link_trace=WilsonLoops<Gimpl>::linkTrace(Umu); 

	if (Fourier) 
	  std::cout << GridLogMessage << "Fourier Iteration "<<i<< " plaq= "<<plaq<< " dmuAmu " << norm2(dmuAmu)<< std::endl;
	else 
	  std::cout << GridLogMessage << " Iteration "<<i<< " plaq= "<<plaq<< " dmuAmu " << norm2(dmuAmu)<< std::endl;
	
	Real Phi  = 1.0 - old_trace / link_trace ;
	Real Omega= 1.0 - trG;

	std::cout << GridLogMessage << " Iteration "<<i<< " Phi= "<<Phi<< " Omega= " << Omega<< " trG " << trG <<std::endl;
	if ( (Omega < Omega_tol) && ( ::fabs(Phi) < Phi_tol) ) {
	  std::cout << GridLogMessage << "Converged ! "<<std::endl;
	  return;
	}

	old_trace = link_trace;

      }
    }
  };
  static Real SteepestDescentStep(std::vector<GaugeMat> &U,GaugeMat &xform,Real & alpha, GaugeMat & dmuAmu) {
    GridBase *grid = U[0]._grid;

    std::vector<GaugeMat> A(Nd,grid);
    GaugeMat g(grid);

    GaugeLinkToLieAlgebraField(U,A);
    ExpiAlphaDmuAmu(A,g,alpha,dmuAmu);


    Real vol = grid->gSites();
    Real trG = TensorRemove(sum(trace(g))).real()/vol/Nc;

    xform = g*xform ;
    SU<Nc>::GaugeTransform(U,g);

    return trG;
  }

  static Real FourierAccelSteepestDescentStep(std::vector<GaugeMat> &U,GaugeMat &xform,Real & alpha, GaugeMat & dmuAmu) {

    GridBase *grid = U[0]._grid;

    Real vol = grid->gSites();

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

      Real TwoPiL =  M_PI * 2.0/ latt_size[mu];
      LatticeCoordinate(pmu,mu);
      pmu = TwoPiL * pmu ;
      psq = psq + 4.0*sin(pmu*0.5)*sin(pmu*0.5); 
    }

    Complex psqMax(16.0);
    Fp =  psqMax*one/psq;

    pokeSite(TComplex(1.0),Fp,coor);

    dmuAmu_p  = dmuAmu_p * Fp; 

    theFFT.FFT_all_dim(dmuAmu,dmuAmu_p,FFT::backward);

    GaugeMat ciadmam(grid);
    Complex cialpha(0.0,-alpha);
    ciadmam = dmuAmu*cialpha;
    SU<Nc>::taExp(ciadmam,g);

    Real trG = TensorRemove(sum(trace(g))).real()/vol/Nc;

    xform = g*xform ;
    SU<Nc>::GaugeTransform(U,g);

    return trG;
  }

  static void ExpiAlphaDmuAmu(const std::vector<GaugeMat> &A,GaugeMat &g,Real & alpha, GaugeMat &dmuAmu) {
    GridBase *grid = g._grid;
    Complex cialpha(0.0,-alpha);
    GaugeMat ciadmam(grid);
    DmuAmu(A,dmuAmu);
    ciadmam = dmuAmu*cialpha;
    SU<Nc>::taExp(ciadmam,g);
  }  
};

}
}
#endif
