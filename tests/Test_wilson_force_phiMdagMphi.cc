    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_wilson_force_phiMdagMphi.cc

    Copyright (C) 2015

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
#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

#define parallel_for PARALLEL_FOR_LOOP for

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();

  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian     RBGrid(latt_size,simd_layout,mpi_layout);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  std::vector<int> seeds({1,2,3,4});

  GridParallelRNG          pRNG(&Grid);
  pRNG.SeedRandomDevice();

  LatticeFermion phi        (&Grid); gaussian(pRNG,phi);
  LatticeFermion Mphi       (&Grid); 
  LatticeFermion Mdagphi       (&Grid); 
  LatticeFermion MphiPrime  (&Grid); 
  LatticeFermion MdagphiPrime  (&Grid); 
  LatticeFermion dMphi      (&Grid); 

  LatticeGaugeField U(&Grid);


  SU3::HotConfiguration(pRNG,U);
  //  SU3::ColdConfiguration(pRNG,U);
  
  ////////////////////////////////////
  // Unmodified matrix element
  ////////////////////////////////////
  RealD mass=-4.0; //kills the diagonal term
  WilsonFermionR Dw     (U,     Grid,RBGrid,mass);
  Dw.M   (phi,Mphi);
  Dw.Mdag(phi,Mdagphi);

  ComplexD S    = innerProduct(Mphi,Mphi); // pdag MdagM p
  ComplexD Sdag = innerProduct(Mdagphi,Mdagphi); // pdag MMdag p

  // get the deriv of phidag MdagM phi with respect to "U"
  LatticeGaugeField UdSdU(&Grid);
  LatticeGaugeField UdSdUdag(&Grid);
  LatticeGaugeField tmp(&Grid);

  Dw.MDeriv(tmp , Mphi,  phi,DaggerNo );  UdSdU=tmp;

  Dw.MDeriv(tmp , Mdagphi,  phi,DaggerYes );  UdSdUdag=tmp;


  LatticeFermion dMdagphi      (&Grid);  dMdagphi=zero;
  LatticeFermion Ftmp      (&Grid);


  //  Dw.MDeriv(UdSdU,Mdagphi,  phi,DaggerYes );// UdSdU =UdSdU +tmp;

  ////////////////////////////////////
  // Modify the gauge field a little in one dir
  ////////////////////////////////////
  RealD dt = 1.0e-3;

  LatticeColourMatrix mommu(&Grid); 
  LatticeGaugeField mom(&Grid); 
  LatticeGaugeField Uprime(&Grid); 

  for(int mu=0;mu<Nd;mu++){

    SU3::GaussianLieAlgebraMatrix(pRNG, mommu); // Traceless antihermitian momentum; gaussian in lie alg

    //    Dw.DoubleStore(Dw.Umu,Uprime); // update U _and_ Udag
    Dw.DhopDirDisp(phi,Ftmp,mu,mu+4,DaggerYes); 
    dMdagphi=dMdagphi+mommu*Ftmp*dt;
    
    PokeIndex<LorentzIndex>(mom,mommu,mu);

    parallel_for(auto i=mom.begin();i<mom.end();i++){
      Uprime[i](mu) =U[i](mu)+ mom[i](mu)*U[i](mu)*dt;
      Dw.Umu[i](mu) =Uprime[i](mu); // update U but _not_ Udag
    }

  }

  Dw.Mdag(phi,MdagphiPrime);
  Dw.M   (phi,MphiPrime);

  std::cout << GridLogMessage << "deltaMdag phi    "<< norm2(dMdagphi) <<std::endl;
  Ftmp=MdagphiPrime - Mdagphi;
  std::cout << GridLogMessage << "diff Mdag phi    "<< norm2(Ftmp) <<std::endl;
  Ftmp = Ftmp - dMdagphi;
  std::cout << GridLogMessage << "err  Mdag phi    "<< norm2(Ftmp) <<std::endl;
  std::cout << dMdagphi<<std::endl;
  Ftmp=MdagphiPrime - Mdagphi;
  std::cout << Ftmp<<std::endl;


  ComplexD Sprime    = innerProduct(Mphi   ,MphiPrime);
  ComplexD Sprimedag = innerProduct(Mdagphi,MdagphiPrime);

  ComplexD deltaSdag = innerProduct(Mdagphi,dMdagphi);
  std::cout << GridLogMessage << "deltaSdag from inner prod of mom* M[u]     "<<deltaSdag<<std::endl;

  //////////////////////////////////////////////
  // Use derivative to estimate dS
  //////////////////////////////////////////////
  LatticeComplex dS(&Grid); dS = zero;
  LatticeComplex dSdag(&Grid); dSdag = zero;
  parallel_for(auto i=mom.begin();i<mom.end();i++){
    for(int mu=0;mu<Nd;mu++){
      //      dS[i]() = dS[i]()+trace(mom[i](mu) * UdSdU[i](mu) - mom[i](mu)* adj( UdSdU[i](mu)) )*dt;
      dS[i]()    =    dS[i]()+trace(mom[i](mu) * UdSdU[i](mu) )*dt;
      dSdag[i]() = dSdag[i]()+trace(mom[i](mu) * UdSdUdag[i](mu) )*dt;
    }
  }
  Complex dSpred    = sum(dS);
  Complex dSdagpred = sum(dSdag);

  std::cout << GridLogMessage << " S      "<<S<<std::endl;
  std::cout << GridLogMessage << " Sprime "<<Sprime<<std::endl;
  std::cout << GridLogMessage << "dS      "<<Sprime-S<<std::endl;
  std::cout << GridLogMessage << "predict dS    "<< dSpred <<std::endl;
  std::cout << "\n\n"<<std::endl;
  std::cout << GridLogMessage << " Sdag      "<<Sdag<<std::endl;
  std::cout << GridLogMessage << " Sprimedag "<<Sprimedag<<std::endl;
  std::cout << GridLogMessage << "dSdag      "<<Sprimedag-Sdag<<std::endl;
  std::cout << GridLogMessage << "predict dSdag    "<< dSdagpred <<std::endl;

  std::cout<< GridLogMessage << "Done" <<std::endl;
  Grid_finalize();
}
