    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_gamma.cc

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
#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

//template<class vobj> class is_pod< iScalar<vobj> >
//{
//
//};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
     
  GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

  GridParallelRNG          pRNG(&Grid);
  pRNG.SeedRandomDevice();

  GridSerialRNG            sRNG;
  sRNG.SeedRandomDevice();

  SpinMatrix ident; ident=zero;
  SpinMatrix rnd  ; random(sRNG,rnd);

  SpinMatrix ll; ll=zero;
  SpinMatrix rr; rr=zero;
  SpinMatrix result;

  SpinVector lv; random(sRNG,lv);
  SpinVector rv; random(sRNG,rv);

  //  std::cout<<GridLogMessage << " Is pod " << std::is_pod<SpinVector>::value  << std::endl;
  //  std::cout<<GridLogMessage << " Is pod double   " << std::is_pod<double>::value  << std::endl;
  //  std::cout<<GridLogMessage << " Is pod ComplexF " << std::is_pod<ComplexF>::value  << std::endl;
  //  std::cout<<GridLogMessage << " Is triv double " << std::has_trivial_default_constructor<double>::value  << std::endl;
  //  std::cout<<GridLogMessage << " Is triv ComplexF " << std::has_trivial_default_constructor<ComplexF>::value  << std::endl;
  //  std::cout<<GridLogMessage << " Is pod Scalar<double> " << std::is_pod<iScalar<double> >::value  << std::endl;
  //  std::cout<<GridLogMessage << " Is pod Scalar<ComplexF> " << std::is_pod<iScalar<ComplexF> >::value  << std::endl;
  //  std::cout<<GridLogMessage << " Is pod Scalar<vComplexF> " << std::is_pod<iScalar<vComplexF> >::value  << std::endl;
  //  std::cout<<GridLogMessage << " Is pod Scalar<vComplexD> " << std::is_pod<iScalar<vComplexD> >::value  << std::endl;
  //  std::cout<<GridLogMessage << " Is pod Scalar<vRealF> " << std::is_pod<iScalar<vRealF> >::value  << std::endl;
  //  std::cout<<GridLogMessage << " Is pod Scalar<vRealD> " << std::is_pod<iScalar<vRealD> >::value  << std::endl;
  //  std::cout<<GridLogMessage << " Is triv Scalar<double> " <<std::has_trivial_default_constructor<iScalar<double> >::value << std::endl;
  //  std::cout<<GridLogMessage << " Is triv Scalar<vComplexD> "<<std::has_trivial_default_constructor<iScalar<vComplexD> >::value  << std::endl;

  for(int a=0;a<Ns;a++){
    ident()(a,a) = ComplexF(1.0);
  }

  const Gamma::GammaMatrix *g = Gamma::GammaMatrices;
  const char **list           = Gamma::GammaMatrixNames;

  result =ll*Gamma(g[0])*rr;
  result =ll*Gamma(g[0]);
  rv = Gamma(g[0])*lv;

  for(int mu=0;mu<12;mu++){

    result = Gamma(g[mu])* ident;

    for(int i=0;i<Ns;i++){

      if(i==0) std::cout<<GridLogMessage << list[mu];
      else     std::cout<<GridLogMessage << list[12];

      std::cout<<"(";
      for(int j=0;j<Ns;j++){
	if ( abs(result()(i,j)())==0 ) {
	  std::cout<< " 0";
	} else if ( abs(result()(i,j)() - Complex(0,1))==0){
	  std::cout<< " i";
	} else if ( abs(result()(i,j)() + Complex(0,1))==0){
	  std::cout<< "-i";
	} else if ( abs(result()(i,j)() - Complex(1,0))==0){
	  std::cout<< " 1";
	} else if ( abs(result()(i,j)() + Complex(1,0))==0){
	  std::cout<< "-1";
	}
	std::cout<<((j==Ns-1) ? ")" : "," );
      }
      std::cout << std::endl;
    }

    std::cout << std::endl;

  }

  std::cout << "Testing Gamma^2 - 1 = 0"<<std::endl;
  for(int mu=0;mu<6;mu++){
    result =  Gamma(g[mu])* ident * Gamma(g[mu]);
    result = result - ident;
    RealD mag = norm2(result);
    std::cout << list[mu]<<" " << mag<<std::endl;
  }

  std::cout << "Testing (MinusGamma + G )M = 0"<<std::endl;
  for(int mu=0;mu<6;mu++){
    result =          rnd * Gamma(g[mu]);
    result = result + rnd * Gamma(g[mu+6]);
    RealD mag = norm2(result);
    std::cout << list[mu]<<" " << mag<<std::endl;
  }

  std::cout << "Testing M(MinusGamma + G )  = 0"<<std::endl;
  for(int mu=0;mu<6;mu++){
    result =           Gamma(g[mu])  *rnd;
    result = result +  Gamma(g[mu+6])*rnd;
    RealD mag = norm2(result);
    std::cout << list[mu]<<" " << mag<<std::endl;
  }

  // Testing spins and reconstructs
  SpinVector     recon; random(sRNG,rv);
  SpinVector      full;
  HalfSpinVector hsp,hsm;

  // Xp
  double mag;
  spProjXp(hsm,rv);
  spReconXp(recon,hsm);
  full = rv + Gamma(Gamma::GammaX) *rv;
  mag = TensorRemove(norm2(full-recon));
  std::cout << "Xp "<< mag<<std::endl;

  // Xm
  spProjXm(hsm,rv);
  spReconXm(recon,hsm);
  full = rv - Gamma(Gamma::GammaX) *rv;
  mag = TensorRemove(norm2(full-recon));
  std::cout << "Xm "<< mag<<std::endl;

  // Yp
  spProjYp(hsm,rv);
  spReconYp(recon,hsm);
  full = rv + Gamma(Gamma::GammaY) *rv;
  mag = TensorRemove(norm2(full-recon));
  std::cout << "Yp "<< mag<<std::endl;

  // Ym
  spProjYm(hsm,rv);
  spReconYm(recon,hsm);
  full = rv - Gamma(Gamma::GammaY) *rv;
  mag = TensorRemove(norm2(full-recon));
  std::cout << "Ym "<< mag<<std::endl;

  // Zp
  spProjZp(hsm,rv);
  spReconZp(recon,hsm);
  full = rv + Gamma(Gamma::GammaZ) *rv;
  mag = TensorRemove(norm2(full-recon));
  std::cout << "Zp "<< mag<<std::endl;

  // Zm
  spProjZm(hsm,rv);
  spReconZm(recon,hsm);
  full = rv - Gamma(Gamma::GammaZ) *rv;
  mag = TensorRemove(norm2(full-recon));
  std::cout << "Zm "<< mag<<std::endl;

  // Tp
  spProjTp(hsm,rv);
  spReconTp(recon,hsm);
  full = rv + Gamma(Gamma::GammaT) *rv;
  mag = TensorRemove(norm2(full-recon));
  std::cout << "Tp "<< mag<<std::endl;

  // Tm
  spProjTm(hsm,rv);
  spReconTm(recon,hsm);
  full = rv - Gamma(Gamma::GammaT) *rv;
  mag = TensorRemove(norm2(full-recon));
  std::cout << "Tm "<< mag<<std::endl;

  // 5p
  spProj5p(hsm,rv);
  spRecon5p(recon,hsm);
  full = rv + Gamma(Gamma::Gamma5) *rv;
  mag = TensorRemove(norm2(full-recon));
  std::cout << "5p "<< mag<<std::endl;

  // 5m
  spProj5m(hsm,rv);
  spRecon5m(recon,hsm);
  full = rv - Gamma(Gamma::Gamma5) *rv;
  mag = TensorRemove(norm2(full-recon));
  std::cout << "5m "<< mag<<std::endl;

  Grid_finalize();
}
