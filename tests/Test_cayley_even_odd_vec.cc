    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_cayley_even_odd.cc

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
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

template<class d>
struct scal {
  d internal;
};

    Gamma::Algebra Gmu [] = {
        Gamma::Algebra::GammaX,
        Gamma::Algebra::GammaY,
        Gamma::Algebra::GammaZ,
        Gamma::Algebra::GammaT
    };

typedef DomainWallFermion<DomainWallVec5dImplR>                      DomainWallVecFermionR;
typedef ZMobiusFermion<ZDomainWallVec5dImplR>                        ZMobiusVecFermionR;
typedef MobiusFermion<DomainWallVec5dImplR>                          MobiusVecFermionR;
typedef MobiusZolotarevFermion<DomainWallVec5dImplR>                 MobiusZolotarevVecFermionR;
typedef ScaledShamirFermion<DomainWallVec5dImplR>                    ScaledShamirVecFermionR;
typedef ShamirZolotarevFermion<DomainWallVec5dImplR>                 ShamirZolotarevVecFermionR;
typedef OverlapWilsonCayleyTanhFermion<DomainWallVec5dImplR>         OverlapWilsonCayleyTanhVecFermionR;
typedef OverlapWilsonCayleyZolotarevFermion<DomainWallVec5dImplR>    OverlapWilsonCayleyZolotarevVecFermionR;

template<class What> 
void  TestWhat(What & Ddwf,
	       GridCartesian         * FGrid,	       GridRedBlackCartesian * FrbGrid,
	       GridCartesian         * UGrid,
	       RealD mass, RealD M5,
	       GridParallelRNG *RNG4,   GridParallelRNG *RNG5);

template<class This,class That> 
void  TestMoo(This & Dw, That &sDw);

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  const int Ls=16;

  std::vector<int> latt4  =GridDefaultLatt();

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(latt4, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  GridCartesian         * sUGrid   = SpaceTimeGrid::makeFourDimDWFGrid(latt4,GridDefaultMpi());
  GridRedBlackCartesian * sUrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(sUGrid);
  GridCartesian         * sFGrid   = SpaceTimeGrid::makeFiveDimDWFGrid(Ls,UGrid);
  GridRedBlackCartesian * sFrbGrid = SpaceTimeGrid::makeFiveDimDWFRedBlackGrid(Ls,UGrid);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});

  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          sRNG4(sUGrid);  sRNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          sRNG5(sFGrid);  sRNG5.SeedFixedIntegers(seeds5);

  LatticeGaugeField Umu(UGrid); random(RNG4,Umu);

  RealD mass=0.1;
  RealD M5  =1.8;

  std::cout<<GridLogMessage<<"**************************************************************"<<std::endl;
  std::cout<<GridLogMessage <<"DomainWallFermion vectorised test"<<std::endl;
  std::cout<<GridLogMessage<<"**************************************************************"<<std::endl;
  DomainWallFermionR Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  DomainWallVecFermionR sDdwf(Umu,*sFGrid,*sFrbGrid,*sUGrid,*sUrbGrid,mass,M5);

  TestMoo(Ddwf,sDdwf);
  TestWhat<DomainWallFermionR>(Ddwf,FGrid,FrbGrid,UGrid,mass,M5,&RNG4,&RNG5);
  TestWhat<DomainWallVecFermionR>(sDdwf,sFGrid,sFrbGrid,sUGrid,mass,M5,&sRNG4,&sRNG5);

  RealD b=1.5;// Scale factor b+c=2, b-c=1
  RealD c=0.5;

  std::cout<<GridLogMessage<<"**************************************************************"<<std::endl;
  std::cout<<GridLogMessage <<"MobiusFermion test"<<std::endl;
  std::cout<<GridLogMessage<<"**************************************************************"<<std::endl;
  MobiusFermionR     Dmob(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c);
  MobiusVecFermionR sDmob(Umu,*sFGrid,*sFrbGrid,*sUGrid,*sUrbGrid,mass,M5,b,c);
  TestMoo(Dmob,sDmob);
  TestWhat<MobiusFermionR>(Dmob,FGrid,FrbGrid,UGrid,mass,M5,&RNG4,&RNG5);
  TestWhat<MobiusVecFermionR>(sDmob,sFGrid,sFrbGrid,sUGrid,mass,M5,&sRNG4,&sRNG5);


  std::cout<<GridLogMessage<<"**************************************************************"<<std::endl;
  std::cout<<GridLogMessage <<"Z-MobiusFermion test"<<std::endl;
  std::cout<<GridLogMessage<<"**************************************************************"<<std::endl;
  std::vector<ComplexD> gamma(Ls,std::complex<double>(1.0,0.0));
  ZMobiusFermionR     zDmob(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,gamma,b,c);
  ZMobiusVecFermionR szDmob(Umu,*sFGrid,*sFrbGrid,*sUGrid,*sUrbGrid,mass,M5,gamma,b,c);
  TestMoo(zDmob,szDmob);
  TestWhat<ZMobiusFermionR>(zDmob,FGrid,FrbGrid,UGrid,mass,M5,&RNG4,&RNG5);
  TestWhat<ZMobiusVecFermionR>(szDmob,sFGrid,sFrbGrid,sUGrid,mass,M5,&sRNG4,&sRNG5);

  std::cout<<GridLogMessage<<"**************************************************************"<<std::endl;
  std::cout<<GridLogMessage <<"MobiusZolotarevFermion test"<<std::endl;
  std::cout<<GridLogMessage<<"**************************************************************"<<std::endl;

  MobiusZolotarevFermionR Dzolo(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c,0.1,2.0);
  MobiusZolotarevVecFermionR sDzolo(Umu,*sFGrid,*sFrbGrid,*sUGrid,*sUrbGrid,mass,M5,b,c,0.1,2.0);

  TestMoo(Dzolo,sDzolo);
  TestWhat<MobiusZolotarevFermionR>(Dzolo,FGrid,FrbGrid,UGrid,mass,M5,&RNG4,&RNG5);
  TestWhat<MobiusZolotarevVecFermionR>(sDzolo,sFGrid,sFrbGrid,sUGrid,mass,M5,&sRNG4,&sRNG5);

  std::cout<<GridLogMessage<<"**************************************************************"<<std::endl;
  std::cout<<GridLogMessage <<"ScaledShamirFermion test"<<std::endl;
  std::cout<<GridLogMessage<<"**************************************************************"<<std::endl;

  ScaledShamirFermionR Dsham(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,2.0);
  ScaledShamirVecFermionR sDsham(Umu,*sFGrid,*sFrbGrid,*sUGrid,*sUrbGrid,mass,M5,2.0);

  TestMoo(Dsham,sDsham);
  TestWhat<ScaledShamirFermionR>(Dsham,FGrid,FrbGrid,UGrid,mass,M5,&RNG4,&RNG5);
  TestWhat<ScaledShamirVecFermionR>(sDsham,sFGrid,sFrbGrid,sUGrid,mass,M5,&sRNG4,&sRNG5);

  std::cout<<GridLogMessage<<"**************************************************************"<<std::endl;
  std::cout<<GridLogMessage <<"ShamirZolotarevFermion test"<<std::endl;
  std::cout<<GridLogMessage<<"**************************************************************"<<std::endl;

  ShamirZolotarevFermionR Dshamz(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,0.1,2.0);
  ShamirZolotarevVecFermionR sDshamz(Umu,*sFGrid,*sFrbGrid,*sUGrid,*sUrbGrid,mass,M5,0.1,2.0);

  TestMoo(Dshamz,sDshamz);
  TestWhat<ShamirZolotarevFermionR>(Dshamz,FGrid,FrbGrid,UGrid,mass,M5,&RNG4,&RNG5);
  TestWhat<ShamirZolotarevVecFermionR>(sDshamz,sFGrid,sFrbGrid,sUGrid,mass,M5,&sRNG4,&sRNG5);

  std::cout<<GridLogMessage<<"**************************************************************"<<std::endl;
  std::cout<<GridLogMessage <<"OverlapWilsonCayleyTanhFermion test"<<std::endl;
  std::cout<<GridLogMessage<<"**************************************************************"<<std::endl;
  OverlapWilsonCayleyTanhFermionR Dov(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,1.0);
  OverlapWilsonCayleyTanhVecFermionR sDov(Umu,*sFGrid,*sFrbGrid,*sUGrid,*sUrbGrid,mass,M5,1.0);

  TestMoo(Dov,sDov);
  TestWhat<OverlapWilsonCayleyTanhFermionR>(Dov,FGrid,FrbGrid,UGrid,mass,M5,&RNG4,&RNG5);
  TestWhat<OverlapWilsonCayleyTanhVecFermionR>(sDov,sFGrid,sFrbGrid,sUGrid,mass,M5,&sRNG4,&sRNG5);

  std::cout<<GridLogMessage<<"**************************************************************"<<std::endl;
  std::cout<<GridLogMessage <<"OverlapWilsonCayleyZolotarevFermion test"<<std::endl;
  std::cout<<GridLogMessage<<"**************************************************************"<<std::endl;

  OverlapWilsonCayleyZolotarevFermionR Dovz(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,0.1,2.0);
  OverlapWilsonCayleyZolotarevVecFermionR sDovz(Umu,*sFGrid,*sFrbGrid,*sUGrid,*sUrbGrid,mass,M5,0.1,2.0);

  TestMoo(Dovz,sDovz);
  TestWhat<OverlapWilsonCayleyZolotarevFermionR>(Dovz,FGrid,FrbGrid,UGrid,mass,M5,&RNG4,&RNG5);
  TestWhat<OverlapWilsonCayleyZolotarevVecFermionR>(sDovz,sFGrid,sFrbGrid,sUGrid,mass,M5,&sRNG4,&sRNG5);

  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;

  Grid_finalize();
}

template<class What> 
void  TestWhat(What & Ddwf, 
	       GridCartesian         * FGrid,	       GridRedBlackCartesian * FrbGrid,
	       GridCartesian         * UGrid,
	       RealD mass, RealD M5,
	       GridParallelRNG *RNG4,
	       GridParallelRNG *RNG5)
{

  LatticeFermion src   (FGrid); random(*RNG5,src);
  LatticeFermion phi   (FGrid); random(*RNG5,phi);
  LatticeFermion chi   (FGrid); random(*RNG5,chi);
  LatticeFermion result(FGrid); result=zero;
  LatticeFermion    ref(FGrid);    ref=zero;
  LatticeFermion    tmp(FGrid);    tmp=zero;
  LatticeFermion    err(FGrid);    tmp=zero;

  LatticeFermion src_e (FrbGrid);
  LatticeFermion src_o (FrbGrid);
  LatticeFermion r_e   (FrbGrid);
  LatticeFermion r_o   (FrbGrid);
  LatticeFermion r_eo  (FGrid);
  LatticeFermion r_eeoo(FGrid);

  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing that Meo + Moe + Moo + Mee = Munprec "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;

  pickCheckerboard(Even,src_e,src);
  pickCheckerboard(Odd,src_o,src);

  Ddwf.Meooe(src_e,r_o);  std::cout<<GridLogMessage<<"Applied Meo"<<std::endl;
  Ddwf.Meooe(src_o,r_e);  std::cout<<GridLogMessage<<"Applied Moe"<<std::endl;
  setCheckerboard(r_eo,r_o);
  setCheckerboard(r_eo,r_e);

  Ddwf.Mooee(src_e,r_e);  std::cout<<GridLogMessage<<"Applied Mee"<<std::endl;
  Ddwf.Mooee(src_o,r_o);  std::cout<<GridLogMessage<<"Applied Moo"<<std::endl;
  setCheckerboard(r_eeoo,r_e);
  setCheckerboard(r_eeoo,r_o);

  r_eo=r_eo+r_eeoo;
  Ddwf.M(src,ref);  

  //  std::cout<<GridLogMessage << r_eo<<std::endl;
  //  std::cout<<GridLogMessage << ref <<std::endl;

  err= ref - r_eo;
  std::cout<<GridLogMessage << "EO norm diff   "<< norm2(err)<< " "<<norm2(ref)<< " " << norm2(r_eo) <<std::endl;
    
  LatticeComplex cerr(FGrid);
  cerr = localInnerProduct(err,err);
  //  std::cout<<GridLogMessage << cerr<<std::endl;

  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Test Ddagger is the dagger of D by requiring                "<<std::endl;
  std::cout<<GridLogMessage<<"=  < phi | Deo | chi > * = < chi | Deo^dag| phi>  "<<std::endl;
  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
  
  LatticeFermion chi_e   (FrbGrid);
  LatticeFermion chi_o   (FrbGrid);

  LatticeFermion dchi_e  (FrbGrid);
  LatticeFermion dchi_o  (FrbGrid);

  LatticeFermion phi_e   (FrbGrid);
  LatticeFermion phi_o   (FrbGrid);

  LatticeFermion dphi_e  (FrbGrid);
  LatticeFermion dphi_o  (FrbGrid);


  pickCheckerboard(Even,chi_e,chi);
  pickCheckerboard(Odd ,chi_o,chi);
  pickCheckerboard(Even,phi_e,phi);
  pickCheckerboard(Odd ,phi_o,phi);

  Ddwf.Meooe(chi_e,dchi_o);
  Ddwf.Meooe(chi_o,dchi_e);
  Ddwf.MeooeDag(phi_e,dphi_o);
  Ddwf.MeooeDag(phi_o,dphi_e);

  ComplexD pDce = innerProduct(phi_e,dchi_e);
  ComplexD pDco = innerProduct(phi_o,dchi_o);
  ComplexD cDpe = innerProduct(chi_e,dphi_e);
  ComplexD cDpo = innerProduct(chi_o,dphi_o);

  std::cout<<GridLogMessage <<"e "<<pDce<<" "<<cDpe <<std::endl;
  std::cout<<GridLogMessage <<"o "<<pDco<<" "<<cDpo <<std::endl;

  std::cout<<GridLogMessage <<"pDce - conj(cDpo) "<< pDce-conj(cDpo) <<std::endl;
  std::cout<<GridLogMessage <<"pDco - conj(cDpe) "<< pDco-conj(cDpe) <<std::endl;

  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Test MeeInv Mee = 1                                         "<<std::endl;
  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;

  pickCheckerboard(Even,chi_e,chi);
  pickCheckerboard(Odd ,chi_o,chi);

  Ddwf.Mooee(chi_e,src_e);
  Ddwf.MooeeInv(src_e,phi_e);

  Ddwf.Mooee(chi_o,src_o);
  Ddwf.MooeeInv(src_o,phi_o);
  
  setCheckerboard(phi,phi_e);
  setCheckerboard(phi,phi_o);

  err = phi-chi;
  std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<< std::endl;

  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Test MeeInvDag MeeDag = 1                                   "<<std::endl;
  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;

  pickCheckerboard(Even,chi_e,chi);
  pickCheckerboard(Odd ,chi_o,chi);

  Ddwf.MooeeDag(chi_e,src_e);
  Ddwf.MooeeInvDag(src_e,phi_e);

  Ddwf.MooeeDag(chi_o,src_o);
  Ddwf.MooeeInvDag(src_o,phi_o);
  
  setCheckerboard(phi,phi_e);
  setCheckerboard(phi,phi_o);

  err = phi-chi;
  std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<< std::endl;

  
}


template<class This,class That> 
void  TestMoo(This & Dw, That &sDw)
{
  GridBase *sgrid= sDw.FermionGrid();
  GridBase *ngrid=  Dw.FermionGrid();

  int Ls = Dw.Ls;

  LatticeFermion ssrc(sgrid);
  LatticeFermion nsrc(ngrid);
  LatticeFermion zz(ngrid); zz=zero;
  LatticeFermion sres(sgrid);
  LatticeFermion nres(ngrid);
  LatticeFermion ndiff(ngrid);
  LatticeFermion sdiff(sgrid);

    Gamma g5( Gamma::Algebra::Gamma5 );

  std::vector<int> seeds({1,2,3,4,5,7,8});
  GridParallelRNG    RNG5(ngrid);  
  RNG5.SeedFixedIntegers(seeds);
  random(RNG5,nsrc);
  //  nsrc = nsrc + g5*nsrc;

  //  Lattice<iScalar<vInteger> > coor(ngrid);
  //  LatticeCoordinate(coor,0);//scoor
  //  nsrc = where(coor==(Integer)0,zz,nsrc);

  std::vector<int> latt4(4);
  for(int d=0;d<4;d++){
    latt4[d] = ngrid->_fdimensions[d+1];
  }

  for(int x=0;x<latt4[0];x++){
  for(int y=0;y<latt4[1];y++){
  for(int z=0;z<latt4[2];z++){
  for(int t=0;t<latt4[3];t++){
  for(int s=0;s<Ls;s++){
    std::vector<int> site({s,x,y,z,t});
    SpinColourVector tmp;
    peekSite(tmp,nsrc,site);
    pokeSite(tmp,ssrc,site);
  }}}}}

  sDw.Mooee(ssrc,sres);
   Dw.Mooee(nsrc,nres);

  sDw.MooeeInternal(ssrc,sdiff,DaggerNo,InverseNo);


  for(int x=0;x<latt4[0];x++){
  for(int y=0;y<latt4[1];y++){
  for(int z=0;z<latt4[2];z++){
  for(int t=0;t<latt4[3];t++){
  for(int s=0;s<Ls;s++){

    std::vector<int> site({s,x,y,z,t});
    SpinColourVector stmp;
    SpinColourVector itmp;
    SpinColourVector dtmp;
    peekSite(stmp,sres,site);
    peekSite(itmp,sdiff,site);

    dtmp=itmp-stmp;
    if ( norm2(dtmp)>1.0e-6) { 
      std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
      std:: cout << x<<" "<<y<<" "<< z<< " "<<t<<"; s= "<<s<<std::endl;
      std:: cout << "stmp "<< stmp <<std::endl;
      std:: cout << "itmp "<< itmp <<std::endl;
      std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
    }
  }}}}}
  sdiff = sdiff -sres;
  std::cout<<GridLogMessage<<" norm MooInternal diff "<<norm2(sdiff)<<std::endl;

  for(int x=0;x<latt4[0];x++){
  for(int y=0;y<latt4[1];y++){
  for(int z=0;z<latt4[2];z++){
  for(int t=0;t<latt4[3];t++){
  for(int s=0;s<Ls;s++){
    std::vector<int> site({s,x,y,z,t});

    SpinColourVector tmp;
    peekSite(tmp,sres,site);
    pokeSite(tmp,ndiff,site);
  }}}}}
  ndiff=ndiff-nres;
  std::cout<<GridLogMessage<<" norm Moo diff "<<norm2(ndiff)<<std::endl;

  sDw.MooeeDag(ssrc,sres);
   Dw.MooeeDag(nsrc,nres);
  sDw.MooeeInternal(ssrc,sdiff,DaggerYes,InverseNo);

  sdiff = sdiff -sres;
  std::cout<<GridLogMessage<<" norm MooInternalDag diff "<<norm2(sdiff)<<std::endl;

  for(int x=0;x<latt4[0];x++){
  for(int y=0;y<latt4[1];y++){
  for(int z=0;z<latt4[2];z++){
  for(int t=0;t<latt4[3];t++){
  for(int s=0;s<Ls;s++){
    std::vector<int> site({s,x,y,z,t});

    SpinColourVector tmp;
    peekSite(tmp,sres,site);
    pokeSite(tmp,ndiff,site);
  }}}}}
  ndiff=ndiff-nres;
  std::cout<<GridLogMessage<<" norm MooeeDag diff "<<norm2(ndiff)<<std::endl;

  sDw.MooeeInv(ssrc,sres);
   Dw.MooeeInv(nsrc,nres);

  for(int x=0;x<latt4[0];x++){
  for(int y=0;y<latt4[1];y++){
  for(int z=0;z<latt4[2];z++){
  for(int t=0;t<latt4[3];t++){
  for(int s=0;s<Ls;s++){
    std::vector<int> site({s,x,y,z,t});

    SpinColourVector tmp;
    peekSite(tmp,sres,site);
    pokeSite(tmp,ndiff,site);
  }}}}}
  ndiff=ndiff-nres;
  std::cout<<GridLogMessage<<" norm MooeeInv diff "<<norm2(ndiff)<<std::endl;

  sDw.MooeeInvDag(ssrc,sres);
   Dw.MooeeInvDag(nsrc,nres);

  for(int x=0;x<latt4[0];x++){
  for(int y=0;y<latt4[1];y++){
  for(int z=0;z<latt4[2];z++){
  for(int t=0;t<latt4[3];t++){
  for(int s=0;s<Ls;s++){
    std::vector<int> site({s,x,y,z,t});
    SpinColourVector tmp;
    peekSite(tmp,sres,site);
    pokeSite(tmp,ndiff,site);
  }}}}}
  ndiff=ndiff-nres;
  std::cout<<GridLogMessage<<" norm MooeeInvDag diff "<<norm2(ndiff)<<std::endl;

  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;


}
