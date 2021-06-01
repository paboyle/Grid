   /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_ddhmc_matrices.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#include <Grid/qcd/action/momentum/DirichletFilter.h>
#include <Grid/qcd/action/fermion/DirichletFermionOperator.h>
#include <Grid/qcd/action/fermion/SchurFactoredFermionOperator.h>

using namespace std;
using namespace Grid;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  const int Ls=8;
  auto latt = GridDefaultLatt();
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  GridCartesian         * UGridF   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
  GridCartesian         * FGridF   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGridF);
  GridRedBlackCartesian * UrbGridF = SpaceTimeGrid::makeFourDimRedBlackGrid(UGridF);
  GridRedBlackCartesian * FrbGridF = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGridF);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});

  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);

  LatticeFermion src   (FGrid); gaussian(RNG5,src);
  LatticeFermion phi   (FGrid); gaussian(RNG5,phi);
  LatticeFermion chi   (FGrid); gaussian(RNG5,chi);
  LatticeFermion result(FGrid); result=Zero();
  LatticeFermion    ref(FGrid);    ref=Zero();
  LatticeFermion    tmp(FGrid);    tmp=Zero();
  LatticeFermion    tmp1(FGrid);
  LatticeFermion    err(FGrid);    tmp=Zero();
  LatticeGaugeField Umu(UGrid); SU<Nc>::HotConfiguration(RNG4,Umu);
  LatticeGaugeFieldF UmuF(UGridF);
  precisionChange(UmuF,Umu);

  RealD mass=0.1;
  RealD M5  =1.8;

  DomainWallFermionR DdwfPeri(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  DomainWallFermionF DdwfPeriF(UmuF,*FGridF,*FrbGridF,*UGridF,*UrbGridF,mass,M5);

  typedef DomainWallFermionR::Impl_t FimplD;
  typedef DomainWallFermionF::Impl_t FimplF;
  typedef DirichletFermionOperator<FimplD> FermOp;
  typedef DirichletFermionOperator<FimplF> FermOpF;
  Coordinate Block({16,16,16,4});

  DomainWallFermionR DdwfPeriTmp(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  DomainWallFermionF DdwfPeriTmpF(UmuF,*FGridF,*FrbGridF,*UGridF,*UrbGridF,mass,M5);
  FermOp  Ddwf(DdwfPeriTmp,Block); 
  FermOpF DdwfF(DdwfPeriTmpF,Block); 
  Ddwf.ImportGauge(Umu);
  DdwfF.ImportGauge(UmuF);
  
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

  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Test MpcDagMpc is Hermitian              "<<std::endl;
  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
  
  gaussian(RNG5,phi);
  gaussian(RNG5,chi);
  pickCheckerboard(Even,chi_e,chi);
  pickCheckerboard(Odd ,chi_o,chi);
  pickCheckerboard(Even,phi_e,phi);
  pickCheckerboard(Odd ,phi_o,phi);
  RealD t1,t2;


  SchurDiagMooeeOperator<FermOp,LatticeFermion> HermOpEO(Ddwf);
  HermOpEO.MpcDagMpc(chi_e,dchi_e);
  HermOpEO.MpcDagMpc(chi_o,dchi_o);

  HermOpEO.MpcDagMpc(phi_e,dphi_e);
  HermOpEO.MpcDagMpc(phi_o,dphi_o);

  pDce = innerProduct(phi_e,dchi_e);
  pDco = innerProduct(phi_o,dchi_o);
  cDpe = innerProduct(chi_e,dphi_e);
  cDpo = innerProduct(chi_o,dphi_o);

  std::cout<<GridLogMessage <<"e "<<pDce<<" "<<cDpe <<std::endl;
  std::cout<<GridLogMessage <<"o "<<pDco<<" "<<cDpo <<std::endl;

  std::cout<<GridLogMessage <<"pDce - conj(cDpo) "<< pDco-conj(cDpo) <<std::endl;
  std::cout<<GridLogMessage <<"pDco - conj(cDpe) "<< pDce-conj(cDpe) <<std::endl;
  
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing one direction at a time "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  tmp = Zero();
  for(int mu=0;mu<Nd;mu++){

    std::vector<TComplex> slice_ref;
    std::vector<TComplex> slice_result;

    // 5D - Ls is innermost This now annoys me.
    DdwfPeri.Mdir(src,ref   ,mu+1,-1);
    Ddwf.Mdir(src,result,mu+1,-1);

    tmp = tmp + result;
    
    auto lip = localInnerProduct(result,result);
    ref = ref - result;
    auto dip = localInnerProduct(ref,ref);

    sliceSum(lip,slice_result,mu+1);
    sliceSum(dip,slice_ref,mu+1);
    for(int t=0;t<latt[mu];t++){
      std::cout << "mu="<<mu<<" result["<<t<<"] "<<slice_result[t]<<" delta "<<slice_ref[t]<<std::endl;
      //      if( (t%Block[mu]) !=0) assert(norm2(slice_ref[t]) < 1.0e-10);
      //      else assert(norm2(slice_result[t]) == 0.0);
    }

    // Opposite dir
    DdwfPeri.Mdir(src,ref   ,mu+1,1);
    Ddwf.Mdir(src,result,mu+1,1);

    tmp = tmp + result;
    
    lip = localInnerProduct(result,result);
    ref = ref - result;
    dip = localInnerProduct(ref,ref);

    sliceSum(lip,slice_result,mu+1);
    sliceSum(dip,slice_ref,mu+1);
    for(int t=0;t<latt[mu];t++){
      std::cout << "mu="<<mu<<" result["<<t<<"] "<<slice_result[t]<<" delta "<<slice_ref[t]<<std::endl;
      //if( (t%Block[mu]) != Block[mu]-1) assert(norm2(slice_ref[t]) < 1.0e-10);
      //else assert(norm2(slice_result[t]) == 0.0);
    }
    
  }
  pickCheckerboard(Even,src_e,src);
  pickCheckerboard(Odd,src_o,src);
  Ddwf.Meooe(src_e,r_o);  std::cout<<GridLogMessage<<"Applied Meo"<<std::endl;
  Ddwf.Meooe(src_o,r_e);  std::cout<<GridLogMessage<<"Applied Moe"<<std::endl;
  setCheckerboard(r_eo,r_o);
  setCheckerboard(r_eo,r_e);
  ref = r_eo - tmp;
  std::cout << " Difference between Moffdiag and sum over directions is "<<norm2(ref)<<std::endl;
  assert(norm2(ref)<1.0e-10);

  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing that POmega+POmegaBar = 1 "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;

  SchurFactoredFermionOperator<FimplD,FimplF> Schur(DdwfPeri,DdwfPeriF,
						    Ddwf,DdwfF,
						    Block);
  
  result = src;
  Schur.ProjectOmega(result);
  tmp = src;
  Schur.ProjectOmegaBar(tmp);
  std::cout << " norm2(src) "<<norm2(src)<< " "<< norm2(result)<<" "<<norm2(tmp)<<std::endl;
  result = result + tmp - src;
  std::cout << " diff = "<<norm2(result)<<std::endl;
  assert(norm2(result)<=1.0e-8);
 
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing that dBoundary+dBoundaryBar+dOmega+dOmegaBar = Munprec "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  Schur.dBoundary    (src,tmp); result=tmp;        std::cout << "dBoundary    "<<norm2(tmp)<<std::endl;
  Schur.dBoundaryBar (src,tmp); result=result+tmp; std::cout << "dBoundaryBar "<<norm2(tmp)<<std::endl;
  Schur.dOmega       (src,tmp); result=result+tmp; std::cout << "dOmega       "<<norm2(tmp)<<std::endl;
  Schur.dOmegaBar    (src,tmp); result=result+tmp; std::cout << "dOmegaBar    "<<norm2(tmp)<<std::endl;

  DdwfPeri.M(src,ref);  
  err= ref - result;
  std::cout<<GridLogMessage << " norm diff   "<< norm2(err)<< " "<<norm2(ref)<< " " << norm2(result) <<std::endl;
  assert(norm2(err)<=1.0e-8);

  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing that (dBoundary+dBoundaryBar+dOmega+dOmegaBar)dag = Mdag "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  Schur.dBoundaryDag    (src,tmp); result=tmp;        std::cout << "dBoundaryDag    "<<norm2(tmp)<<std::endl;
  Schur.dBoundaryBarDag (src,tmp); result=result+tmp; std::cout << "dBoundaryBarDag "<<norm2(tmp)<<std::endl;
  Schur.dOmegaDag       (src,tmp); result=result+tmp; std::cout << "dOmegaDag       "<<norm2(tmp)<<std::endl;
  Schur.dOmegaBarDag    (src,tmp); result=result+tmp; std::cout << "dOmegaBarDag    "<<norm2(tmp)<<std::endl;

  DdwfPeri.Mdag(src,ref);  
  err= ref - result;
  std::cout<<GridLogMessage << " norm diff   "<< norm2(err)<< " "<<norm2(ref)<< " " << norm2(result) <<std::endl;
  assert(norm2(err)<=1.0e-8);

  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing that <chi|dBoundary|phi> =  <phi|dBoundaryDag|chi>^* "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
 
  Schur.dBoundary(phi,tmp);       std::cout << "<chi|dBoundary|phi>"<<innerProduct(chi,tmp)<<std::endl;
  Schur.dBoundaryDag(chi,tmp);    std::cout << "<phi|dBoundaryDag|chi>"<<innerProduct(phi,tmp)<<std::endl;

  Schur.dBoundaryBar(phi,tmp);    std::cout << "<chi|dBoundaryBar|phi>"<<innerProduct(chi,tmp)<<std::endl;
  Schur.dBoundaryBarDag(chi,tmp); std::cout << "<phi|dBoundaryBarDag|chi>"<<innerProduct(phi,tmp)<<std::endl;

  Schur.dOmega(phi,tmp);       std::cout << "<chi|dOmega|phi>"<<innerProduct(chi,tmp)<<std::endl;
  Schur.dOmegaDag(chi,tmp);    std::cout << "<phi|dOmegaDag|chi>"<<innerProduct(phi,tmp)<<std::endl;

  Schur.dOmegaBar(phi,tmp);    std::cout << "<chi|dOmegaBar|phi>"<<innerProduct(chi,tmp)<<std::endl;
  Schur.dOmegaBarDag(chi,tmp); std::cout << "<phi|dOmegaBarDag|chi>"<<innerProduct(phi,tmp)<<std::endl;
  
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing that dBoundary ProjectBoundary = dBoundary "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;

  tmp = src;
  Schur.ProjectBoundary(tmp);
  Schur.dBoundary(tmp,result);
  Schur.dBoundary(src,tmp);
  result=result - tmp;
  std::cout << " diff = "<<norm2(result)<< " result "<<norm2(tmp)<<" "<<norm2(src)<<std::endl;
  assert(norm2(result)<=1.0e-8);

  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing that dBoundaryBar ProjectBoundaryBar = dBoundaryBar "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;

  tmp = src;
  Schur.ProjectBoundaryBar(tmp);
  Schur.dBoundaryBar(tmp,result);
  Schur.dBoundaryBar(src,tmp);
  result=result - tmp;
  std::cout << " diff = "<<norm2(result)<< " result "<<norm2(tmp)<<std::endl;
  assert(norm2(result)<=1.0e-8);

  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing that dOmega dOmegaInv = 1 "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;

  tmp = src;
  Schur.ProjectOmega(tmp);
  Schur.dOmega(tmp,tmp1);
  Schur.dOmegaInv(tmp1,result);
  tmp=tmp-result;
  std::cout << " diff = "<<norm2(tmp)<< " result "<<norm2(result)<<std::endl;
  assert(norm2(tmp)<=1.0e-8);

  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing that dOmegaBar dOmegaBarInv = 1 "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;

  tmp = src;
  Schur.ProjectOmegaBar(tmp);
  Schur.dOmegaBar(tmp,tmp1);
  Schur.dOmegaBarInv(tmp1,result);

  tmp=tmp-result;
  std::cout << " diff = "<<norm2(tmp)<< " result "<<norm2(result)<<std::endl;
  assert(norm2(tmp)<=1.0e-8);

  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing that dOmegaDag dOmegaDagInv = 1 "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;

  tmp = src;
  Schur.ProjectOmega(tmp);
  Schur.dOmegaDag(tmp,tmp1);
  Schur.dOmegaDagInv(tmp1,result);
  tmp=tmp-result;
  std::cout << " diff = "<<norm2(tmp)<< " result "<<norm2(result)<<std::endl;
  assert(norm2(tmp)<=1.0e-8);

  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing that dOmegaBarDag dOmegaBarDagInv = 1 "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;

  tmp = src;
  Schur.ProjectOmegaBar(tmp);
  Schur.dOmegaBarDag(tmp,tmp1);
  Schur.dOmegaBarDagInv(tmp1,result);
  tmp=tmp-result;
  std::cout << " diff = "<<norm2(tmp)<< " result "<<norm2(result)<<std::endl;
  assert(norm2(tmp)<=1.0e-8);


  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing that R RInv = PboundaryBar "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;

  LatticeFermion Rphi   (FGrid);
  LatticeFermion Rdagchi(FGrid);

  tmp = phi;
  Schur.R(tmp,Rphi);
  Schur.RInv(Rphi,result);

  tmp = phi;
  Schur.ProjectBoundaryBar(tmp);
  //  std::cout << "Project Boundary Bar" << tmp<< std::endl;
  
  tmp=tmp-result;

  std::cout << " diff = "<<norm2(tmp)<< " result "<<norm2(result)<<std::endl;
  assert(norm2(tmp)<1.0e-8);

  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing that Rdag RInvdag = PboundaryBar "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;

  tmp = chi;
  Schur.RDag(tmp,Rdagchi);
  Schur.RDagInv(Rdagchi,result);

  tmp = chi;
  Schur.ProjectBoundaryBar(tmp);
  tmp=tmp-result;

  std::cout << " diff = "<<norm2(tmp)<< " result "<<norm2(result)<<std::endl;
  assert(norm2(tmp)<1.0e-8);
  
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing that <chi|R|phi> = <phi|Rdag|chi>* "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;

  std::cout << "<chi|R|phi>"<<innerProduct(chi,Rphi)<<std::endl;
  std::cout << "<phi|Rdag|chi>"<<innerProduct(phi,Rdagchi)<<std::endl;
  
  Grid_finalize();
}
