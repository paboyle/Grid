/*************************************************************************************
Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_cayley_cg.cc

Copyright (C) 2022

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Fabian Joswig <fabian.joswig@ed.ac.uk>

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


template<class What>
void  TestConserved(What & Dw,
		    LatticeGaugeField &Umu,
		    GridCartesian         * UGrid,	       GridRedBlackCartesian * UrbGrid,
		    GridParallelRNG *RNG4);

  Gamma::Algebra Gmu [] = {
    Gamma::Algebra::GammaX,
    Gamma::Algebra::GammaY,
    Gamma::Algebra::GammaZ,
    Gamma::Algebra::GammaT,
    Gamma::Algebra::Gamma5
  };

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
								   GridDefaultSimd(Nd,vComplex::Nsimd()),
								   GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG4(UGrid);
  std::vector<int> seeds4({1,2,3,4}); RNG4.SeedFixedIntegers(seeds4);

  LatticeGaugeField Umu(UGrid);
  if( argc > 1 && argv[1][0] != '-' )
  {
    std::cout<<GridLogMessage <<"Loading configuration from "<<argv[1]<<std::endl;
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, argv[1]);
  }
  else
  {
    std::cout<<GridLogMessage <<"Using hot configuration"<<std::endl;
    SU<Nc>::HotConfiguration(RNG4,Umu);
  }

  typename WilsonCloverFermionR::ImplParams params;
  WilsonAnisotropyCoefficients anis;
  RealD mass = 0.1;
  RealD csw_r = 1.0;
  RealD csw_t = 1.0;

  std::cout<<GridLogMessage <<"=================================="<<std::endl;
  std::cout<<GridLogMessage <<"WilsonFermion test"<<std::endl;
  std::cout<<GridLogMessage <<"=================================="<<std::endl;
  WilsonFermionR Dw(Umu,*UGrid,*UrbGrid,mass,params);
  TestConserved<WilsonFermionR>(Dw,Umu,UGrid,UrbGrid,&RNG4);

  std::cout<<GridLogMessage <<"=================================="<<std::endl;
  std::cout<<GridLogMessage <<"WilsonCloverFermion test"<<std::endl;
  std::cout<<GridLogMessage <<"=================================="<<std::endl;
  WilsonCloverFermionR Dwc(Umu, *UGrid, *UrbGrid, mass, csw_r, csw_t, anis, params);
  TestConserved<WilsonCloverFermionR>(Dwc,Umu,UGrid,UrbGrid,&RNG4);

  std::cout<<GridLogMessage <<"=================================="<<std::endl;
  std::cout<<GridLogMessage <<"CompactWilsonCloverFermion test"<<std::endl;
  std::cout<<GridLogMessage <<"=================================="<<std::endl;
  CompactWilsonCloverFermionR Dwcc(Umu, *UGrid, *UrbGrid, mass, csw_r, csw_t, 1.0, anis, params);
  TestConserved<CompactWilsonCloverFermionR>(Dwcc,Umu,UGrid,UrbGrid,&RNG4);

  std::cout<<GridLogMessage <<"=================================="<<std::endl;
  std::cout<<GridLogMessage <<"WilsonExpCloverFermion test"<<std::endl;
  std::cout<<GridLogMessage <<"=================================="<<std::endl;
  WilsonExpCloverFermionR Dewc(Umu, *UGrid, *UrbGrid, mass, csw_r, csw_t, anis, params);
  TestConserved<WilsonExpCloverFermionR>(Dewc,Umu,UGrid,UrbGrid,&RNG4);

  std::cout<<GridLogMessage <<"=================================="<<std::endl;
  std::cout<<GridLogMessage <<"CompactWilsonExpCloverFermion test"<<std::endl;
  std::cout<<GridLogMessage <<"=================================="<<std::endl;
  CompactWilsonExpCloverFermionR Dewcc(Umu, *UGrid, *UrbGrid, mass, csw_r, csw_t, 1.0, anis, params);
  TestConserved<CompactWilsonExpCloverFermionR>(Dewcc,Umu,UGrid,UrbGrid,&RNG4);

  Grid_finalize();
}



template<class Action>
void  TestConserved(Action & Dw,
		    LatticeGaugeField &Umu,
		    GridCartesian         * UGrid,	       GridRedBlackCartesian * UrbGrid,
		    GridParallelRNG *RNG4)
{
  LatticePropagator phys_src(UGrid);
  LatticePropagator seqsrc(UGrid);
  LatticePropagator prop4(UGrid);
  LatticePropagator Vector_mu(UGrid);
  LatticeComplex    SV (UGrid);
  LatticeComplex    VV (UGrid);
  LatticePropagator seqprop(UGrid);

  SpinColourMatrix kronecker; kronecker=1.0;
  Coordinate coor({0,0,0,0});
  phys_src=Zero();
  pokeSite(kronecker,phys_src,coor);

  ConjugateGradient<LatticeFermion> CG(1.0e-16,100000);
  SchurRedBlackDiagTwoSolve<LatticeFermion> schur(CG);
  ZeroGuesser<LatticeFermion> zpg;
  for(int s=0;s<Nd;s++){
    for(int c=0;c<Nc;c++){
      LatticeFermion src4  (UGrid);
      PropToFerm<Action>(src4,phys_src,s,c);

      LatticeFermion result4(UGrid); result4=Zero();
      schur(Dw,src4,result4,zpg);
      std::cout<<GridLogMessage<<"spin "<<s<<" color "<<c<<" norm2(sourc4d) "<<norm2(src4)
               <<" norm2(result4d) "<<norm2(result4)<<std::endl;
      FermToProp<Action>(prop4,result4,s,c);
    }
  }

  auto curr = Current::Vector;
  const int mu_J=0;
  const int t_J=0;

  LatticeComplex    ph (UGrid); ph=1.0;

  Dw.SeqConservedCurrent(prop4,
			   seqsrc,
			   phys_src,
			   curr,
			   mu_J,
			   t_J,
			   t_J,// whole lattice
			   ph);

  for(int s=0;s<Nd;s++){
    for(int c=0;c<Nc;c++){

      LatticeFermion src4  (UGrid);
      PropToFerm<Action>(src4,seqsrc,s,c);

      LatticeFermion result4(UGrid); result4=Zero();
      schur(Dw,src4,result4,zpg);

      FermToProp<Action>(seqprop,result4,s,c);
    }
  }

  Gamma g5(Gamma::Algebra::Gamma5);
  Gamma gT(Gamma::Algebra::GammaT);

  std::vector<TComplex> sumSV;
  std::vector<TComplex> sumVV;

  Dw.ContractConservedCurrent(prop4,prop4,Vector_mu,phys_src,Current::Vector,Tdir);

  SV       = trace(Vector_mu);        // Scalar-Vector conserved current
  VV       = trace(gT*Vector_mu);     // (local) Vector-Vector conserved current

  // Spatial sum
  sliceSum(SV,sumSV,Tdir);
  sliceSum(VV,sumVV,Tdir);

  const int Nt{static_cast<int>(sumSV.size())};

  std::cout<<GridLogMessage<<"Vector Ward identity by timeslice (~ 0)"<<std::endl;
  for(int t=0;t<Nt;t++){
    std::cout<<GridLogMessage <<" t "<<t<<" SV "<<real(TensorRemove(sumSV[t]))<<" VV "<<real(TensorRemove(sumVV[t]))<<std::endl;
    assert(abs(real(TensorRemove(sumSV[t]))) < 1e-10);
    assert(abs(real(TensorRemove(sumVV[t]))) < 1e-2);
  }

  ///////////////////////////////
  // 3pt vs 2pt check
  ///////////////////////////////
  {
    Gamma::Algebra        gA = Gamma::Algebra::Identity;
    Gamma                 g(gA);

    LatticePropagator cur(UGrid);
    LatticePropagator tmp(UGrid);
    LatticeComplex c(UGrid);
    SpinColourMatrix qSite;
    peekSite(qSite, seqprop, coor);

    Complex               test_S, test_V, check_S, check_V;

    std::vector<TComplex> check_buf;

    test_S = trace(qSite*g);
    test_V = trace(qSite*g*Gamma::gmu[mu_J]);

    Dw.ContractConservedCurrent(prop4,prop4,cur,phys_src,curr,mu_J);

    c = trace(cur*g);
    sliceSum(c, check_buf, Tp);
    check_S = TensorRemove(check_buf[t_J]);

    auto gmu=Gamma::gmu[mu_J];
    c = trace(cur*g*gmu);
    sliceSum(c, check_buf, Tp);
    check_V = TensorRemove(check_buf[t_J]);


    std::cout<<GridLogMessage << std::setprecision(14)<<"Test S  = " << abs(test_S)   << std::endl;
    std::cout<<GridLogMessage << "Test V  = " << abs(test_V) << std::endl;
    std::cout<<GridLogMessage << "Check S = " << abs(check_S) << std::endl;
    std::cout<<GridLogMessage << "Check V = " << abs(check_V) << std::endl;

    // Check difference = 0
    check_S = check_S - test_S;
    check_V = check_V - test_V;

    std::cout<<GridLogMessage << "Consistency check for sequential conserved " <<std::endl;
    std::cout<<GridLogMessage << "Diff S  = " << abs(check_S) << std::endl;
    assert(abs(check_S) < 1e-8);
    std::cout<<GridLogMessage << "Diff V  = " << abs(check_V) << std::endl;
    assert(abs(check_V) < 1e-8);
  }

}
