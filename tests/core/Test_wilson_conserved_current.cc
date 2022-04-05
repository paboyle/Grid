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
#include <Grid/qcd/action/fermion/Reconstruct5Dprop.h>

using namespace std;
using namespace Grid;


template<class What>
void  TestConserved(What & Ddwf,
		    LatticeGaugeField &Umu,
		    GridCartesian         * FGrid,	       GridRedBlackCartesian * FrbGrid,
		    GridCartesian         * UGrid,	       GridRedBlackCartesian * UrbGrid,
		    RealD mass, RealD M5,
		    GridParallelRNG *RNG4,
		    GridParallelRNG *RNG5,
            What *Ddwfrev=nullptr);

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

  const int Ls=10;
  std::vector < ComplexD  > omegas;
  std::vector < ComplexD  > omegasrev(Ls);

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
								   GridDefaultSimd(Nd,vComplex::Nsimd()),
								   GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);


  GridCartesian         * UGridF   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
								    GridDefaultSimd(Nd,vComplexF::Nsimd()),
								    GridDefaultMpi());
  GridRedBlackCartesian * UrbGridF = SpaceTimeGrid::makeFourDimRedBlackGrid(UGridF);
  GridCartesian         * FGridF   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGridF);
  GridRedBlackCartesian * FrbGridF = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGridF);


  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
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

  RealD mass=0.3;
  RealD M5  =1.0;
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  std::cout<<GridLogMessage <<"WilsonFermion test"<<std::endl;
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  DomainWallFermionR Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  TestConserved<DomainWallFermionR>(Ddwf,Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);

  std::cout<<GridLogMessage <<"======================"<<std::endl;
  std::cout<<GridLogMessage <<"WilsonCloverFermion test"<<std::endl;
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  MobiusFermionR Dmob(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c);
  TestConserved<MobiusFermionR>(Dmob,Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);

  Grid_finalize();
}



template<class Action>
void  TestConserved(Action & Ddwf,
		    LatticeGaugeField &Umu,
		    GridCartesian         * FGrid,	       GridRedBlackCartesian * FrbGrid,
		    GridCartesian         * UGrid,	       GridRedBlackCartesian * UrbGrid,
		    RealD mass, RealD M5,
		    GridParallelRNG *RNG4,
		    GridParallelRNG *RNG5,
                    Action * Ddwfrev)
{
  LatticePropagator phys_src(UGrid);
  LatticePropagator seqsrc(FGrid);
  LatticePropagator prop5(FGrid);
  LatticePropagator prop5rev(FGrid);
  LatticePropagator prop4(UGrid);
  LatticePropagator Axial_mu(UGrid);
  LatticePropagator Vector_mu(UGrid);
  LatticeComplex    PA (UGrid);
  LatticeComplex    SV (UGrid);
  LatticeComplex    VV (UGrid);
  LatticeComplex    PJ5q(UGrid);
  LatticeComplex    PP (UGrid);
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

      LatticeFermion src5  (FGrid);
      Ddwf.ImportPhysicalFermionSource(src4,src5);

      LatticeFermion result5(FGrid); result5=Zero();
      schur(Ddwf,src5,result5,zpg);
      std::cout<<GridLogMessage<<"spin "<<s<<" color "<<c<<" norm2(sourc5d) "<<norm2(src5)
               <<" norm2(result5d) "<<norm2(result5)<<std::endl;
      FermToProp<Action>(prop5,result5,s,c);

      LatticeFermion result4(UGrid);
      Ddwf.ExportPhysicalFermionSolution(result5,result4);
      FermToProp<Action>(prop4,result4,s,c);

      if( Ddwfrev ) {
        Ddwfrev->ImportPhysicalFermionSource(src4,src5);
        result5 = Zero();
        schur(*Ddwfrev,src5,result5,zpg);
      }
      FermToProp<Action>(prop5rev,result5,s,c);
    }
  }

  auto curr = Current::Vector;
  const int mu_J=0;
  const int t_J=0;

  LatticeComplex    ph (UGrid); ph=1.0;

  Ddwf.SeqConservedCurrent(prop5,
			   seqsrc,
			   phys_src,
			   curr,
			   mu_J,
			   t_J,
			   t_J,// whole lattice
			   ph);

  for(int s=0;s<Nd;s++){
    for(int c=0;c<Nc;c++){

      LatticeFermion src5  (FGrid);
      PropToFerm<Action>(src5,seqsrc,s,c);

      LatticeFermion result5(FGrid); result5=Zero();
      schur(Ddwf,src5,result5,zpg);

      LatticeFermion result4(UGrid);
      Ddwf.ExportPhysicalFermionSolution(result5,result4);
      FermToProp<Action>(seqprop,result4,s,c);
    }
  }

  Gamma g5(Gamma::Algebra::Gamma5);
  Gamma gT(Gamma::Algebra::GammaT);

  std::vector<TComplex> sumPA;
  std::vector<TComplex> sumSV;
  std::vector<TComplex> sumVV;
  std::vector<TComplex> sumPP;
  std::vector<TComplex> sumPJ5q;

  Ddwf.ContractConservedCurrent(prop5rev,prop5,Axial_mu,phys_src,Current::Axial,Tdir);
  Ddwf.ContractConservedCurrent(prop5rev,prop5,Vector_mu,phys_src,Current::Vector,Tdir);
  Ddwf.ContractJ5q(prop5,PJ5q);

  PA       = trace(g5*Axial_mu);      // Pseudoscalar-Axial conserved current
  SV       = trace(Vector_mu);        // Scalar-Vector conserved current
  VV       = trace(gT*Vector_mu);     // (local) Vector-Vector conserved current
  PP       = trace(adj(prop4)*prop4); // Pseudoscalar density

  // Spatial sum
  sliceSum(PA,sumPA,Tdir);
  sliceSum(SV,sumSV,Tdir);
  sliceSum(VV,sumVV,Tdir);
  sliceSum(PP,sumPP,Tdir);
  sliceSum(PJ5q,sumPJ5q,Tdir);

  const int Nt{static_cast<int>(sumPA.size())};
  std::cout<<GridLogMessage<<"Vector Ward identity by timeslice (~ 0)"<<std::endl;
  for(int t=0;t<Nt;t++){
    std::cout<<GridLogMessage <<" t "<<t<<" SV "<<real(TensorRemove(sumSV[t]))<<" VV "<<real(TensorRemove(sumVV[t]))<<std::endl;
  }
  std::cout<<GridLogMessage<<"Axial Ward identity by timeslice (defect ~ 0)"<<std::endl;
  for(int t=0;t<Nt;t++){
    const RealD DmuPAmu{real(TensorRemove(sumPA[t]-sumPA[(t-1+Nt)%Nt]))};
    std::cout<<GridLogMessage<<" t "<<t<<" DmuPAmu "<<DmuPAmu
             <<" PP "<<real(TensorRemove(sumPP[t]))<<" PJ5q "<<real(TensorRemove(sumPJ5q[t]))
             <<" Ward Identity defect " <<(DmuPAmu - 2.*real(TensorRemove(Ddwf.mass*sumPP[t] + sumPJ5q[t])))<<std::endl;
  }

  ///////////////////////////////
  // 3pt vs 2pt check
  ///////////////////////////////
  {
    Gamma::Algebra        gA = (curr == Current::Axial) ? Gamma::Algebra::Gamma5 : Gamma::Algebra::Identity;
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

    Ddwf.ContractConservedCurrent(prop5rev,prop5,cur,phys_src,curr,mu_J);

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
    std::cout<<GridLogMessage << "Diff V  = " << abs(check_V) << std::endl;
  }

}
