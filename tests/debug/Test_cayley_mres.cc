/*************************************************************************************
    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_cayley_cg.cc

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
#include <Grid/qcd/action/fermion/Reconstruct5Dprop.h>

using namespace std;
using namespace Grid;


template<class What> 
void  TestConserved(What & Ddwf, What & Ddwfrev, 
		    LatticeGaugeField &Umu,
		    GridCartesian         * FGrid,	       GridRedBlackCartesian * FrbGrid,
		    GridCartesian         * UGrid,	       GridRedBlackCartesian * UrbGrid,
		    RealD mass, RealD M5,
		    GridParallelRNG *RNG4,
		    GridParallelRNG *RNG5);

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

#if 1
  omegas.push_back( std::complex<double>(1.45806438985048,-0) );
  omegas.push_back( std::complex<double>(0.830951166685955,-0) );
  omegas.push_back( std::complex<double>(0.341985020453729,-0) );
  omegas.push_back( std::complex<double>(0.126074299502912,-0) );
  //  omegas.push_back( std::complex<double>(0.0686324988446592,0.0550658530827402) );
  //  omegas.push_back( std::complex<double>(0.0686324988446592,-0.0550658530827402) );
  omegas.push_back( std::complex<double>(0.0686324988446592,0));
  omegas.push_back( std::complex<double>(0.0686324988446592,0));
  omegas.push_back( std::complex<double>(0.0990136651962626,-0) );
  omegas.push_back( std::complex<double>(0.21137902619029,-0) );
  omegas.push_back( std::complex<double>(0.542352409156791,-0) );
  omegas.push_back( std::complex<double>(1.18231318389348,-0) );
#else 
  omegas.push_back( std::complex<double>(0.8,0.0));
  omegas.push_back( std::complex<double>(1.1,0.0));
  omegas.push_back( std::complex<double>(1.2,0.0));
  omegas.push_back( std::complex<double>(1.3,0.0));
  omegas.push_back( std::complex<double>(0.5,0.2));
  omegas.push_back( std::complex<double>(0.5,-0.2));
  omegas.push_back( std::complex<double>(0.8,0.0));
  omegas.push_back( std::complex<double>(1.1,0.0));
  omegas.push_back( std::complex<double>(1.2,0.0));
  omegas.push_back( std::complex<double>(1.3,0.0));
#endif

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


  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeGaugeField Umu(UGrid);
  SU3::ColdConfiguration(Umu);
  //  SU3::HotConfiguration(RNG4,Umu);

  RealD mass=0.3;
  RealD M5  =1.0;
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  std::cout<<GridLogMessage <<"DomainWallFermion test"<<std::endl;
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  DomainWallFermionR Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  TestConserved<DomainWallFermionR>(Ddwf,Ddwf,Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);

  RealD b=1.5;// Scale factor b+c=2, b-c=1
  RealD c=0.5;
  //  std::vector<ComplexD> gamma(Ls,ComplexD(1.0,0.0));

  std::cout<<GridLogMessage <<"======================"<<std::endl;
  std::cout<<GridLogMessage <<"MobiusFermion test"<<std::endl;
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  MobiusFermionR Dmob(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c);
  TestConserved<MobiusFermionR>(Dmob,Dmob,Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);

  std::cout<<GridLogMessage <<"======================"<<std::endl;
  std::cout<<GridLogMessage <<"ScaledShamirFermion test"<<std::endl;
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  ScaledShamirFermionR Dsham(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,2.0);
  TestConserved<ScaledShamirFermionR>(Dsham,Dsham,Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);

  std::cout<<GridLogMessage <<"======================"<<std::endl;
  std::cout<<GridLogMessage <<"ZMobiusFermion test"<<std::endl;
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  for(int s=0;s<Ls;s++) omegasrev[s]=conjugate(omegas[Ls-1-s]);
  //  for(int s=0;s<Ls;s++) omegasrev[s]=omegas[Ls-1-s];
  ZMobiusFermionR ZDmob(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,omegas,b,c);
  ZMobiusFermionR ZDmobrev(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,omegasrev,b,c);

  TestConserved<ZMobiusFermionR>(ZDmob,ZDmobrev,Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);

  Grid_finalize();
}



template<class Action> 
void  TestConserved(Action & Ddwf, 
		    Action & Ddwfrev, 
		    LatticeGaugeField &Umu,
		    GridCartesian         * FGrid,	       GridRedBlackCartesian * FrbGrid,
		    GridCartesian         * UGrid,	       GridRedBlackCartesian * UrbGrid,
		    RealD mass, RealD M5,
		    GridParallelRNG *RNG4,
		    GridParallelRNG *RNG5)
{
  int Ls=Ddwf.Ls;

  LatticePropagator phys_src(UGrid); 

  std::vector<LatticeColourMatrix> U(4,UGrid);
  
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
  
  MdagMLinearOperator<Action,LatticeFermion> HermOp(Ddwf);
  MdagMLinearOperator<Action,LatticeFermion> HermOprev(Ddwfrev);
  ConjugateGradient<LatticeFermion> CG(1.0e-16,100000);
  for(int s=0;s<Nd;s++){
    for(int c=0;c<Nc;c++){
      LatticeFermion src4  (UGrid); 
      PropToFerm<Action>(src4,phys_src,s,c);

      LatticeFermion src5  (FGrid); 
      Ddwf.ImportPhysicalFermionSource(src4,src5);

      LatticeFermion result5(FGrid); result5=Zero();

      // CGNE
      LatticeFermion Mdagsrc5  (FGrid); 
      Ddwf.Mdag(src5,Mdagsrc5);
      CG(HermOp,Mdagsrc5,result5);
      FermToProp<Action>(prop5,result5,s,c);

      LatticeFermion result4(UGrid);
      Ddwf.ExportPhysicalFermionSolution(result5,result4);
      FermToProp<Action>(prop4,result4,s,c);

      Ddwfrev.ImportPhysicalFermionSource(src4,src5);
      Ddwfrev.Mdag(src5,Mdagsrc5);
      CG(HermOprev,Mdagsrc5,result5);
      FermToProp<Action>(prop5rev,result5,s,c);
    }
  }

#if 1
  auto curr = Current::Axial;
  const int mu_J=Nd-1;
#else
  auto curr = Current::Vector;
  const int mu_J=0;
#endif
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

      // CGNE
      LatticeFermion Mdagsrc5  (FGrid); 
      Ddwf.Mdag(src5,Mdagsrc5);
      CG(HermOp,Mdagsrc5,result5);

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
  
  PA       = trace(g5*Axial_mu);
  SV       = trace(Vector_mu);
  VV       = trace(gT*Vector_mu);
  PP       = trace(adj(prop4)*prop4);
  
  // Spatial sum
  sliceSum(PA,sumPA,Tdir);
  sliceSum(SV,sumSV,Tdir);
  sliceSum(VV,sumVV,Tdir);
  sliceSum(PP,sumPP,Tdir);
  sliceSum(PJ5q,sumPJ5q,Tdir);

  int Nt=sumPA.size();
  for(int t=0;t<Nt;t++){
    std::cout <<" SV "<<real(TensorRemove(sumSV[t]));
    std::cout <<" VV "<<real(TensorRemove(sumVV[t]))<<std::endl;
  }
  for(int t=0;t<Nt;t++){
    std::cout <<" PAc "<<real(TensorRemove(sumPA[t]));
    std::cout <<" PJ5q "<<real(TensorRemove(sumPJ5q[t]));
    std::cout <<" Ward Identity defect " <<real(TensorRemove(sumPA[t]-sumPA[(t-1+Nt)%Nt] - 2.0*(Ddwf.mass*sumPP[t] + sumPJ5q[t]) ))<<"\n";
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

      /*
#if 0
template<class Action> 
void  TestConserved1(Action & Ddwf, Action & Ddwfrev, 
		       LatticeGaugeField &Umu,
		       GridCartesian         * FGrid,	       GridRedBlackCartesian * FrbGrid,
		       GridCartesian         * UGrid,	       GridRedBlackCartesian * UrbGrid,
		       RealD mass, RealD M5,
		       GridParallelRNG *RNG4,
		       GridParallelRNG *RNG5)
{
  int Ls=Ddwf.Ls;

  LatticePropagator phys_src(UGrid); 

  std::vector<LatticeColourMatrix> U(4,UGrid);
  
  LatticePropagator prop5(FGrid); 
  LatticePropagator prop5rev(FGrid); 
  LatticePropagator prop4(UGrid); 
  LatticePropagator Axial_mu(UGrid); 
  LatticeComplex    PA (UGrid); 
  LatticeComplex    PAxyz(UGrid); 
  LatticeComplex    PJ5q(UGrid);
  LatticeComplex    PP (UGrid);
  std::vector<LatticePropagator> prop(Ls,UGrid);
  std::vector<LatticePropagator> proprev(Ls,UGrid);

  SpinColourMatrix kronecker; kronecker=1.0;
  std::cout << kronecker << std::endl;
  phys_src=Zero();
  pokeSite(kronecker,phys_src,Coordinate({0,0,0,0}));
  
  MdagMLinearOperator<Action,LatticeFermion> HermOp(Ddwf);
  MdagMLinearOperator<Action,LatticeFermion> HermOprev(Ddwfrev);
  ConjugateGradient<LatticeFermion> CG(1.0e-12,10000);

  for(int s=0;s<Nd;s++){
    for(int c=0;c<Nc;c++){
      LatticeFermion src4  (UGrid); 
      PropToFerm<Action>(src4,phys_src,s,c);

      LatticeFermion src5  (FGrid); 
      Ddwf.ImportPhysicalFermionSource(src4,src5);

      LatticeFermion result5(FGrid); result5=Zero();

      // CGNE
      LatticeFermion Mdagsrc5  (FGrid); 
      Ddwf.Mdag(src5,Mdagsrc5);
      CG(HermOp,Mdagsrc5,result5);
      FermToProp<Action>(prop5,result5,s,c);

      LatticeFermion result4(UGrid);
      Ddwf.ExportPhysicalFermionSolution(result5,result4);
      FermToProp<Action>(prop4,result4,s,c);

      Ddwfrev.Mdag(src5,Mdagsrc5);
      CG(HermOprev,Mdagsrc5,result5);
      FermToProp<Action>(prop5rev,result5,s,c);
    }
  }

  for(int s=0;s<Ls;s++){
    ExtractSlice(prop[s], prop5, s , 0);
    ExtractSlice(proprev[s], prop5rev, s , 0);
  }

  Gamma g5(Gamma::Algebra::Gamma5);
  LatticeComplex    C(UGrid); 
  std::vector<LatticeComplex> PAmu(Nd,UGrid); 
  LatticePropagator p5d(UGrid); 
  LatticePropagator us_p5d(UGrid); 
  LatticePropagator gp5d(UGrid); 
  LatticePropagator gus_p5d(UGrid); 

#define Pp(Q) (0.5*(Q+g5*Q))
#define Pm(Q) (0.5*(Q-g5*Q))
#define Q_4d(Q) (Pm((Q)[0]) + Pp((Q)[Ls-1]))
#define TopRowWithSource(Q) (phys_src + (1.0-mass)*Q_4d(Q))

  std::vector<LatticePropagator> L_Q(Ls,UGrid); L_Q=proprev; // shorthand name
  std::vector<LatticePropagator> R_Q(Ls,UGrid); R_Q=prop; // shorthand name

  LatticePropagator L_TmLsGq0(UGrid); 
  LatticePropagator L_TmLsTmp(UGrid);
  LatticePropagator R_TmLsGq0(UGrid); 
  LatticePropagator R_TmLsTmp(UGrid);
  {
    LatticePropagator TermA(UGrid);
    LatticePropagator TermB(UGrid);
    LatticePropagator TermC(UGrid);
    LatticePropagator TermD(UGrid);
    TermA = (Pp(Q_4d(L_Q)));
    TermB = (Pm(Q_4d(L_Q)));
    TermC = (Pm(TopRowWithSource(L_Q)));
    TermD = (Pp(TopRowWithSource(L_Q)));

    L_TmLsGq0 = (TermD - TermA + TermB);
    L_TmLsTmp = (TermC - TermB + TermA);

    TermA = (Pp(Q_4d(R_Q)));
    TermB = (Pm(Q_4d(R_Q)));
    TermC = (Pm(TopRowWithSource(R_Q)));
    TermD = (Pp(TopRowWithSource(R_Q)));

    R_TmLsGq0 = (TermD - TermA + TermB);
    R_TmLsTmp = (TermC - TermB + TermA);
  }

  std::vector<LatticePropagator> R_TmLsGq(Ls,UGrid);
  std::vector<LatticePropagator> L_TmLsGq(Ls,UGrid);
  for(int s=0;s<Ls;s++){
    R_TmLsGq[s] = (Pm((R_Q)[(s)]) + Pp((R_Q)[((s)-1+Ls)%Ls]));
    L_TmLsGq[s] = (Pm((L_Q)[(s)]) + Pp((L_Q)[((s)-1+Ls)%Ls]));
  }
  
  for(int mu=0;mu<Nd;mu++){
    PAmu[mu]=Zero();
    Gamma gmu=Gamma(Gmu[mu]);

    for(int s=0;s<Ls;s++){

      int sp = (s+1)%Ls;
      int sr = Ls-1-s;
      int srp= (sr+1)%Ls;

      // Mobius parameters
      auto b=Ddwf.bs[s];
      auto c=Ddwf.cs[s];
      assert(Ddwfrev.bs[sr]==Ddwf.bs[s]);
      assert(Ddwfrev.cs[sr]==Ddwf.cs[s]);

      LatticePropagator tmp(UGrid); 

      if (s == 0) {
	p5d    =(b*Pm(L_TmLsGq[Ls-1])+ c*Pp(L_TmLsGq[Ls-1]) + b*Pp(L_TmLsTmp)   + c*Pm(L_TmLsTmp     ));
	tmp    =(b*Pm(R_TmLsGq0)     + c*Pp(R_TmLsGq0     ) + b*Pp(R_TmLsGq[1]) + c*Pm(R_TmLsGq[1]));
	us_p5d = peekLorentz(Umu,mu)*Cshift(tmp,mu,1);
      } else if (s == Ls-1) {
	p5d    =(b*Pm(L_TmLsGq0)     + c*Pp(L_TmLsGq0     ) + b*Pp(L_TmLsGq[1]) + c*Pm(L_TmLsGq[1]));
	tmp    =(b*Pm(R_TmLsGq[Ls-1])+ c*Pp(R_TmLsGq[Ls-1]) + b*Pp(R_TmLsTmp)   + c*Pm(R_TmLsTmp   ));
	us_p5d = peekLorentz(Umu,mu)*Cshift(tmp,mu,1);
      } else {
	p5d    =(b*Pm(L_TmLsGq[sr]) + c*Pp(L_TmLsGq[sr])+ b*Pp(L_TmLsGq[srp])+ c*Pm(L_TmLsGq[srp]));
	tmp    =(b*Pm(R_TmLsGq[s])  + c*Pp(R_TmLsGq[s]) + b*Pp(R_TmLsGq[sp ])+ c*Pm(R_TmLsGq[sp]));
	us_p5d = peekLorentz(Umu,mu)*Cshift(tmp,mu,1);
      }
              
      gp5d=g5*p5d;
      gus_p5d=gmu*us_p5d;
      auto bpc = 0.5/(b+c);
      C = bpc*localInnerProduct(gp5d,gus_p5d);
      C-= bpc*localInnerProduct(gp5d,us_p5d);

      if (s == 0) {
	p5d    =(b*Pm(R_TmLsGq0)     + c*Pp(R_TmLsGq0  )    + b*Pp(R_TmLsGq[1]) + c*Pm(R_TmLsGq[1]));
	tmp    =(b*Pm(L_TmLsGq[Ls-1])+ c*Pp(L_TmLsGq[Ls-1]) + b*Pp(L_TmLsTmp)   + c*Pm(L_TmLsTmp  ));
	us_p5d = peekLorentz(Umu,mu)*Cshift(tmp,mu,1);
      } else if (s == Ls-1) {
	p5d    =(b*Pm(R_TmLsGq[Ls-1])+ c*Pp(R_TmLsGq[Ls-1]) + b*Pp(R_TmLsTmp)   + c*Pm(R_TmLsTmp  ));
	tmp    =(b*Pm(L_TmLsGq0)     + c*Pp(L_TmLsGq0  )    + b*Pp(L_TmLsGq[1]) + c*Pm(L_TmLsGq[1]));
	us_p5d = peekLorentz(Umu,mu)*Cshift(tmp,mu,1);
      } else {
	p5d    =(b*Pm(R_TmLsGq[s])  + c*Pp(R_TmLsGq[s])  + b*Pp(R_TmLsGq[sp ])+ c*Pm(R_TmLsGq[sp]));
	tmp    =(b*Pm(L_TmLsGq[sr]) + c*Pp(L_TmLsGq[sr]) + b*Pp(L_TmLsGq[srp])+ c*Pm(L_TmLsGq[srp]));
	us_p5d = peekLorentz(Umu,mu)*Cshift(tmp,mu,1);
      }

      gp5d=gmu*p5d;
      gus_p5d=g5*us_p5d;

      bpc = 0.5/(b+c);
      C+= bpc*localInnerProduct(gus_p5d,gp5d);
      C+= bpc*localInnerProduct(gus_p5d,p5d);

      if (s < Ls/2) PAmu[mu] -= C;
      else          PAmu[mu] += C;
              
    }
  }

  std::cout << "done "<<std::endl;
  LatticePropagator psi(UGrid);
  psi = (prop[Ls/2-1]+g5*prop[Ls/2-1] +prop[Ls/2]  -g5*prop[Ls/2]   )*0.5;
  PJ5q=localInnerProduct(psi,psi);
  std::cout << " J5qref "<<norm2(PJ5q)<<std::endl;

  std::cout << " DmuAmu "<<std::endl;
  LatticeComplex Defect(UGrid);
  Defect = Zero();
  for(int mu=0;mu<Nd;mu++) {
    Defect = Defect + PAmu[mu]-Cshift(PAmu[mu],mu,-1);
  }
  Ddwf.ContractConservedCurrent(prop5rev,prop5,Axial_mu,phys_src,Current::Axial,Tdir);
  PA       = trace(g5*Axial_mu);
  PP       = trace(adj(prop4)*prop4);

  Defect = Defect - 2.0*Ddwf.mass* PP;
  Defect = Defect - 2.0*PJ5q;

  std::vector<TComplex> sumPAref;
  std::vector<TComplex> sumPA;
  std::vector<TComplex> sumPP;
  std::vector<TComplex> sumPJ5qref;
  std::vector<TComplex> sumPJ5q;
  std::vector<TComplex> sumDefect;

  // Spatial sum
  sliceSum(PAmu[Tdir],sumPAref,Tdir);
  sliceSum(PA,sumPA,Tdir);
  sliceSum(PJ5q,sumPJ5q,Tdir);
  sliceSum(PP,sumPP,Tdir);
  sliceSum(Defect,sumDefect,Tdir);
  
  Ddwf.ContractJ5q(prop5,PJ5q);
  sliceSum(PJ5q,sumPJ5qref,Tdir);

  int Nt=sumPA.size();
  for(int t=0;t<Nt;t++){
    std::cout <<t<<" PAc reference "<<real(TensorRemove(sumPAref[t]));
    std::cout    <<" PAc action    "<<real(TensorRemove(sumPA[t]));
    std::cout    <<" PJ5q ref      "<<real(TensorRemove(sumPJ5qref[t]));
    std::cout    <<" PJ5q action   "<<real(TensorRemove(sumPJ5q[t]));
    std::cout <<"WTI defects "<<real(TensorRemove(sumPAref[t]-sumPAref[(t-1+Nt)%Nt] - 2.0*(Ddwf.mass*sumPP[t] + sumPJ5q[t]) ))<<",";
    std::cout <<real(TensorRemove(sumPA[t]-sumPA[(t-1+Nt)%Nt] - 2.0*(Ddwf.mass*sumPP[t] + sumPJ5q[t]) ))<<"\n";
  }
}
#endif
      // Verify solution with independent true residual
      LatticeGaugeField Umu5d(FGrid); 
      std::vector<LatticeColourMatrix> U(4,FGrid);
      {
	auto Umu5d_v = Umu5d.View();
	auto Umu_v = Umu.View();
	for(int ss=0;ss<Umu.Grid()->oSites();ss++){
	  for(int s=0;s<Ls;s++){
	    Umu5d_v[Ls*ss+s] = Umu_v[ss];
	  }
	}
      }
      for(int mu=0;mu<Nd;mu++){
	U[mu] = PeekIndex<LorentzIndex>(Umu5d,mu);
      }
      LatticeFermion ref(FGrid);
      LatticeFermion tmp(FGrid);
      ref = Zero();
      for(int mu=0;mu<Nd;mu++){
        tmp = U[mu]*Cshift(result5,mu+1,1);
	ref=ref + tmp - Gamma(Gmu[mu])*tmp;

	tmp =adj(U[mu])*result5;
	tmp =Cshift(tmp,mu+1,-1);
	ref=ref + tmp + Gamma(Gmu[mu])*tmp;
      }
      ref = -0.5*ref;
      // Dperp
      {
	RealD diag = 5.0 - Ddwf.M5;
	mass = Ddwf.mass;
	auto psi=result5.View();
	auto chi=tmp.View();
	thread_for(sss,UGrid->oSites(),{
	  uint64_t ss= sss*Ls;
	  typedef vSpinColourVector spinor;
	  spinor tmp1, tmp2;
	  for(int s=0;s<Ls;s++){
	    uint64_t idx_u = ss+((s+1)%Ls);
	    uint64_t idx_l = ss+((s+Ls-1)%Ls);
	    spProj5m(tmp1,psi(idx_u));
	    spProj5p(tmp2,psi(idx_l));
	    double pu = (s==(Ls-1)) ? mass: -1.0;
	    double pl = (s==0) ? mass: -1.0;
	    chi[ss+s]=diag*psi(ss+s)+pu*tmp1+pl*tmp2;
	  }
 	});
      }
      ref = ref + tmp;
      ref = ref - src5;
      std::cout << "residual "<< norm2(ref)<< std::endl;
      std::cout << "src      "<< norm2(src5)<< std::endl;
      std::cout << "result   "<< norm2(result5)<< std::endl;
      */
