/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_general_coarse_hdcg.cc

    Copyright (C) 2023

Author: Peter Boyle <pboyle@bnl.gov>

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
#include <Grid/algorithms/iterative/ImplicitlyRestartedBlockLanczos.h>
#include <Grid/algorithms/iterative/ImplicitlyRestartedBlockLanczosCoarse.h>
#include <Grid/algorithms/iterative/AdefMrhs.h>
#include <Grid/algorithms/iterative/PowerSpectrum.h>
#include <Grid/algorithms/iterative/BlockConjugateGradient.h>

using namespace std;
using namespace Grid;

template<class aggregation>
void SaveFineEvecs(aggregation &Agg,std::string file)
{
#ifdef HAVE_LIME
  emptyUserRecord record;
  ScidacWriter WR(Agg[0].Grid()->IsBoss());
  WR.open(file);
  for(int b=0;b<Agg.size();b++){
    WR.writeScidacFieldRecord(Agg[b],record,0,Grid::BinaryIO::BINARYIO_LEXICOGRAPHIC);
  }
  WR.close();
#endif
}
template<class aggregation>
void SaveBasis(aggregation &Agg,std::string file)
{
#ifdef HAVE_LIME
  emptyUserRecord record;
  ScidacWriter WR(Agg.FineGrid->IsBoss());
  WR.open(file);
  for(int b=0;b<Agg.subspace.size();b++){
    WR.writeScidacFieldRecord(Agg.subspace[b],record,0,Grid::BinaryIO::BINARYIO_LEXICOGRAPHIC);
    //    WR.writeScidacFieldRecord(Agg.subspace[b],record);
  }
  WR.close();
#endif
}
template<class aggregation>
void LoadBasis(aggregation &Agg, std::string file)
{
#ifdef HAVE_LIME
  emptyUserRecord record;
  ScidacReader RD ;
  RD.open(file);
  for(int b=0;b<Agg.subspace.size();b++){
    RD.readScidacFieldRecord(Agg.subspace[b],record,Grid::BinaryIO::BINARYIO_LEXICOGRAPHIC);
    //    RD.readScidacFieldRecord(Agg.subspace[b],record,0);
  }    
  RD.close();
#endif
}

template<class aggregation>
void LoadBasisSkip(aggregation &Agg, std::string file,int N,LatticeFermionF & tmp)
{
#ifdef HAVE_LIME
  emptyUserRecord record;
  ScidacReader RD ;
  
  RD.open(file);
  for(int b=0;b<Agg.subspace.size();b++){
    for(int n=0;n<N;n++){
      RD.readScidacFieldRecord(tmp,record,Grid::BinaryIO::BINARYIO_LEXICOGRAPHIC);
      if(n==0) precisionChange(Agg.subspace[b],tmp);
    }
    //    RD.readScidacFieldRecord(Agg.subspace[b],record,0);
  }    
  RD.close();
#endif
}
template<class aggregation>
void LoadBasisSum(aggregation &Agg, std::string file,int N,LatticeFermionF & tmp)
{
#ifdef HAVE_LIME
  emptyUserRecord record;
  ScidacReader RD ;
  
  LatticeFermionF sum(tmp.Grid());
  RD.open(file);
  for(int b=0;b<Agg.subspace.size();b++){
    sum=Zero();
    for(int n=0;n<N;n++){
      RD.readScidacFieldRecord(tmp,record,Grid::BinaryIO::BINARYIO_LEXICOGRAPHIC);
      sum=sum+tmp;
    }
    precisionChange(Agg.subspace[b],sum);
    //    RD.readScidacFieldRecord(Agg.subspace[b],record,0);
  }    
  RD.close();
#endif
}

template<class CoarseVector>
void SaveEigenvectors(std::vector<RealD>            &eval,
		      std::vector<CoarseVector>     &evec,
		      std::string evec_file,
		      std::string eval_file)
{
#ifdef HAVE_LIME
  emptyUserRecord record;
  ScidacWriter WR(evec[0].Grid()->IsBoss());
  WR.open(evec_file);
  for(int b=0;b<evec.size();b++){
    WR.writeScidacFieldRecord(evec[b],record,0,0);
  }
  WR.close();
  XmlWriter WRx(eval_file);
  write(WRx,"evals",eval);
#endif
}
template<class CoarseVector>
void LoadEigenvectors(std::vector<RealD>            &eval,
		      std::vector<CoarseVector>     &evec,
		      std::string evec_file,
		      std::string eval_file)
{
#ifdef HAVE_LIME
    XmlReader RDx(eval_file);
    read(RDx,"evals",eval);
    emptyUserRecord record;

    Grid::ScidacReader RD ;
    RD.open(evec_file);
    assert(evec.size()==eval.size());
    for(int k=0;k<eval.size();k++) {
      RD.readScidacFieldRecord(evec[k],record);
    }
    RD.close();
#endif
}

// Want Op in CoarsenOp to call MatPcDagMatPc
template<class Field>
class HermOpAdaptor : public LinearOperatorBase<Field>
{
  LinearOperatorBase<Field> & wrapped;
public:
  HermOpAdaptor(LinearOperatorBase<Field> &wrapme) : wrapped(wrapme)  {};
  void Op     (const Field &in, Field &out)   { wrapped.HermOp(in,out);  }
  void HermOp(const Field &in, Field &out)    { wrapped.HermOp(in,out); }
  void AdjOp     (const Field &in, Field &out){ wrapped.HermOp(in,out);  }
  void OpDiag (const Field &in, Field &out)                  {    assert(0);  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {    assert(0);  }
  void OpDirAll  (const Field &in, std::vector<Field> &out)  {    assert(0);  };
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){    assert(0);  }
};

template<class Field> class FixedCGPolynomial : public LinearFunction<Field>
{
public:
  using LinearFunction<Field>::operator();
  typedef LinearOperatorBase<Field> FineOperator;
  FineOperator   & _SmootherOperator;
  ConjugateGradientPolynomial<Field>  CG;
  int iters;
  bool record;
  int replay_count;
  FixedCGPolynomial(int _iters, FineOperator &SmootherOperator) :
    _SmootherOperator(SmootherOperator),
    iters(_iters),
    record(true),
    CG(0.0,_iters,false)
  {
    std::cout << GridLogMessage<<" FixedCGPolynomial order "<<iters<<std::endl;
    replay_count = 0;
  };
  void operator() (const Field &in, Field &out) 
  {
#if 1
    GridBase *grid = in.Grid();
    Field Mx0(grid);
    Field r0(grid);
    Field Minvr0(grid);

    _SmootherOperator.HermOp(out,Mx0);

    r0 = in - Mx0;

    Minvr0 = Zero();
    Minvr0.Checkerboard()=in.Checkerboard();
    
    if ( record ) {
      std::cout << " FixedCGPolynomial recording polynomial "<<std::endl;
      CG.Solve(_SmootherOperator,r0,Minvr0);
      record = false;
      /*
      std::cout << "P(x) = 0 "<<std::endl;
      for(int i=0;i<CG.polynomial.size();i++){
	std::cout<<" + "<< CG.polynomial[i]<<" * (x**"<<i<<")"<<std::endl;
	}
      */
      Field tmp(Minvr0.Grid());
      CG.CGsequenceHermOp(_SmootherOperator,r0,tmp);
      tmp = tmp - Minvr0;
      std::cout << " CGsequence error "<<norm2(tmp)<<" / "<<norm2(out)<<std::endl;
    } else {
      std::cout << " FixedCGPolynomial replaying polynomial "<<std::endl;
      CG.CGsequenceHermOp(_SmootherOperator,r0,Minvr0);
      if ( replay_count %5== 0 ) record=true;
      replay_count++;
    }
    out = out + Minvr0;
    _SmootherOperator.HermOp(out,r0);
    r0 = r0 - in;
    RealD rr=norm2(r0);
    RealD ss=norm2(in);
    std::cout << " FixedCGPolynomial replayed polynomial resid "<<::sqrt(rr/ss)<<std::endl;
#else
    out = Zero();
    out.Checkerboard()=in.Checkerboard();
    if ( record ) {
      std::cout << " FixedCGPolynomial recording polynomial "<<std::endl;
      CG.Solve(_SmootherOperator,in,out);
      record = false;
      std::cout << "P(x) = 0 "<<std::endl;
      for(int i=0;i<CG.polynomial.size();i++){
	std::cout<<" + "<< CG.polynomial[i]<<" * (x**"<<i<<")"<<std::endl;
      }
      Field tmp(in.Grid());
      CG.CGsequenceHermOp(_SmootherOperator,in,tmp);
      tmp = tmp - out;
      std::cout << " CGsequence error "<<norm2(tmp)<<" / "<<norm2(out)<<std::endl;
    } else {
      std::cout << " FixedCGPolynomial replaying polynomial "<<std::endl;
      CG.CGsequenceHermOp(_SmootherOperator,in,out);
      if ( replay_count %5== 5 ) record=true;
      replay_count++;
    }
#endif
    
  }
  void operator() (const std::vector<Field> &in, std::vector<Field> &out)
  {
    for(int i=0;i<out.size();i++){
      out[i]=Zero();
    }
    int blockDim = 0;//not used for BlockCGVec
    BlockConjugateGradient<Field>    BCGV  (BlockCGrQVec,blockDim,0.0,iters,false);
    BCGV(_SmootherOperator,in,out);
  }
  
};
template<class Field> class CGSmoother : public LinearFunction<Field>
{
public:
  using LinearFunction<Field>::operator();
  typedef LinearOperatorBase<Field> FineOperator;
  FineOperator   & _SmootherOperator;
  int iters;
  CGSmoother(int _iters, FineOperator &SmootherOperator) :
    _SmootherOperator(SmootherOperator),
    iters(_iters)
  {
    std::cout << GridLogMessage<<" Mirs smoother order "<<iters<<std::endl;
  };
  void operator() (const Field &in, Field &out) 
  {
    ConjugateGradient<Field>  CG(0.0,iters,false); // non-converge is just fine in a smoother

    out=Zero();

    CG(_SmootherOperator,in,out);
  }
};


RealD InverseApproximation(RealD x){
  return 1.0/x;
}
template<class Field> class ChebyshevSmoother : public LinearFunction<Field>
{
public:
  using LinearFunction<Field>::operator();
  typedef LinearOperatorBase<Field> FineOperator;
  FineOperator   & _SmootherOperator;
  Chebyshev<Field> Cheby;
  ChebyshevSmoother(RealD _lo,RealD _hi,int _ord, FineOperator &SmootherOperator) :
    _SmootherOperator(SmootherOperator),
    Cheby(_lo,_hi,_ord,InverseApproximation)
  {
    std::cout << GridLogMessage<<" Chebyshev smoother order "<<_ord<<" ["<<_lo<<","<<_hi<<"]"<<std::endl;
  };
  void operator() (const Field &in, Field &out) 
  {
    //    Field r(out.Grid());
    Cheby(_SmootherOperator,in,out);
    //    _SmootherOperator.HermOp(out,r);
    //    r=r-in;
    //    RealD rr=norm2(r);
    //    RealD ss=norm2(in);
    //    std::cout << GridLogMessage<<" Chebyshev smoother resid "<<::sqrt(rr/ss)<<std::endl;
  }
};

template<class Field> class ChebyshevInverter : public LinearFunction<Field>
{
public:
  using LinearFunction<Field>::operator();
  typedef LinearOperatorBase<Field> FineOperator;
  FineOperator   & _Operator;
  Chebyshev<Field> Cheby;
  ChebyshevInverter(RealD _lo,RealD _hi,int _ord, FineOperator &Operator) :
    _Operator(Operator),
    Cheby(_lo,_hi,_ord,InverseApproximation)
  {
    std::cout << GridLogMessage<<" Chebyshev Inverter order "<<_ord<<" ["<<_lo<<","<<_hi<<"]"<<std::endl;
  };
  void operator() (const Field &in, Field &out) 
  {
    Field r(in.Grid());
    Field AinvR(in.Grid());
    _Operator.HermOp(out,r);
    r = in - r; // b - A x
    Cheby(_Operator,r,AinvR); // A^{-1} ( b - A x ) ~ A^{-1} b - x
    out = out + AinvR;
    _Operator.HermOp(out,r);
    r = in - r; // b - A x
    RealD rr = norm2(r);
    RealD ss = norm2(in);
    std::cout << "ChebshevInverse resid " <<::sqrt(rr/ss)<<std::endl;
  }
};



int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int sample=1;
  if( GridCmdOptionExists(argv,argv+argc,"--sample") ){
    std::string arg;
    arg = GridCmdOptionPayload(argv,argv+argc,"--sample");
    GridCmdOptionInt(arg,sample);
  }
  
  const int Ls=24;
  const int nbasis = 64;
  const int cb = 0 ;
  RealD mass=0.00078;

  if( GridCmdOptionExists(argv,argv+argc,"--mass") ){
    std::string arg;
    arg = GridCmdOptionPayload(argv,argv+argc,"--mass");
    GridCmdOptionFloat(arg,mass);
  }

  RealD M5=1.8;
  RealD b=1.5;
  RealD c=0.5;

  std::cout << GridLogMessage << " *************************** " <<std::endl;
  std::cout << GridLogMessage << " Mass " <<mass<<std::endl;
  std::cout << GridLogMessage << " M5   " <<M5<<std::endl;
  std::cout << GridLogMessage << " Ls   " <<Ls<<std::endl;
  std::cout << GridLogMessage << " b    " <<b<<std::endl;
  std::cout << GridLogMessage << " c    " <<c<<std::endl;
  std::cout << GridLogMessage << " *************************** " <<std::endl;
  
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
								   GridDefaultSimd(Nd,vComplex::Nsimd()),
								   GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  //////////////////////////////////////////
  // Single precision grids -- lanczos + smoother
  //////////////////////////////////////////
  GridCartesian         * UGridF   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
								   GridDefaultSimd(Nd,vComplexF::Nsimd()),
								   GridDefaultMpi());
  GridRedBlackCartesian * UrbGridF = SpaceTimeGrid::makeFourDimRedBlackGrid(UGridF);
  GridCartesian         * FGridF   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGridF);
  GridRedBlackCartesian * FrbGridF = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGridF);

  ///////////////////////// Configuration /////////////////////////////////
  LatticeGaugeField Umu(UGrid);

  FieldMetaData header;
  std::string file("ckpoint_lat.1000");
  NerscIO::readConfiguration(Umu,header,file);

  //////////////////////// Fermion action //////////////////////////////////
  MobiusFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c);
  SchurDiagMooeeOperator<MobiusFermionD, LatticeFermion> HermOpEO(Ddwf);

  std::cout << "**************************************"<<std::endl;
  std::cout << "         Fine Power method            "<<std::endl;
  std::cout << "**************************************"<<std::endl;

  {
    LatticeFermionD pm_src(FrbGrid);
    pm_src = ComplexD(1.0);
    PowerMethod<LatticeFermionD>       fPM;
    fPM(HermOpEO,pm_src);
  }

  if(0)
  {

    std::cout << "**************************************"<<std::endl;
    std::cout << "         Fine Lanczos           "<<std::endl;
    std::cout << "**************************************"<<std::endl;

    typedef LatticeFermionF FermionField;
    LatticeGaugeFieldF UmuF(UGridF);
    precisionChange(UmuF,Umu);
    MobiusFermionF DdwfF(UmuF,*FGridF,*FrbGridF,*UGridF,*UrbGridF,mass,M5,b,c);
    SchurDiagMooeeOperator<MobiusFermionF, LatticeFermionF> HermOpEOF(DdwfF);

    const int Fine_Nstop = 200;
    const int Fine_Nk = 200;
    const int Fine_Np = 200;
    const int Fine_Nm = Fine_Nk+Fine_Np;
    const int Fine_MaxIt= 10;

    RealD Fine_resid = 1.0e-4;
    std::cout << GridLogMessage << "Fine Lanczos "<<std::endl;
    std::cout << GridLogMessage << "Nstop "<<Fine_Nstop<<std::endl;
    std::cout << GridLogMessage << "Nk "<<Fine_Nk<<std::endl;
    std::cout << GridLogMessage << "Np "<<Fine_Np<<std::endl;
    std::cout << GridLogMessage << "resid "<<Fine_resid<<std::endl;

    Chebyshev<FermionField> Cheby(0.002,92.0,401);
    //    Chebyshev<FermionField> Cheby(0.1,92.0,401);
    FunctionHermOp<FermionField> OpCheby(Cheby,HermOpEOF);
    PlainHermOp<FermionField> Op     (HermOpEOF);
    ImplicitlyRestartedLanczos<FermionField> IRL(OpCheby,Op,Fine_Nstop,Fine_Nk,Fine_Nm,Fine_resid,Fine_MaxIt);
    std::vector<RealD>          Fine_eval(Fine_Nm);
    FermionField                Fine_src(FrbGridF); 
    Fine_src = ComplexF(1.0);
    std::vector<FermionField> Fine_evec(Fine_Nm,FrbGridF);

    int Fine_Nconv;
    std::cout << GridLogMessage <<" Calling IRL.calc single prec"<<std::endl;
    IRL.calc(Fine_eval,Fine_evec,Fine_src,Fine_Nconv);

    std::string evec_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/Subspace.phys48.evecF");
    SaveFineEvecs(Fine_evec,evec_file);
  }


  //////////////////////////////////////////
  // Construct a coarsened grid with 4^4 cell
  //////////////////////////////////////////
  Coordinate Block({4,4,6,4});
  Coordinate clatt = GridDefaultLatt();
  for(int d=0;d<clatt.size();d++){
    clatt[d] = clatt[d]/Block[d];
  }

  GridCartesian *Coarse4d =  SpaceTimeGrid::makeFourDimGrid(clatt,
							    GridDefaultSimd(Nd,vComplex::Nsimd()),
							    GridDefaultMpi());;
  GridCartesian *Coarse5d =  SpaceTimeGrid::makeFiveDimGrid(1,Coarse4d);

  ///////////////////////// RNGs /////////////////////////////////
  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  std::vector<int> cseeds({5,6,7,8});

  GridParallelRNG          RNG5(FGrid);   RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);   RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          CRNG(Coarse5d);CRNG.SeedFixedIntegers(cseeds);

  
  typedef HermOpAdaptor<LatticeFermionD> HermFineMatrix;
  HermFineMatrix FineHermOp(HermOpEO);

  ////////////////////////////////////////////////////////////
  ///////////// Coarse basis and Little Dirac Operator ///////
  ////////////////////////////////////////////////////////////
  typedef GeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> LittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector CoarseVector;

  NextToNextToNextToNearestStencilGeometry5D geom(Coarse5d);

  typedef Aggregation<vSpinColourVector,vTComplex,nbasis> Subspace;
  Subspace Aggregates(Coarse5d,FrbGrid,cb);

  ////////////////////////////////////////////////////////////
  // Need to check about red-black grid coarsening
  ////////////////////////////////////////////////////////////
  std::string subspace_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/Subspace.phys48.mixed.2500.60");
  //  //  std::string subspace_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/Subspace.phys48.new.62");
  std::string refine_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/Subspace.phys48.evecF");
  //  std::string refine_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/Refine.phys48.mixed.2500.60");
  std::string ldop_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/LittleDiracOp.phys48.mixed.60");
  std::string evec_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/evecs.scidac");
  std::string eval_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/eval.xml");
  bool load_agg=true;
  bool load_refine=true;
  bool load_mat=false;
  bool load_evec=false;

  int refine=1;
  if ( load_agg ) {
    if ( !(refine) || (!load_refine) ) { 
      LoadBasis(Aggregates,subspace_file);
    }
  } else {
    //    Aggregates.CreateSubspaceMultishift(RNG5,HermOpEO,
    //    					0.0003,1.0e-5,2000); // Lo, tol, maxit
    //    Aggregates.CreateSubspaceChebyshev(RNG5,HermOpEO,nbasis,95.,0.01,1500);// <== last run
    Aggregates.CreateSubspaceChebyshevNew(RNG5,HermOpEO,95.); 
    SaveBasis(Aggregates,subspace_file);
  }

  std::cout << "**************************************"<<std::endl;
  std::cout << "Building MultiRHS Coarse operator"<<std::endl;
  std::cout << "**************************************"<<std::endl;
  ConjugateGradient<CoarseVector>  coarseCG(4.0e-2,20000,true);
    
  const int nrhs=12;
    
  Coordinate mpi=GridDefaultMpi();
  Coordinate rhMpi ({1,1,mpi[0],mpi[1],mpi[2],mpi[3]});
  Coordinate rhLatt({nrhs,1,clatt[0],clatt[1],clatt[2],clatt[3]});
  Coordinate rhSimd({vComplex::Nsimd(),1, 1,1,1,1});
    
  GridCartesian *CoarseMrhs = new GridCartesian(rhLatt,rhSimd,rhMpi); 
  typedef MultiGeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> MultiGeneralCoarsenedMatrix_t;
  MultiGeneralCoarsenedMatrix_t mrhs(geom,CoarseMrhs);

  std::cout << "**************************************"<<std::endl;
  std::cout << "         Coarse Lanczos               "<<std::endl;
  std::cout << "**************************************"<<std::endl;

  typedef HermitianLinearOperator<MultiGeneralCoarsenedMatrix_t,CoarseVector> MrhsHermMatrix;
  Chebyshev<CoarseVector>      IRLCheby(0.005,42.0,301);  // 1 iter
  MrhsHermMatrix MrhsCoarseOp     (mrhs);

  //  CoarseVector pm_src(CoarseMrhs);
  //  pm_src = ComplexD(1.0);
  //  PowerMethod<CoarseVector>       cPM;   cPM(MrhsCoarseOp,pm_src);

  int Nk=192;
  int Nm=384;
  int Nstop=Nk;
  int Nconv_test_interval=1;
  
  ImplicitlyRestartedBlockLanczosCoarse<CoarseVector> IRL(MrhsCoarseOp,
							  Coarse5d,
							  CoarseMrhs,
							  nrhs,
							  IRLCheby,
							  Nstop,
							  Nconv_test_interval,
							  nrhs,
							  Nk,
							  Nm,
							  1e-5,10);

  int Nconv;
  std::vector<RealD>            eval(Nm);
  std::vector<CoarseVector>     evec(Nm,Coarse5d);
  std::vector<CoarseVector>     c_src(nrhs,Coarse5d);

  ///////////////////////
  // Deflation guesser object
  ///////////////////////
  MultiRHSDeflation<CoarseVector> MrhsGuesser;

  //////////////////////////////////////////
  // Block projector for coarse/fine
  //////////////////////////////////////////
  MultiRHSBlockProject<LatticeFermionD> MrhsProjector;

  //////////////////////////
  // Extra HDCG parameters
  //////////////////////////
  int maxit=300;
  ConjugateGradient<CoarseVector>  CG(5.0e-2,maxit,false);
  ConjugateGradient<CoarseVector>  CGstart(5.0e-2,maxit,false);
  RealD lo=2.0;
  int ord = 7;
  //  int ord = 11;

  int blockDim = 0;//not used for BlockCG
  BlockConjugateGradient<CoarseVector>    BCG  (BlockCGrQ,blockDim,5.0e-5,maxit,true);

  DoNothingGuesser<CoarseVector> DoNothing;
  //  HPDSolver<CoarseVector> HPDSolveMrhs(MrhsCoarseOp,CG,DoNothing);
  //  HPDSolver<CoarseVector> HPDSolveMrhsStart(MrhsCoarseOp,CGstart,DoNothing);
  //  HPDSolver<CoarseVector> HPDSolveMrhs(MrhsCoarseOp,BCG,DoNothing);
  //  HPDSolver<CoarseVector> HPDSolveMrhsRefine(MrhsCoarseOp,BCG,DoNothing);
  //  FixedCGPolynomial<CoarseVector>  HPDSolveMrhs(maxit,MrhsCoarseOp);

  ChebyshevInverter<CoarseVector> HPDSolveMrhs(1.0e-2,40.0,120,MrhsCoarseOp);  //
  //  ChebyshevInverter<CoarseVector> HPDSolveMrhs(1.0e-2,40.0,110,MrhsCoarseOp);  // 114 iter with Chebysmooth and BlockCG
  //  ChebyshevInverter<CoarseVector> HPDSolveMrhs(1.0e-2,40.0,120,MrhsCoarseOp); // 138 iter with Chebysmooth
  //  ChebyshevInverter<CoarseVector> HPDSolveMrhs(1.0e-2,40.0,200,MrhsCoarseOp); // 139 iter
  //  ChebyshevInverter<CoarseVector> HPDSolveMrhs(3.0e-3,40.0,200,MrhsCoarseOp); // 137 iter, CG smooth, flex
  //  ChebyshevInverter<CoarseVector> HPDSolveMrhs(1.0e-3,40.0,200,MrhsCoarseOp); // 146 iter, CG smooth, flex
  //  ChebyshevInverter<CoarseVector> HPDSolveMrhs(3.0e-4,40.0,200,MrhsCoarseOp); // 156 iter, CG smooth, flex

  /////////////////////////////////////////////////
  // Mirs smoother
  /////////////////////////////////////////////////
  ShiftedHermOpLinearOperator<LatticeFermionD> ShiftedFineHermOp(HermOpEO,lo);
  //  FixedCGPolynomial<LatticeFermionD> CGsmooth(ord,ShiftedFineHermOp) ;
  //  CGSmoother<LatticeFermionD> CGsmooth(ord,ShiftedFineHermOp) ;
  ChebyshevSmoother<LatticeFermionD> CGsmooth(2.0,92.0,8,HermOpEO) ;
  
  if ( load_refine ) {
    //LoadBasis(Aggregates,refine_file);
    LatticeFermionF conv_tmp(FrbGridF);
    LoadBasisSum(Aggregates,refine_file,sample,conv_tmp);
  } else {
    Aggregates.RefineSubspace(HermOpEO,0.001,1.0e-3,3000); // 172 iters
    SaveBasis(Aggregates,refine_file);
  }
  Aggregates.Orthogonalise();

  std::cout << "**************************************"<<std::endl;
  std::cout << "Coarsen after refine"<<std::endl;
  std::cout << "**************************************"<<std::endl;
  mrhs.CoarsenOperator(FineHermOp,Aggregates,Coarse5d);

  std::cout << "**************************************"<<std::endl;
  std::cout << " Recompute coarse evecs  "<<std::endl;
  std::cout << "**************************************"<<std::endl;
  evec.resize(Nm,Coarse5d);
  eval.resize(Nm);
  for(int r=0;r<nrhs;r++){
    random(CRNG,c_src[r]);
  }
 IRL.calc(eval,evec,c_src,Nconv,LanczosType::irbl);

  std::cout << "**************************************"<<std::endl;
  std::cout << " Reimport coarse evecs  "<<std::endl;
  std::cout << "**************************************"<<std::endl;
  MrhsGuesser.ImportEigenBasis(evec,eval);

  std::cout << "**************************************"<<std::endl;
  std::cout << " Setting up mRHS HDCG"<<std::endl;
  std::cout << "**************************************"<<std::endl;
  MrhsProjector.Allocate(nbasis,FrbGrid,Coarse5d);
  MrhsProjector.ImportBasis(Aggregates.subspace);
      
  std::cout << "**************************************"<<std::endl;
  std::cout << "Calling mRHS HDCG"<<std::endl;
  std::cout << "**************************************"<<std::endl;
  TwoLevelADEF2mrhs<LatticeFermion,CoarseVector>
    HDCGmrhs(1.0e-8, 300,
	     FineHermOp,
	     CGsmooth,
	     HPDSolveMrhs,    // Used in M1
	     HPDSolveMrhs,          // Used in Vstart
	     MrhsProjector,
	     MrhsGuesser,
	     CoarseMrhs);
    
  std::vector<LatticeFermionD> src_mrhs(nrhs,FrbGrid);
  std::vector<LatticeFermionD> res_mrhs(nrhs,FrbGrid);
  LatticeFermionD result_accurate(FrbGrid);
  LatticeFermionD result_sloppy(FrbGrid);
  LatticeFermionD error(FrbGrid);
  LatticeFermionD residual(FrbGrid);

  for(int r=0;r<nrhs;r++){
    random(RNG5,src_mrhs[r]);
    res_mrhs[r]=Zero();
  }
  HDCGmrhs(src_mrhs,res_mrhs);
  result_accurate = res_mrhs[0];

#if 0
  std::vector<RealD> tols({1.0e-3,1.0e-4,1.0e-5});

  std::vector<RealD>   bins({1.0e-3,1.0e-2,1.0e-1,1.0,10.0,100.0});
  std::vector<int>   orders({6000  ,4000  ,1000  ,500,500 ,500});
  PowerSpectrum GraphicEqualizer(bins,orders);

  for(auto tol : tols) {
    
    TwoLevelADEF2mrhs<LatticeFermion,CoarseVector>
      HDCGmrhsSloppy(tol, 500,
		     FineHermOp,
		     CGsmooth,
		     HPDSolveMrhs,    // Used in M1
		     HPDSolveMrhs,    // Used in Vstart
		     MrhsProjector,
		     MrhsGuesser,
		     CoarseMrhs);
  
    //  Solve again to 10^-5
    for(int r=0;r<nrhs;r++){
      res_mrhs[r]=Zero();
    }
    HDCGmrhsSloppy(src_mrhs,res_mrhs);
    
    result_sloppy = res_mrhs[0];
    error = result_sloppy - result_accurate;
    FineHermOp.HermOp(result_sloppy,residual);
    residual = residual - src_mrhs[0];
    
    std::cout << "**************************************"<<std::endl;
    std::cout << GridLogMessage << " Converged to tolerance "<< tol<<std::endl;
    std::cout << GridLogMessage << " Absolute error "<<norm2(error)<<std::endl;
    std::cout << GridLogMessage << " Residual       "<<norm2(residual)<<std::endl;
    std::cout << "**************************************"<<std::endl;

    std::cout << "**************************************"<<std::endl;
    std::cout << GridLogMessage << " PowerSpectrum of error   "<<std::endl;
    std::cout << "**************************************"<<std::endl;
    GraphicEqualizer(FineHermOp,error);
    std::cout << "**************************************"<<std::endl;
    std::cout << GridLogMessage << " PowerSpectrum of residual   "<<std::endl;
    std::cout << "**************************************"<<std::endl;
    GraphicEqualizer(FineHermOp,residual);

  };
#endif
  
  // Standard CG
#if 0
  {
  std::cout << "**************************************"<<std::endl;
  std::cout << "Calling red black CG"<<std::endl;
  std::cout << "**************************************"<<std::endl;
      
    LatticeFermion result(FrbGrid); result=Zero();
    LatticeFermion    src(FrbGrid); random(RNG5,src);
    result=Zero();

    ConjugateGradient<LatticeFermionD>  CGfine(1.0e-8,30000,false);
    CGfine(HermOpEO, src, result);
  }
#endif  
  Grid_finalize();
  return 0;
}
