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

using namespace std;
using namespace Grid;

template<class Coarsened>
void SaveOperator(Coarsened &Operator,std::string file)
{
#ifdef HAVE_LIME
  emptyUserRecord record;
  ScidacWriter WR(Operator.Grid()->IsBoss());
  assert(Operator._A.size()==Operator.geom.npoint);
  WR.open(file);
  for(int p=0;p<Operator._A.size();p++){
    auto tmp = Operator.Cell.Extract(Operator._A[p]);
    WR.writeScidacFieldRecord(tmp,record,0,0);
    //    WR.writeScidacFieldRecord(tmp,record,0,BINARYIO_LEXICOGRAPHIC);
  }
  WR.close();
#endif
}
template<class Coarsened>
void LoadOperator(Coarsened &Operator,std::string file)
{
#ifdef HAVE_LIME
  emptyUserRecord record;
  Grid::ScidacReader RD ;
  RD.open(file);
  assert(Operator._A.size()==Operator.geom.npoint);
  for(int p=0;p<Operator.geom.npoint;p++){
    conformable(Operator._A[p].Grid(),Operator.CoarseGrid());
    //    RD.readScidacFieldRecord(Operator._A[p],record,BINARYIO_LEXICOGRAPHIC);
    RD.readScidacFieldRecord(Operator._A[p],record,0);
  }    
  RD.close();
  Operator.ExchangeCoarseLinks();
#endif
}
template<class Coarsened>
void ReLoadOperator(Coarsened &Operator,std::string file)
{
#ifdef HAVE_LIME
  emptyUserRecord record;
  Grid::ScidacReader RD ;
  RD.open(file);
  assert(Operator._A.size()==Operator.geom.npoint);
  for(int p=0;p<Operator.geom.npoint;p++){
    auto tmp=Operator.Cell.Extract(Operator._A[p]);
    RD.readScidacFieldRecord(tmp,record,0);
    Operator._A[p] = Operator.Cell.ExchangePeriodic(tmp);
  }    
  RD.close();
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
    //WR.writeScidacFieldRecord(Agg.subspace[b],record,0,BINARYIO_LEXICOGRAPHIC);
    WR.writeScidacFieldRecord(Agg.subspace[b],record,0,0);
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
    //    RD.readScidacFieldRecord(Agg.subspace[b],record,BINARYIO_LEXICOGRAPHIC);
    RD.readScidacFieldRecord(Agg.subspace[b],record,0);
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

RealD InverseApproximation(RealD x){
  return 1.0/x;
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
    Field tmp(in.Grid());
    tmp = in;
    Cheby(_SmootherOperator,tmp,out);
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


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls=24;
  const int nbasis = 62;
  //  const int nbasis = 56;
  //  const int nbasis = 44;
  //  const int nbasis = 36;
  const int cb = 0 ;
  RealD mass=0.00078;
  RealD M5=1.8;
  RealD b=1.5;
  RealD c=0.5;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
								   GridDefaultSimd(Nd,vComplex::Nsimd()),
								   GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  // Construct a coarsened grid with 4^4 cell
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

  ///////////////////////// Configuration /////////////////////////////////
  LatticeGaugeField Umu(UGrid);
  MemoryManager::Print();

  FieldMetaData header;
  std::string file("ckpoint_lat.1000");
  NerscIO::readConfiguration(Umu,header,file);
  MemoryManager::Print();

  //////////////////////// Fermion action //////////////////////////////////
  MobiusFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c);

  SchurDiagMooeeOperator<MobiusFermionD, LatticeFermion> HermOpEO(Ddwf);

  typedef HermOpAdaptor<LatticeFermionD> HermFineMatrix;
  HermFineMatrix FineHermOp(HermOpEO);

  // Run power method on FineHermOp
  //  PowerMethod<LatticeFermion>       PM;   PM(HermOpEO,src);
 
  ////////////////////////////////////////////////////////////
  ///////////// Coarse basis and Little Dirac Operator ///////
  ////////////////////////////////////////////////////////////
  typedef GeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> LittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector CoarseVector;

  NextToNextToNextToNearestStencilGeometry5D geom(Coarse5d);
  
  // Warning: This routine calls PVdagM.Op, not PVdagM.HermOp
  typedef Aggregation<vSpinColourVector,vTComplex,nbasis> Subspace;
  Subspace Aggregates(Coarse5d,FrbGrid,cb);

  ////////////////////////////////////////////////////////////
  // Need to check about red-black grid coarsening
  ////////////////////////////////////////////////////////////
  LittleDiracOperator LittleDiracOp(geom,FrbGrid,Coarse5d);

  std::string subspace_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/Subspace.phys48.rat.18node.62");
  std::string refine_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/Refine.phys48.rat.18node.62");
  std::string ldop_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/LittleDiracOp.phys48.rat.18node.62");
  std::string evec_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/evecs.scidac");
  std::string eval_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/eval.xml");
  bool load_agg=true;
  bool load_refine=true;
  bool load_mat=true;
  bool load_evec=false;
  MemoryManager::Print();

  int refine=1;
  if ( load_agg ) {
    if ( !(refine) || (!load_refine) ) { 
      LoadBasis(Aggregates,subspace_file);
    }
  } else {
    Aggregates.CreateSubspaceMultishift(RNG5,HermOpEO,
					0.0003,1.0e-5,2000); // Lo, tol, maxit

    //    Aggregates.CreateSubspaceChebyshev(RNG5,HermOpEO,nbasis,95.,0.01,1500); <== last run
    SaveBasis(Aggregates,subspace_file);
  }

  if(refine){
    if ( load_refine ) {
      LoadBasis(Aggregates,refine_file);
    } else {
      // HDCG used Pcg to refine
      Aggregates.RefineSubspace(HermOpEO,0.001,1.0e-3,3000);
      SaveBasis(Aggregates,refine_file);
    }
  }

  Aggregates.Orthogonalise();
  if ( load_mat ) {
    LoadOperator(LittleDiracOp,ldop_file);
  } else {
    LittleDiracOp.CoarsenOperator(FineHermOp,Aggregates);
    //    SaveOperator(LittleDiracOp,ldop_file);
  }
  
  // I/O test:
  CoarseVector c_src(Coarse5d);   random(CRNG,c_src);
  CoarseVector c_res(Coarse5d); 
  CoarseVector c_ref(Coarse5d);

  if (0){
    ///////////////////////////////////////////////////
    // Test the operator
    ///////////////////////////////////////////////////
    CoarseVector c_proj(Coarse5d);
    LatticeFermionD    tmp(FrbGrid);
    LatticeFermionD    prom(FrbGrid);
    
    blockPromote(c_src,prom,Aggregates.subspace);

    FineHermOp.HermOp(prom,tmp);

    std::cout<<GridLogMessage<<" Calling big dirac op "<<norm2(tmp)<<std::endl;
    blockProject(c_proj,tmp,Aggregates.subspace);

    std::cout<<GridLogMessage<<" Calling little Dirac Op "<<std::endl;

    LittleDiracOp.M(c_src,c_res);

    std::cout<<GridLogMessage<<"Little dop : "<<norm2(c_res)<<std::endl;
    std::cout<<GridLogMessage<<"Big dop in subspace : "<<norm2(c_proj)<<std::endl;

    c_proj = c_proj - c_res;
    std::cout<<GridLogMessage<<" ldop error: "<<norm2(c_proj)<<std::endl;
  }

  //////////////////////////////////////
  // mrhs coarse operator
  //  Create a higher dim coarse grid
  //////////////////////////////////////////////////////////////////////////////////////

  std::cout << "**************************************"<<std::endl;
  std::cout << "Building MultiRHS Coarse operator"<<std::endl;
  std::cout << "**************************************"<<std::endl;
  ConjugateGradient<CoarseVector>  coarseCG(4.0e-2,20000,true);
    
  const int nrhs=vComplex::Nsimd()*3;
    
  Coordinate mpi=GridDefaultMpi();
  Coordinate rhMpi ({1,1,mpi[0],mpi[1],mpi[2],mpi[3]});
  Coordinate rhLatt({nrhs,1,clatt[0],clatt[1],clatt[2],clatt[3]});
  Coordinate rhSimd({vComplex::Nsimd(),1, 1,1,1,1});
    
  GridCartesian *CoarseMrhs = new GridCartesian(rhLatt,rhSimd,rhMpi); 
  //  MultiGeneralCoarsenedMatrix mrhs(LittleDiracOp,CoarseMrhs);
  typedef MultiGeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> MultiGeneralCoarsenedMatrix_t;
  MultiGeneralCoarsenedMatrix_t mrhs(geom,CoarseMrhs);
  //  mrhs.CopyMatrix(LittleDiracOp);
  //  mrhs.SetMatrix(LittleDiracOp.);
  mrhs.CoarsenOperator(FineHermOp,Aggregates,Coarse5d);
  //  mrhs.CheckMatrix(LittleDiracOp);
  
  //////////////////////////////////////////
  // Build a coarse lanczos
  //////////////////////////////////////////
  std::cout << "**************************************"<<std::endl;
  std::cout << "Building Coarse Lanczos               "<<std::endl;
  std::cout << "**************************************"<<std::endl;

  typedef HermitianLinearOperator<LittleDiracOperator,CoarseVector> HermMatrix;
  HermMatrix CoarseOp     (LittleDiracOp);

  int Nk=192;
  int Nm=256;
  int Nstop=Nk;
  
  Chebyshev<CoarseVector>      IRLCheby(0.005,40.0,201);  // 1 iter
  FunctionHermOp<CoarseVector> IRLOpCheby(IRLCheby,CoarseOp);
  PlainHermOp<CoarseVector>    IRLOp    (CoarseOp);
  
  ImplicitlyRestartedLanczos<CoarseVector> IRL(IRLOpCheby,IRLOp,Nstop,Nk,Nm,1e-5,10);

  int Nconv;
  std::vector<RealD>            eval(Nm);
  std::vector<CoarseVector>     evec(Nm,Coarse5d);

  PowerMethod<CoarseVector>       cPM;   cPM(CoarseOp,c_src);

  if ( load_evec ) {
    eval.resize(Nstop);
    evec.resize(Nstop,Coarse5d);
    LoadEigenvectors(eval,evec,evec_file,eval_file);
  } else { 
    IRL.calc(eval,evec,c_src,Nconv);
    assert(Nstop==eval.size());
    SaveEigenvectors(eval,evec,evec_file,eval_file);
  }

  DeflatedGuesser<CoarseVector> DeflCoarseGuesser(evec,eval);

  MultiRHSDeflation<CoarseVector> MrhsGuesser;
  
  //////////////////////////////////////////
  // Build a coarse space solver
  //////////////////////////////////////////
  int maxit=30000;
  ConjugateGradient<CoarseVector>  CG(1.0e-10,maxit,false);
  ConjugateGradient<LatticeFermionD>  CGfine(1.0e-8,30000,false);
  ZeroGuesser<CoarseVector> CoarseZeroGuesser;
  
  HPDSolver<CoarseVector> HPDSolve(CoarseOp,CG,DeflCoarseGuesser);
  c_res=Zero();

  /////////// MRHS test .////////////
  typedef HermitianLinearOperator<MultiGeneralCoarsenedMatrix_t,CoarseVector> MrhsHermMatrix;
  MrhsHermMatrix MrhsCoarseOp     (mrhs);

#if 1
  { 
    CoarseVector rh_res(CoarseMrhs);
    CoarseVector rh_guess(CoarseMrhs);
    CoarseVector rh_src(CoarseMrhs);

    rh_res= Zero();
    rh_guess= Zero();

    std::cout << "*************************"<<std::endl;
    std::cout << " MrhsGuesser importing"<<std::endl;
    std::cout << "*************************"<<std::endl;
    MrhsGuesser.ImportEigenBasis(evec,eval);
    std::vector<CoarseVector> BlasGuess(nrhs,Coarse5d);
    std::vector<CoarseVector> BlasSource(nrhs,Coarse5d);
    for(int r=0;r<nrhs;r++){
      random(CRNG,BlasSource[r]);
    }

    MrhsGuesser.DeflateSources(BlasSource,BlasGuess);

    for(int r=0;r<nrhs;r++){
      std::cout << "*************************"<<std::endl;
      std::cout << "**** DeflCoarseGuesser &&&&& "<<std::endl;
      std::cout << "*************************"<<std::endl;
      c_src=BlasSource[r];
      DeflCoarseGuesser(c_src,c_res);
      std::cout << "Deflated guess      "<< norm2(c_res)<<std::endl;
      std::cout << "Blas deflated guess "<< norm2(BlasGuess[r])<<std::endl;
      std::cout << "*************************"<<std::endl;
      BlasGuess[r] = BlasGuess[r] - c_res;
      std::cout << "Diff " <<norm2(BlasGuess[r])<<std::endl;
      std::cout << "*************************"<<std::endl;
      InsertSlice(c_res,rh_res,r,0);
      InsertSlice(c_res,rh_guess,r,0);
      InsertSlice(c_src,rh_src,r,0);
    }

    std::cout << " Calling the multiRHS coarse CG"<<std::endl;
    coarseCG(MrhsCoarseOp,rh_src,rh_res);

    //redo with block CG ?
    for(int r=0;r<nrhs;r++){
      std::cout << " compare to single RHS "<<r<<"/"<<nrhs<<std::endl;
      ExtractSlice(c_src,rh_src,r,0);
      ExtractSlice(c_res,rh_res,r,0);
      ExtractSlice(c_ref,rh_guess,r,0);
      coarseCG(CoarseOp,c_src,c_ref);
      std::cout << " mrhs [" <<r <<"] "<< norm2(c_res)<<std::endl;
      std::cout << " srhs [" <<r <<"] "<< norm2(c_ref)<<std::endl;
      c_ref=c_ref-c_res;
      RealD diff =norm2(c_ref)/norm2(c_src);
      std::cout << r << " diff " << diff<<std::endl;
      assert(diff < 1.0e-1);
    }
  }
#endif

  //////////////////////////////////////
  // fine solve
  //////////////////////////////////////
  
  std::vector<RealD> los({2.0});
  std::vector<int> ords({7}); 

 /*
 Powerlaw setup 62 vecs
slurm-1494943.out:Grid : Message : 4874.186617 s : HDCG: Pcg converged in 171 iterations and 1706.548006 s 1.0 32
slurm-1494943.out:Grid : Message : 6490.121648 s : HDCG: Pcg converged in 194 iterations and 1616.219654 s 1.0 16

 Cheby setup: 56vecs
 -- CG smoother O(16): 487
 
Power law setup, 56 vecs -- lambda^-5
slurm-1494383.out:Grid : Message : 4377.173265 s : HDCG: Pcg converged in 204 iterations and 1153.548935 s 1.0 32

Power law setup, 56 vecs -- lambda^-3

slurm-1494242.out:Grid : Message : 4370.464814 s : HDCG: Pcg converged in 204 iterations and 1143.494776 s  1.0 32
slurm-1494242.out:Grid : Message : 5432.414820 s : HDCG: Pcg converged in 237 iterations and 1061.455882 s  1.0 16
slurm-1494242.out:Grid : Message : 6588.727977 s : HDCG: Pcg converged in 205 iterations and 1156.565210 s  0.5 32

 Power law setup, 56 vecs -- lambda^-4
 -- CG smoother    O(16): 290
 -- Cheby smoother O(16): 218 -- getting close to the deflation level I expect 169 from BFM paper @O(7) smoother and 64 nbasis

Conclusion: higher order smoother is doing better. Much better. Use a Krylov smoother instead Mirs as in BFM version.
 */
				      //
  MemoryManager::Print();
  for(int l=0;l<los.size();l++){

    RealD lo = los[l];

    for(int o=0;o<ords.size();o++){

      ConjugateGradient<CoarseVector>  CGsloppy(4.0e-2,maxit,false);
      HPDSolver<CoarseVector> HPDSolveSloppy(CoarseOp,CGsloppy,DeflCoarseGuesser);
      
      //    ChebyshevSmoother<LatticeFermionD,HermFineMatrix > Smoother(lo,92,10,FineHermOp); // 36 best case
      ChebyshevSmoother<LatticeFermionD > ChebySmooth(lo,95,ords[o],FineHermOp);  // 311

      RealD MirsShift = lo;
      ShiftedHermOpLinearOperator<LatticeFermionD> ShiftedFineHermOp(HermOpEO,MirsShift);
      CGSmoother<LatticeFermionD> CGsmooth(ords[o],ShiftedFineHermOp) ;
  
      //////////////////////////////////////////
      // Build a HDCG solver
      //////////////////////////////////////////
      TwoLevelADEF2<LatticeFermion,CoarseVector,Subspace>
	HDCG(1.0e-8, 700,
	     FineHermOp,
	     CGsmooth,
	     HPDSolveSloppy,
	     HPDSolve,
	     Aggregates);
      //      result=Zero();
      //      std::cout << "Calling HDCG single RHS"<<std::endl;
      //      HDCG(src,result);

      //////////////////////////////////////////
      // Build a HDCG mrhs solver
      //////////////////////////////////////////
#if 1
  MemoryManager::Print();
      DoNothingGuesser<CoarseVector> DoNothing;
      HPDSolver<CoarseVector> HPDSolveMrhs(MrhsCoarseOp,CG,DoNothing);
      HPDSolver<CoarseVector> HPDSolveMrhsSloppy(MrhsCoarseOp,CGsloppy,DoNothing);
      TwoLevelADEF2mrhs<LatticeFermion,CoarseVector,Subspace>
	HDCGmrhs(1.0e-8, 500,
		 FineHermOp,
		 CGsmooth,
		 //		 HPDSolveSloppy, // Never used
		 //		 HPDSolve,       // Used in Vstart
		 HPDSolveMrhsSloppy,    // Used in M1
		 HPDSolveMrhs,          // Used in Vstart
		 DeflCoarseGuesser, // single RHS guess used in M1
		 CoarseMrhs,        // Grid needed to Mrhs grid
		 Aggregates);

      std::cout << "Calling mRHS HDCG"<<std::endl;
      FrbGrid->Barrier();
      
      std::vector<LatticeFermionD> src_mrhs(nrhs,FrbGrid);
      std::cout << " mRHS source"<<std::endl;
      std::vector<LatticeFermionD> res_mrhs(nrhs,FrbGrid);
      std::cout << " mRHS result"<<std::endl;

  random(RNG5,src_mrhs[0]);
  for(int r=0;r<nrhs;r++){
	if(r>0)src_mrhs[r]=src_mrhs[0];
	res_mrhs[r]=Zero();
	std::cout << "Setup mrhs source "<<r<<std::endl;
  }
  std::cout << "Calling the mRHS HDCG"<<std::endl;
  MemoryManager::Print();
  HDCGmrhs(src_mrhs,res_mrhs);
  MemoryManager::Print();
#endif
    }
  }

  // Standard CG
#if 1
  {
    LatticeFermion result(FrbGrid); result=Zero();
    LatticeFermion    src(FrbGrid); random(RNG5,src);
    result=Zero();
    CGfine(HermOpEO, src, result);
  }
#endif  
  Grid_finalize();
  return 0;
}
