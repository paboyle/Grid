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

using namespace std;
using namespace Grid;

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

  FieldMetaData header;
  std::string file("ckpoint_lat.1000");
  NerscIO::readConfiguration(Umu,header,file);

  //////////////////////// Fermion action //////////////////////////////////
  MobiusFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c);

  SchurDiagMooeeOperator<MobiusFermionD, LatticeFermion> HermOpEO(Ddwf);

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
  //  std::string subspace_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/Subspace.phys48.new.62");
  std::string refine_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/Refine.phys48.mixed.2500.60");
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

  CoarseVector pm_src(CoarseMrhs);
  pm_src = ComplexD(1.0);
  PowerMethod<CoarseVector>       cPM;   cPM(MrhsCoarseOp,pm_src);

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
  int maxit=3000;
  ConjugateGradient<CoarseVector>  CG(5.0e-2,maxit,false);
  RealD lo=2.0;
  int ord = 7;

  DoNothingGuesser<CoarseVector> DoNothing;
  HPDSolver<CoarseVector> HPDSolveMrhs(MrhsCoarseOp,CG,DoNothing);
  HPDSolver<CoarseVector> HPDSolveMrhsRefine(MrhsCoarseOp,CG,DoNothing);

  /////////////////////////////////////////////////
  // Mirs smoother
  /////////////////////////////////////////////////
  RealD MirsShift = lo;
  ShiftedHermOpLinearOperator<LatticeFermionD> ShiftedFineHermOp(HermOpEO,MirsShift);
  CGSmoother<LatticeFermionD> CGsmooth(ord,ShiftedFineHermOp) ;

  
  if ( load_refine ) {
    LoadBasis(Aggregates,refine_file);
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
  std::cout << "Calling mRHS HDCG"<<std::endl;
  std::cout << "**************************************"<<std::endl;
  MrhsProjector.Allocate(nbasis,FrbGrid,Coarse5d);
  MrhsProjector.ImportBasis(Aggregates.subspace);
      
  TwoLevelADEF2mrhs<LatticeFermion,CoarseVector>
    HDCGmrhs(1.0e-8, 500,
	     FineHermOp,
	     CGsmooth,
	     HPDSolveMrhs,    // Used in M1
	     HPDSolveMrhs,          // Used in Vstart
	     MrhsProjector,
	     MrhsGuesser,
	     CoarseMrhs);
    
  std::vector<LatticeFermionD> src_mrhs(nrhs,FrbGrid);
  std::vector<LatticeFermionD> res_mrhs(nrhs,FrbGrid);
  
  for(int r=0;r<nrhs;r++){
    random(RNG5,src_mrhs[r]);
    res_mrhs[r]=Zero();
  }
  
  HDCGmrhs(src_mrhs,res_mrhs);

  // Standard CG
#if 1
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
