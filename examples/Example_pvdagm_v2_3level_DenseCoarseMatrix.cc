/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Example_pvdagm_v2_3level_DenseCoarseMatrix.cc

    Copyright (C) 2026

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

//
// PVdagM three level multigrid on the V2 coarse operator.
//
// STAGE ONE: grids, types, subspace, and the L1 coarsening only. The L2
// chain, the dense bottom and the solves are not here yet.
//
// Differences from Example_pvdagm_mrhs_3level_DenseCoarseMatrix.cc:
//
//  * The coarse space is UNVECTORISED (sComplexD). The fine space stays
//    vectorised. MultiRHSBlockProject carries the mixed layout.
//
//  * One operator, not two. V1 needed GeneralCoarsenedMatrix to coarsen and
//    MultiGeneralCoarsenedMatrix to apply, bridged by CopyMatrix. V2 does
//    both, and single versus multiRHS is SetGrid on the same object with the
//    matrix elements built once.
//
//  * Nrhs is unconstrained. V1 required nrhs % vComplex::Nsimd() == 0 because
//    its multiRHS grid carried the SIMD in the rhs direction.
//
//  * CoarsenOperator takes the subspace vectors, not an Aggregation. It block
//    orthonormalises them IN PLACE -- the vectors are far too large to copy
//    defensively -- so rawNull is taken first and the RAW vectors are what
//    define the L2 null space. Do not insert an Orthogonalise() anywhere:
//    projecting a block-orthonormal vector onto its own block-orthonormalised
//    aggregation gives e_k, and the near null content is silently gone. The
//    ||<psi|psi> - I||_F guard below is what catches that.
//
// Env: LATT LS MASS NBASIS(compile time) NRHS BLOCK COARSEN_BATCH
//      HOT_START CONFIG SUBSPACE_FILE V1_CHECK
//

#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidualNonHermitian.h>

#include <memory>

using namespace std;
using namespace Grid;

// Compile time so it can be cut down for laptop runs: -DNBASIS=8
#ifndef NBASIS
#define NBASIS 60
#endif

RealD mass      = 0.00078;
int   Nrhs      = 12;
int   Ls        = 24;
int   CoarsenBatch = 9;
std::vector<int> lat_size({48,48,48,96});

void ParseEnvironment(void)
{
  if(getenv("MASS"))          mass        = atof(getenv("MASS"));
  if(getenv("NRHS"))          Nrhs        = atoi(getenv("NRHS"));
  if(getenv("LS"))            Ls          = atoi(getenv("LS"));
  if(getenv("COARSEN_BATCH")) CoarsenBatch= atoi(getenv("COARSEN_BATCH"));
  if(getenv("LATT")){
    Coordinate l;
    GridCmdOptionIntVector(std::string(getenv("LATT")),l);
    GRID_ASSERT(l.size()==4);
    for(int d=0;d<4;d++) lat_size[d]=l[d];
  }

  std::cout << GridLogMessage << "PARAM: LATT               "
	    << lat_size[0]<<"."<<lat_size[1]<<"."<<lat_size[2]<<"."<<lat_size[3] << std::endl;
  std::cout << GridLogMessage << "PARAM: LS                 " << Ls           << std::endl;
  std::cout << GridLogMessage << "PARAM: MASS               " << mass         << std::endl;
  std::cout << GridLogMessage << "PARAM: NBASIS             " << NBASIS       << std::endl;
  std::cout << GridLogMessage << "PARAM: NRHS               " << Nrhs         << std::endl;
  std::cout << GridLogMessage << "PARAM: COARSEN_BATCH      " << CoarsenBatch << std::endl;
}

template <class Field>
void saveSubspace(std::vector<Field> &subspace, std::string const fname){
#ifdef HAVE_LIME
  Grid::emptyUserRecord record;
  Grid::ScidacWriter SW(subspace[0].Grid()->IsBoss());
  SW.open(fname);
  for (int k = 0; k < (int)subspace.size(); k++) SW.writeScidacFieldRecord(subspace[k], record);
  SW.close();
#endif
}
template <class Field>
void loadSubspace(std::vector<Field> &subspace, std::string const fname){
#ifdef HAVE_LIME
  Grid::emptyUserRecord record;
  Grid::ScidacReader SR;
  SR.open(fname);
  for (int k = 0; k < (int)subspace.size(); k++) SR.readScidacFieldRecord(subspace[k], record);
  SR.close();
#endif
}

//////////////////////////////////////////////////////////////////////
// A = PV^dag M (non-Hermitian)
//////////////////////////////////////////////////////////////////////
template<class Matrix,class Field>
class PVdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat; Matrix &_PV;
public:
  PVdagMLinearOperator(Matrix &Mat,Matrix &PV): _Mat(Mat),_PV(PV) {};
  void OpDiag (const Field &in, Field &out) { assert(0); }
  void OpDir  (const Field &in, Field &out,int dir,int disp) { assert(0); }
  void OpDirAll  (const Field &in, std::vector<Field> &out){ assert(0); };
  void Op     (const Field &in, Field &out){ Field tmp(in.Grid()); _Mat.M(in,tmp); _PV.Mdag(tmp,out); }
  void AdjOp  (const Field &in, Field &out){ Field tmp(in.Grid()); _PV.M(in,tmp); _Mat.Mdag(tmp,out); }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){ HermOp(in,out); ComplexD d=innerProduct(in,out); n1=real(d); n2=norm2(out); }
  void HermOp(const Field &in, Field &out){ Field tmp(in.Grid()); Op(in,tmp); AdjOp(tmp,out); }
};

//////////////////////////////////////////////////////////////////////
// ||<v|v> - I||_F over a set of coarse vectors. ~0.23 means the raw near
// null content survived the projection; ~sqrt(N_sites) means block
// orthonormal vectors leaked in and every image collapsed to e_k.
//////////////////////////////////////////////////////////////////////
template<class CoarseField>
RealD GramDefect(std::vector<CoarseField> &v)
{
  RealD s2=0.0;
  for(int i=0;i<(int)v.size();i++){
    for(int j=0;j<(int)v.size();j++){
      ComplexD sij=TensorRemove(innerProduct(v[i],v[j]));
      ComplexD d=sij-(i==j?ComplexD(1.0):ComplexD(0.0));
      s2+=real(d)*real(d)+imag(d)*imag(d);
    }
  }
  return std::sqrt(s2);
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);
  ParseEnvironment();

  RealD M5=1.8, b=1.5, c=0.5;
  const int nbasis=NBASIS;
  const int nrhs=Nrhs;
  const int batch=CoarsenBatch;

  Coordinate mpi = GridDefaultMpi();
  Coordinate fsimd= GridDefaultSimd(Nd,vComplex::Nsimd());

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(lat_size,fsimd,mpi);
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  // Level 1 blocking (default 2^4)
  Coordinate clatt = lat_size;
  Coordinate Block({2,2,2,2});
  if ( getenv("BLOCK") ){ GridCmdOptionIntVector(std::string(getenv("BLOCK")),Block); GRID_ASSERT(Block.size()==4); }
  for(int d=0;d<4;d++){ GRID_ASSERT(lat_size[d]%Block[d]==0); clatt[d]=lat_size[d]/Block[d]; }
  std::cout << GridLogMessage << "Block  " << Block  << "  coarse lattice " << clatt << std::endl;

  //////////////////////////////////////////////////////////////////////
  // The coarse space is unvectorised. The 5D coarse grid is built here
  // rather than through SpaceTimeGrid so the SIMD layout is ours.
  //////////////////////////////////////////////////////////////////////
  Coordinate c5latt({1,clatt[0],clatt[1],clatt[2],clatt[3]});
  Coordinate c5simd({1,1,1,1,1});
  Coordinate c5mpi ({1,mpi[0],mpi[1],mpi[2],mpi[3]});
  GridCartesian *Coarse5d = new GridCartesian(c5latt,c5simd,c5mpi);

  // 6D coarse multiRHS grid: rhs is dim 0, undistributed and unvectorised.
  // No divisibility constraint on nrhs, unlike V1.
  Coordinate cmlatt({nrhs,1,clatt[0],clatt[1],clatt[2],clatt[3]});
  Coordinate cmsimd({1,1,1,1,1,1});
  Coordinate cmmpi ({1,1,mpi[0],mpi[1],mpi[2],mpi[3]});
  GridCartesian *CoarseMrhs = new GridCartesian(cmlatt,cmsimd,cmmpi);

  // 6D coarse grid at the coarsening batch, used only while CoarsenOperator
  // runs. The matrix elements survive the change back to nrhs.
  Coordinate cblatt({batch,1,clatt[0],clatt[1],clatt[2],clatt[3]});
  GridCartesian *CoarseBatch = new GridCartesian(cblatt,cmsimd,cmmpi);

  // 6D fine grid carrying the coarsening batch: fine SIMD layout preserved
  Coordinate fmlatt({batch,Ls,lat_size[0],lat_size[1],lat_size[2],lat_size[3]});
  Coordinate fmsimd({1,1,fsimd[0],fsimd[1],fsimd[2],fsimd[3]});
  Coordinate fmmpi ({1,1,mpi[0],mpi[1],mpi[2],mpi[3]});
  GridCartesian *FineMrhs = new GridCartesian(fmlatt,fmsimd,fmmpi);

  std::cout << GridLogMessage << "Nsimd  fine " << FGrid->Nsimd()
	    << "   coarse " << Coarse5d->Nsimd() << std::endl;

  GridParallelRNG RNG4(UGrid); RNG4.SeedFixedIntegers({1,2,3,4});
  GridParallelRNG RNG5(FGrid); RNG5.SeedFixedIntegers({5,6,7,8});

  //////////////////////////////////////////////////////////////////////
  // Gauge field
  //////////////////////////////////////////////////////////////////////
  LatticeGaugeField Umu(UGrid);
  if ( getenv("HOT_START") ) {
    std::cout << GridLogMessage << "Hot start gauge field" << std::endl;
    SU<Nc>::HotConfiguration(RNG4,Umu);
  } else {
    std::string file("/ccs/home/poare/ckpoint_lat.1000");
    if ( getenv("CONFIG") ) file = std::string(getenv("CONFIG"));
    std::cout << GridLogMessage << "Reading gauge field " << file << std::endl;
    FieldMetaData header;
    NerscIO::readConfiguration(Umu,header,file);
  }

  MobiusFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c);
  MobiusFermionD Dpv (Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,1.0, M5,b,c);

  typedef PVdagMLinearOperator<MobiusFermionD,LatticeFermionD> PVdagM_t;
  PVdagM_t PVdagM(Ddwf,Dpv);

  //////////////////////////////////////////////////////////////////////
  // Level 1 types: unvectorised coarse scalar
  //////////////////////////////////////////////////////////////////////
  typedef sTComplexD                                                          CComplexS;
  typedef MultiGeneralCoarsenedOperatorV2<vSpinColourVector,CComplexS,nbasis> CoarseOperator;
  typedef CoarseOperator::CoarseVector                                        CoarseVector;
  typedef Aggregation<vSpinColourVector,CComplexS,nbasis>                     Subspace;

  NextToNearestStencilGeometry5D geom(Coarse5d);

  //////////////////////////////////////////////////////////////////////
  // Subspace: load RAW (no Orthogonalise!), or generate.
  //
  // The Aggregation is scaffolding for CreateSubspaceGCR only. That runs
  // entirely on the fine grid and ends in GlobalOrthonormalise, which is a
  // whole-lattice Gram-Schmidt, so the coarse grid it holds is never
  // dereferenced and may be the unvectorised one.
  //////////////////////////////////////////////////////////////////////
  std::string subspace_file = "subspace_nb" + std::to_string(nbasis) + ".scidac";
  if ( getenv("SUBSPACE_FILE") ) subspace_file = std::string(getenv("SUBSPACE_FILE"));
  uint64_t file_exists=0;
  if ( UGrid->IsBoss() ){ std::ifstream f(subspace_file); file_exists=f.good()?1:0; }
  UGrid->GlobalSum(file_exists);

  const int cb=0;
  Subspace AggregatesGCR(Coarse5d,FGrid,cb);
  if ( file_exists ){
    std::cout << GridLogMessage << "*** Loading subspace from disk (kept RAW) ***" << std::endl;
    loadSubspace(AggregatesGCR.subspace, subspace_file);
  } else {
    std::cout << GridLogMessage << "*** GCR subspace generation ***" << std::endl;
    AggregatesGCR.CreateSubspaceGCR(RNG5,PVdagM,nbasis);
    saveSubspace(AggregatesGCR.subspace, subspace_file);
  }

  // RAW copy BEFORE CoarsenOperator block-orthonormalises in place.
  std::vector<LatticeFermionD> rawNull(nbasis,FGrid);
  for(int k=0;k<nbasis;k++) rawNull[k]=AggregatesGCR.subspace[k];

  //////////////////////////////////////////////////////////////////////
  // L1 coarsening. The fine operator is single RHS, so it is promoted to
  // the 6D batch grid; a natively multiRHS fine operator would substitute
  // here with no other change.
  //////////////////////////////////////////////////////////////////////
  CoarseOperator CoarseOpPV(geom,Coarse5d);
  CoarseOpPV.SetGrid(CoarseBatch);

  std::cout << GridLogMessage << "*** L1 CoarsenOperator, batch "<<batch<<" ***" << std::endl;
  if ( getenv("MRHS_COARSEN") ) {
    // Promote the single RHS operator and pack the batch: costs an
    // ExtractSlice/InsertSlice pair per rhs. Here for the A/B only; this is
    // the path a natively multiRHS fine operator would take.
    MrhsPromotedOperator<LatticeFermionD> MrhsPVdagM(PVdagM,FGrid,batch);
    CoarseOpPV.CoarsenOperator(MrhsPVdagM,FineMrhs,AggregatesGCR.subspace,Coarse5d);
  } else {
    // PVdagM is single RHS: apply it directly, batch on the coarse side.
    CoarseOpPV.CoarsenOperator(PVdagM,AggregatesGCR.subspace,Coarse5d,batch);
  }

  // Stay on the batch grid: the L2 coarsening drives this operator at the
  // batch. It is switched to the solve Nrhs once L2 is built.

  //////////////////////////////////////////////////////////////////////
  // psi_coarse = P^dag (RAW fine null) -> Galerkin images that carry the
  // near null content, and are free: A_c (P psi) = P A psi.
  //////////////////////////////////////////////////////////////////////
  MultiRHSBlockProject<LatticeFermionD> MrhsProjector;
  MrhsProjector.Allocate(nbasis,FGrid,Coarse5d);
  MrhsProjector.ImportBasis(AggregatesGCR.subspace);   // block orthonormal basis

  std::vector<CoarseVector> psi_coarse(nbasis,Coarse5d);
  MrhsProjector.blockProject(rawNull,psi_coarse);      // RAW vectors in
  rawNull.clear(); rawNull.shrink_to_fit();

  {
    RealD defect = GramDefect(psi_coarse);
    RealD leak   = std::sqrt((double)Coarse5d->gSites());
    std::cout << GridLogMessage << "GUARD: ||<psi_coarse|psi_coarse> - I||_F = " << defect
	      << "   (~0.23 good; ~sqrt(N_coarse)=" << leak << " = e_k leak)" << std::endl;
    GRID_ASSERT( defect < leak );
  }

  //////////////////////////////////////////////////////////////////////
  // Optional cross check of the coarse matrix elements against the V1
  // path, which needs a vectorised coarse space. Block Gram-Schmidt is
  // idempotent, so V1 may re-orthonormalise the same vectors in place
  // without a second copy of the subspace.
  //////////////////////////////////////////////////////////////////////
  if ( getenv("V1_CHECK") ) {

    typedef GeneralCoarsenedMatrix     <vSpinColourVector,vTComplex,nbasis> LittleDiracOperator;
    typedef MultiGeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> MrhsLittleDiracOperator;
    typedef Aggregation                <vSpinColourVector,vTComplex,nbasis> SubspaceV;

    Coordinate v5latt({1,clatt[0],clatt[1],clatt[2],clatt[3]});
    Coordinate v5simd({1,fsimd[0],fsimd[1],fsimd[2],fsimd[3]});
    GridCartesian *Coarse5dV = new GridCartesian(v5latt,v5simd,c5mpi);

    int nrhs_v1 = vComplex::Nsimd();
    Coordinate vmlatt({nrhs_v1,1,clatt[0],clatt[1],clatt[2],clatt[3]});
    Coordinate vmsimd({vComplex::Nsimd(),1,1,1,1,1});
    GridCartesian *CoarseMrhsV = new GridCartesian(vmlatt,vmsimd,cmmpi);

    NextToNearestStencilGeometry5D geomV(Coarse5dV);

    SubspaceV AggV(Coarse5dV,FGrid,cb);
    for(int k=0;k<nbasis;k++) AggV.subspace[k]=AggregatesGCR.subspace[k];

    LittleDiracOperator LittleDiracOpPV(geomV,FGrid,Coarse5dV);
    std::cout << GridLogMessage << "*** V1 CoarsenOperator (cross check) ***" << std::endl;
    LittleDiracOpPV.CoarsenOperator(PVdagM,AggV);

    MrhsLittleDiracOperator mrhsV1(geomV,CoarseMrhsV);
    mrhsV1.CopyMatrix(LittleDiracOpPV);

    // BLAS_A is written by GridtoBLAS in lSite order and both sides carry
    // the same scalar_object, so the two are directly comparable.
    typedef MrhsLittleDiracOperator::calcMatrix calcMatrix;
    int npoint = geom.npoint;
    RealD num=0.0, den=0.0;
    for(int p=0;p<npoint;p++){
      int64_t sites = mrhsV1.BLAS_A[p].size();
      GRID_ASSERT(sites == (int64_t)CoarseOpPV.BLAS_A[p].size());
      std::vector<calcMatrix> h1(sites),h2(sites);
      acceleratorCopyFromDevice(&mrhsV1.BLAS_A[p][0],    &h1[0],sites*sizeof(calcMatrix));
      acceleratorCopyFromDevice(&CoarseOpPV.BLAS_A[p][0],&h2[0],sites*sizeof(calcMatrix));
      ComplexD *w1=(ComplexD *)&h1[0];
      ComplexD *w2=(ComplexD *)&h2[0];
      int64_t words = sites*sizeof(calcMatrix)/sizeof(ComplexD);
      for(int64_t i=0;i<words;i++){
	ComplexD d=w1[i]-w2[i];
	num += real(d)*real(d)+imag(d)*imag(d);
	den += real(w1[i])*real(w1[i])+imag(w1[i])*imag(w1[i]);
      }
    }
    std::cout << GridLogMessage << "V1_CHECK: |A_V1|^2 = " << den << std::endl;
    std::cout << GridLogMessage << "V1_CHECK: |A_V1 - A_V2|^2 / |A_V1|^2 = " << num/den << std::endl;
    GRID_ASSERT( den > 0.0 );
    GRID_ASSERT( num/den < 1.0e-18 );
  }

  //////////////////////////////////////////////////////////////////////
  // STAGE TWO: L2 -> L3.
  //
  // The fine operator here is V2 at L1, which is natively multiRHS, so the
  // multiRHS driver applies with no promotion adapter: its D+1 grid IS the
  // batch grid the L1 operator is currently set to.
  //////////////////////////////////////////////////////////////////////
  Coordinate cclatt = clatt;
  Coordinate Block2({8,4,3,6});
  if ( getenv("BLOCK2") ){ GridCmdOptionIntVector(std::string(getenv("BLOCK2")),Block2); GRID_ASSERT(Block2.size()==4); }
  for(int d=0;d<4;d++){ GRID_ASSERT(clatt[d]%Block2[d]==0); cclatt[d]=clatt[d]/Block2[d]; }
  std::cout << GridLogMessage << "Block2 " << Block2 << "  coarse-coarse lattice " << cclatt << std::endl;

  Coordinate cc5latt({1,cclatt[0],cclatt[1],cclatt[2],cclatt[3]});
  GridCartesian *CoarseCoarse5d = new GridCartesian(cc5latt,c5simd,c5mpi);

  Coordinate ccmlatt({nrhs,1,cclatt[0],cclatt[1],cclatt[2],cclatt[3]});
  GridCartesian *CoarseCoarseMrhs = new GridCartesian(ccmlatt,cmsimd,cmmpi);

  Coordinate ccblatt({batch,1,cclatt[0],cclatt[1],cclatt[2],cclatt[3]});
  GridCartesian *CoarseCoarseBatch = new GridCartesian(ccblatt,cmsimd,cmmpi);

  // Coarsening deepens the tensor nest by one iScalar
  typedef CoarseVector::vector_object                                         CoarseSiteObj;
  typedef iScalar<CComplexS>                                                  CComplexS2;
  typedef MultiGeneralCoarsenedOperatorV2<CoarseSiteObj,CComplexS2,nbasis>    CoarseCoarseOperator;
  typedef CoarseCoarseOperator::CoarseVector                                  CoarseCoarseVector;

  NextToNearestStencilGeometry5D geom2(CoarseCoarse5d);

  // RAW copy of the coarse null vectors, for the same reason as rawNull:
  // the L2 CoarsenOperator block-orthonormalises its subspace in place, and
  // the L3 basis must be defined by the vectors that still carry content.
  std::vector<CoarseVector> rawPsi(nbasis,Coarse5d);
  for(int k=0;k<nbasis;k++) rawPsi[k]=psi_coarse[k];

  CoarseCoarseOperator CoarseOpL2(geom2,CoarseCoarse5d);
  CoarseOpL2.SetGrid(CoarseCoarseBatch);

  NonHermitianLinearOperator<CoarseOperator,CoarseVector> LinOpCoarse(CoarseOpPV);

  std::cout << GridLogMessage << "*** L2 CoarsenOperator, batch "<<batch<<" ***" << std::endl;
  CoarseOpL2.CoarsenOperator(LinOpCoarse,CoarseBatch,psi_coarse,CoarseCoarse5d);

  //////////////////////////////////////////////////////////////////////
  // Both operators to the solve Nrhs. The matrix elements are Nrhs
  // independent and survive the change.
  //////////////////////////////////////////////////////////////////////
  CoarseOpPV.SetGrid(CoarseMrhs);
  CoarseOpL2.SetGrid(CoarseCoarseMrhs);
  std::cout << GridLogMessage << "L1 operator at Nrhs " << CoarseOpPV.Nrhs()
	    << ",  L2 operator at Nrhs " << CoarseOpL2.Nrhs() << std::endl;

  //////////////////////////////////////////////////////////////////////
  // psi_cc from the RAW coarse null vectors, and the same guard
  //////////////////////////////////////////////////////////////////////
  MultiRHSBlockProject<CoarseVector> MrhsProjectorL2;
  MrhsProjectorL2.Allocate(nbasis,Coarse5d,CoarseCoarse5d);
  MrhsProjectorL2.ImportBasis(psi_coarse);            // block orthonormal basis

  {
    std::vector<CoarseCoarseVector> psi_cc(nbasis,CoarseCoarse5d);
    MrhsProjectorL2.blockProject(rawPsi,psi_cc);      // RAW vectors in

    RealD defect = GramDefect(psi_cc);
    RealD leak   = std::sqrt((double)CoarseCoarse5d->gSites());
    std::cout << GridLogMessage << "GUARD: ||<psi_cc|psi_cc> - I||_F = " << defect
	      << "   (~0.23 good; ~sqrt(N_cc)=" << leak << " = e_k leak)" << std::endl;
    GRID_ASSERT( defect < leak );
  }
  rawPsi.clear(); rawPsi.shrink_to_fit();

  //////////////////////////////////////////////////////////////////////
  // Both coarse operators apply on their solve grids
  //////////////////////////////////////////////////////////////////////
  {
    GridParallelRNG cRNG(Coarse5d); cRNG.SeedFixedIntegers({3,4,5,6});
    CoarseVector cin(CoarseMrhs), cout_(CoarseMrhs);
    random(cRNG,cin);
    CoarseOpPV.M(cin,cout_);
    std::cout << GridLogMessage << "L1 apply |in|^2 = " << norm2(cin)
	      << "   |M in|^2 = " << norm2(cout_) << std::endl;
    GRID_ASSERT( norm2(cout_) > 0.0 );

    GridParallelRNG ccRNG(CoarseCoarse5d); ccRNG.SeedFixedIntegers({7,8,9,10});
    CoarseCoarseVector ccin(CoarseCoarseMrhs), ccout(CoarseCoarseMrhs);
    random(ccRNG,ccin);
    CoarseOpL2.M(ccin,ccout);
    std::cout << GridLogMessage << "L2 apply |in|^2 = " << norm2(ccin)
	      << "   |M in|^2 = " << norm2(ccout) << std::endl;
    GRID_ASSERT( norm2(ccout) > 0.0 );
  }

  std::cout << GridLogMessage << "*** stage two complete: L1 and L2 coarse operators built ***" << std::endl;

  Grid_finalize();
}
