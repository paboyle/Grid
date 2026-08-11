/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Example_pvdagm_census.cc

    Copyright (C) 2026

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

// Spectral census of the coarsened PVdagM operator A_c.
//
// Measures the three sets that discriminate between the candidate explanations
// for slow coarse-grid Krylov convergence:
//
//   0. Fine Ritz diagonal of RAW subspace vectors (pre-block-orthog).
//      NB CoarsenOperator block-orthogonalises subspace[] IN PLACE; all
//      nulliness/deflation bases must be built from a raw copy.
//   1. Adjoint correctness check    <y|A x> == <A^dag y|x>      (fail-fast)
//   2. Raw-vector coarse images vs A_c: RQ (must equal CENSUS 0 by Galerkin),
//      ||A_c psi_c||/||psi_c||, and representability error
//   3. sigma_max^2 = lambda_max(A_c^dag A_c)  via power method
//   4. Low singular values          Chebyshev-filtered IRL on A_c^dag A_c
//      -> sigma_min census = pseudospectrum of A_c evaluated at the origin
//   5. Half-plane margin            lambda_min/max of H = (A_c + A_c^dag)/2
//      -> min Re W(A_c);  positive-real check (Eisenstat-Elman-Schultz bound)
//
// Interpretation:
//   sigma_min ~ min|lambda|, ~nbasis tiny then gap : effectively normal, bipartite
//   sigma_min ~ min|lambda|, dense low tail        : normal but rank-starved
//   sigma_min << min|lambda|                       : non-normal near origin
//   lambda_min(H) < 0                              : half-plane condition violated
//
// Requires the dagger code path in GeneralCoarsenedMatrix:
//   _Adag allocated, PopulateAdag active, _Adag exchanged, hermitian=0.
//
// Env vars:
//   MASS               fermion mass                     (default 0.00078)
//   SUBSPACE_FILE      subspace cache path
//   CoarseSolverShift  shift baked into coarsening      (default 0.0: pure Galerkin)
//   CENSUS_NSTOP       converged low modes wanted       (default 60)
//   CENSUS_NK          Lanczos Nk                       (default 96)
//   CENSUS_NM          Lanczos Nm                       (default 192)
//   CENSUS_TOL         Lanczos residual                 (default 1e-5)
//   CENSUS_MAXIT       Lanczos max restarts             (default 50)
//   CHEBY_LO           filter low edge in sigma^2       (default 4.0)
//   CHEBY_HI           filter high edge; 0 = auto from power method x1.1
//   CHEBY_ORDER        filter order                     (default 401)
//                      filter gain at 0 ~ cosh(order*2*sqrt(lo/hi)); with
//                      hi~2200, lo=4, order=401 => gain ~ 1e14. lo=0.01 at
//                      order 201 gives gain ~1.4 (stagnation).

#include <Grid/Grid.h>
#include <Grid/Grid_Eigen_Dense.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>

using namespace std;
using namespace Grid;

RealD mass              = 0.00078;
RealD CoarseSolverShift = 0.0;
int   CensusNstop       = 60;
int   CensusNk          = 96;
int   CensusNm          = 192;
RealD CensusTol         = 1.0e-5;
int   CensusMaxIt       = 50;
RealD ChebyLo           = 4.0;   // sigma^2 cutoff: amplifies sigma < 2. Filter gain ~ cosh(order*2*sqrt(lo/hi))
RealD ChebyHi           = 0.0;   // 0 => auto: 1.1 * power-method sigma_max^2
int   ChebyOrder        = 401;
RealD CGdeflTol         = 1.0e-8; // CENSUS 6 deflated-CG tolerance
int   CGdeflMaxIt       = 4000;   // CENSUS 6 deflated-CG max iterations
int   DeflRank          = 0;      // CENSUS 6 deflation rank; 0 => all available per basis

void ParseEnvironment(void)
{
  if(getenv("MASS"))              mass              = atof(getenv("MASS"));
  if(getenv("CoarseSolverShift")) CoarseSolverShift = atof(getenv("CoarseSolverShift"));
  if(getenv("CENSUS_NSTOP"))      CensusNstop       = atoi(getenv("CENSUS_NSTOP"));
  if(getenv("CENSUS_NK"))         CensusNk          = atoi(getenv("CENSUS_NK"));
  if(getenv("CENSUS_NM"))         CensusNm          = atoi(getenv("CENSUS_NM"));
  if(getenv("CENSUS_TOL"))        CensusTol         = atof(getenv("CENSUS_TOL"));
  if(getenv("CENSUS_MAXIT"))      CensusMaxIt       = atoi(getenv("CENSUS_MAXIT"));
  if(getenv("CHEBY_LO"))          ChebyLo           = atof(getenv("CHEBY_LO"));
  if(getenv("CHEBY_HI"))          ChebyHi           = atof(getenv("CHEBY_HI"));
  if(getenv("CHEBY_ORDER"))       ChebyOrder        = atoi(getenv("CHEBY_ORDER"));
  if(getenv("CGDEFL_TOL"))        CGdeflTol         = atof(getenv("CGDEFL_TOL"));
  if(getenv("CGDEFL_MAXIT"))      CGdeflMaxIt       = atoi(getenv("CGDEFL_MAXIT"));
  if(getenv("DEFL_RANK"))         DeflRank          = atoi(getenv("DEFL_RANK"));

  std::cout << GridLogMessage << "PARAM: MASS              " << mass              << std::endl;
  std::cout << GridLogMessage << "PARAM: CoarseSolverShift " << CoarseSolverShift << std::endl;
  std::cout << GridLogMessage << "PARAM: CENSUS_NSTOP      " << CensusNstop       << std::endl;
  std::cout << GridLogMessage << "PARAM: CENSUS_NK         " << CensusNk          << std::endl;
  std::cout << GridLogMessage << "PARAM: CENSUS_NM         " << CensusNm          << std::endl;
  std::cout << GridLogMessage << "PARAM: CENSUS_TOL        " << CensusTol         << std::endl;
  std::cout << GridLogMessage << "PARAM: CENSUS_MAXIT      " << CensusMaxIt       << std::endl;
  std::cout << GridLogMessage << "PARAM: CHEBY_LO          " << ChebyLo           << std::endl;
  std::cout << GridLogMessage << "PARAM: CHEBY_HI          " << ChebyHi           << std::endl;
  std::cout << GridLogMessage << "PARAM: CHEBY_ORDER       " << ChebyOrder        << std::endl;
}

template <class Field>
void saveSubspace(std::vector<Field> &subspace, std::string const fname){
#ifdef HAVE_LIME
  std::cout << Grid::GridLogMessage << "Saving subspace (" << subspace.size() << " vectors) to: " << fname << std::endl;
  Grid::emptyUserRecord record;
  Grid::ScidacWriter SW(subspace[0].Grid()->IsBoss());
  SW.open(fname);
  for (int k = 0; k < (int)subspace.size(); k++)
    SW.writeScidacFieldRecord(subspace[k], record);
  SW.close();
#endif
}

template <class Field>
void loadSubspace(std::vector<Field> &subspace, std::string const fname){
#ifdef HAVE_LIME
  std::cout << Grid::GridLogMessage << "Loading subspace (" << subspace.size() << " vectors) from: " << fname << std::endl;
  Grid::emptyUserRecord record;
  Grid::ScidacReader SR;
  SR.open(fname);
  for (int k = 0; k < (int)subspace.size(); k++)
    SR.readScidacFieldRecord(subspace[k], record);
  SR.close();
#endif
}

template<class Matrix,class Field>
class PVdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
  Matrix &_PV;
public:
  PVdagMLinearOperator(Matrix &Mat,Matrix &PV): _Mat(Mat),_PV(PV) {};
  void OpDiag (const Field &in, Field &out) {    assert(0);  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {    assert(0);  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){    assert(0);  };
  void Op     (const Field &in, Field &out){
    Field tmp(in.Grid());
    _Mat.M(in,tmp);
    _PV.Mdag(tmp,out);
  }
  void AdjOp     (const Field &in, Field &out){
    Field tmp(in.Grid());
    _PV.M(in,tmp);
    _Mat.Mdag(tmp,out);
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){
    HermOp(in,out);
    ComplexD dot = innerProduct(in,out);
    n1=real(dot);
    n2=norm2(out);
  }
  void HermOp(const Field &in, Field &out){
    Field tmp(in.Grid());
    Op(in,tmp);
    AdjOp(tmp,out);
  }
};

template<class Matrix,class Field>
class ShiftedPVdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
  Matrix &_PV;
public:
  RealD shift;
  ShiftedPVdagMLinearOperator(RealD _shift,Matrix &Mat,Matrix &PV): shift(_shift),_Mat(Mat),_PV(PV){};
  void OpDiag (const Field &in, Field &out) {    assert(0);  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {    assert(0);  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){    assert(0);  };
  void Op     (const Field &in, Field &out){
    Field tmp(in.Grid());
    _Mat.M(in,tmp);
    _PV.Mdag(tmp,out);
    out = out + shift * in;
  }
  void AdjOp     (const Field &in, Field &out){
    Field tmp(in.Grid());
    _PV.M(tmp,out);
    _Mat.Mdag(in,tmp);
    out = out + shift * in;
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){    assert(0);  }
  void HermOp(const Field &in, Field &out){
    Field tmp(in.Grid());
    Op(in,tmp);
    AdjOp(tmp,out);
  }
};

// H = (A + A^dag)/2 : Hermitian part of the coarse operator.
// lambda_min(H) = min Re W(A) is the half-plane margin; the EES GCR
// convergence theorem requires it positive.
template<class Matrix,class Field>
class HermitianPartOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
public:
  HermitianPartOperator(Matrix &Mat): _Mat(Mat) {};
  void OpDiag (const Field &in, Field &out) {    assert(0);  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {    assert(0);  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){    assert(0);  };
  void Op     (const Field &in, Field &out){ HermOp(in,out); }
  void AdjOp  (const Field &in, Field &out){ HermOp(in,out); }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){
    HermOp(in,out);
    ComplexD dot = innerProduct(in,out);
    n1=real(dot);
    n2=norm2(out);
  }
  void HermOp(const Field &in, Field &out){
    Field tmp(in.Grid());
    _Mat.M(in,out);
    _Mat.Mdag(in,tmp);
    out = 0.5*(out + tmp);
  }
};

// s*I - Op : power method on this gives s - lambda_min(Op) for Hermitian Op.
template<class Field>
class ShiftedNegatedOperator : public LinearOperatorBase<Field> {
  LinearOperatorBase<Field> &_Op;
  RealD s;
public:
  ShiftedNegatedOperator(RealD _s, LinearOperatorBase<Field> &Op): _Op(Op), s(_s) {};
  void OpDiag (const Field &in, Field &out) {    assert(0);  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {    assert(0);  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){    assert(0);  };
  void Op     (const Field &in, Field &out){ HermOp(in,out); }
  void AdjOp  (const Field &in, Field &out){ HermOp(in,out); }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){
    HermOp(in,out);
    ComplexD dot = innerProduct(in,out);
    n1=real(dot);
    n2=norm2(out);
  }
  void HermOp(const Field &in, Field &out){
    _Op.HermOp(in,out);
    out = s*in - out;
  }
};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);
  ParseEnvironment();

  const int Ls=24;
  RealD M5=1.8;
  RealD b=1.5;
  RealD c=0.5;
  const int nbasis = 60;

  std::cout << GridLogMessage << "Census of coarse PVdagM: mass=" << mass << " Ls=" << Ls << " nbasis=" << nbasis << std::endl;

  std::vector<int> lat_size {48, 48, 48, 96};

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(lat_size, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  // Blocking: default matches Example_pvdagm.cc; override with e.g. BLOCK=2.2.2.2
  Coordinate clatt = lat_size;
  Coordinate Block({4,4,6,4});
  if ( getenv("BLOCK") ) {
    GridCmdOptionIntVector(std::string(getenv("BLOCK")),Block);
    GRID_ASSERT(Block.size()==4);
  }
  for(int d=0;d<clatt.size();d++){
    GRID_ASSERT(lat_size[d] % Block[d] == 0);
    clatt[d] = lat_size[d]/Block[d];
  }
  std::cout << GridLogMessage << "Block " << Block << "  coarse lattice " << clatt << std::endl;

  GridCartesian *Coarse4d =  SpaceTimeGrid::makeFourDimGrid(clatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridCartesian *Coarse5d =  SpaceTimeGrid::makeFiveDimGrid(1,Coarse4d);

  GridParallelRNG RNG5(FGrid);    RNG5.SeedFixedIntegers({5,6,7,8});
  GridParallelRNG RNG4(UGrid);    RNG4.SeedFixedIntegers({1,2,3,4});
  GridParallelRNG CRNG(Coarse5d); CRNG.SeedFixedIntegers({5,6,7,8});

  LatticeGaugeField Umu(UGrid);
  std::cout << GridLogMessage << "Reading gauge field" << std::endl;
  FieldMetaData header;
  std::string file("/ccs/home/poare/ckpoint_lat.1000");
  NerscIO::readConfiguration(Umu,header,file);

  MobiusFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c);
  MobiusFermionD Dpv (Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,1.0, M5,b,c);

  typedef PVdagMLinearOperator<MobiusFermionD,LatticeFermionD>        PVdagM_t;
  typedef ShiftedPVdagMLinearOperator<MobiusFermionD,LatticeFermionD> ShiftedPVdagM_t;
  typedef GeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis>  LittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector                           CoarseVector;
  typedef Aggregation<vSpinColourVector,vTComplex,nbasis>             Subspace;

  PVdagM_t        PVdagM(Ddwf,Dpv);
  ShiftedPVdagM_t ShiftedPVdagM(CoarseSolverShift,Ddwf,Dpv);

  NextToNearestStencilGeometry5D geom(Coarse5d);

  //////////////////////////////////////////////////////////////////////
  // Subspace: load from cache or generate
  //////////////////////////////////////////////////////////////////////
  std::string subspace_file = "/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/subspace_nb"
                              + std::to_string(nbasis) + ".scidac";
  if ( getenv("SUBSPACE_FILE") ) subspace_file = std::string(getenv("SUBSPACE_FILE"));

  uint64_t file_exists = 0;
  if ( UGrid->IsBoss() ) {
    std::ifstream f(subspace_file);
    file_exists = f.good() ? 1 : 0;
  }
  UGrid->GlobalSum(file_exists);

  const int cb = 0;
  Subspace AggregatesGCR(Coarse5d,FGrid,cb);

  if ( file_exists ) {
    std::cout << GridLogMessage << "*** Loading subspace from disk ***" << std::endl;
    loadSubspace(AggregatesGCR.subspace, subspace_file);
  } else {
    std::cout << GridLogMessage << "*** GCR subspace generation ***" << std::endl;
    AggregatesGCR.CreateSubspaceGCR(RNG5,PVdagM,nbasis);
    saveSubspace(AggregatesGCR.subspace, subspace_file);
  }

  //////////////////////////////////////////////////////////////////////
  // Keep the RAW (pre-block-orthogonalisation) near-null vectors.
  // CoarsenOperator block-orthogonalises subspace[] IN PLACE, after which
  // subspace[k] is the orthonormal basis phi_k and Project(phi_k) = e_k,
  // the block-constant unit vector -- NOT a near-null direction.
  // All nulliness measurements and any deflation basis must use raw[].
  //////////////////////////////////////////////////////////////////////
  std::vector<LatticeFermionD> raw(nbasis,FGrid);
  for(int k=0;k<nbasis;k++) raw[k] = AggregatesGCR.subspace[k];

  //////////////////////////////////////////////////////////////////////
  // CENSUS 0: fine-grid Ritz diagonal on the loaded/generated raw vectors.
  // Expect Re <psi|A|psi>/<psi|psi> ~ the nulliness achieved at generation
  // (~2e-3). O(0.1-10) values mean the cache holds orthogonalised vectors
  // and must be regenerated.
  //////////////////////////////////////////////////////////////////////
  std::cout << GridLogMessage << "=================================================" << std::endl;
  std::cout << GridLogMessage << "CENSUS 0: fine Ritz diagonal of raw subspace vectors" << std::endl;
  std::cout << GridLogMessage << "=================================================" << std::endl;
  {
    LatticeFermionD Ap(FGrid);
    for(int k=0;k<nbasis;k++){
      PVdagM.Op(raw[k],Ap);
      RealD    n2psi = norm2(raw[k]);
      ComplexD rq    = innerProduct(raw[k],Ap)/n2psi;
      std::cout << GridLogMessage << "CENSUS: raw[" << k << "]  <psi|A|psi>/<psi|psi> = " << rq
                << "   ||A psi||/||psi|| = " << std::sqrt(norm2(Ap)/n2psi) << std::endl;
    }
  }

  //////////////////////////////////////////////////////////////////////
  // Coarsen. hermitian=0 is REQUIRED: enables PopulateAdag so that
  // Mdag applies A^dag rather than silently aliasing to A.
  //////////////////////////////////////////////////////////////////////
  LittleDiracOperator LittleDiracOpPV(geom,FGrid,Coarse5d,0);
  if ( CoarseSolverShift != 0.0 ) {
    std::cout << GridLogMessage << "Coarsening SHIFTED operator, shift=" << CoarseSolverShift << std::endl;
    LittleDiracOpPV.CoarsenOperator(ShiftedPVdagM, AggregatesGCR);
  } else {
    std::cout << GridLogMessage << "Coarsening pure Galerkin operator (no shift)" << std::endl;
    LittleDiracOpPV.CoarsenOperator(PVdagM, AggregatesGCR);
  }

  CoarseVector c_x(Coarse5d);
  CoarseVector c_y(Coarse5d);
  CoarseVector c_t1(Coarse5d);
  CoarseVector c_t2(Coarse5d);

  //////////////////////////////////////////////////////////////////////
  // CENSUS 1: adjoint correctness (fail fast)
  //   <y|A x> == <A^dag y|x>  for random x,y
  //////////////////////////////////////////////////////////////////////
  std::cout << GridLogMessage << "=================================================" << std::endl;
  std::cout << GridLogMessage << "CENSUS 1: adjoint correctness of dagger code path" << std::endl;
  std::cout << GridLogMessage << "=================================================" << std::endl;
  random(CRNG,c_x);
  random(CRNG,c_y);
  LittleDiracOpPV.M   (c_x,c_t1);       // A x
  LittleDiracOpPV.Mdag(c_y,c_t2);       // A^dag y
  ComplexD ip1 = innerProduct(c_y,c_t1);   // <y|A x>
  ComplexD ip2 = innerProduct(c_t2,c_x);   // <A^dag y|x>
  RealD reldiff = abs(ip1-ip2)/abs(ip1);
  std::cout << GridLogMessage << "CENSUS: <y|Ax>       = " << ip1 << std::endl;
  std::cout << GridLogMessage << "CENSUS: <Adag y|x>   = " << ip2 << std::endl;
  std::cout << GridLogMessage << "CENSUS: rel diff     = " << reldiff << "  (expect ~1e-14; FAIL if O(1))" << std::endl;
  GRID_ASSERT(reldiff < 1.0e-8);

  // Coarse near-null ("global") vectors psi_c[k] = P^dag raw[k], stored for the
  // Ritz-matrix + deflation study in CENSUS 6 (filled in CENSUS 2's projection
  // loop below, before raw[]/subspace[] are freed).
  std::vector<CoarseVector> psi_c(nbasis,Coarse5d);

  //////////////////////////////////////////////////////////////////////
  // CENSUS 2: nulliness of the RAW vectors' coarse images against A_c.
  //   psi_c[k] = P^dag raw[k].  Galerkin guarantees the Rayleigh quotient
  //   equals CENSUS 0's fine value exactly (raw[k] is in span of its own
  //   chopped pieces) -- agreement is a machine-precision validation of
  //   the coarsening.  ||A_c psi_c||/||psi_c|| is the sigma-relevant norm.
  //   The representability column ||raw - P psi_c||/||raw|| must be ~eps.
  //////////////////////////////////////////////////////////////////////
  std::cout << GridLogMessage << "=================================================" << std::endl;
  std::cout << GridLogMessage << "CENSUS 2: raw-vector coarse images against coarse operator" << std::endl;
  std::cout << GridLogMessage << "=================================================" << std::endl;
  {
    LatticeFermionD back(FGrid);
    for(int k=0;k<nbasis;k++){
      AggregatesGCR.ProjectToSubspace(c_x, raw[k]);
      psi_c[k] = c_x;                       // store coarse near-null vector for CENSUS 6
      AggregatesGCR.PromoteFromSubspace(c_x, back);
      back = back - raw[k];
      RealD represent = std::sqrt(norm2(back)/norm2(raw[k]));
      LittleDiracOpPV.M(c_x, c_t1);
      RealD    n2psi = norm2(c_x);
      RealD    n2Apsi= norm2(c_t1);
      ComplexD rq    = innerProduct(c_x,c_t1) / n2psi;
      std::cout << GridLogMessage << "CENSUS: psi_c[" << k << "]  <psi|A|psi>/<psi|psi> = " << rq
                << "   ||A psi||/||psi|| = " << std::sqrt(n2Apsi/n2psi)
                << "   represent_err = " << represent << std::endl;
    }
  }

  // Fine subspace + raw copy are needed only through CENSUS 2; CENSUS 3-5 are
  // entirely coarse (LittleDiracOpPV only), and the CENSUS 4 evec save writes the
  // coarse vectors directly.  Release the ~2*nbasis fine 5D fields (~14 GB/GCD at
  // 2^4) HERE, before the order-ChebyOrder Lanczos whose padded coarse temporaries
  // otherwise push host memory over the top on top of _A + _Adag (the AccCache
  // CpuPtr!=NULL abort seen mid-iteration).
  // Direct orthonormality check of the fine near-null vectors (GlobalOrthonormalise
  // in CreateSubspaceGCR).  raw is freed just below, so this runs here, not CENSUS 6.
  // If this is ~0 but the coarse Gram S (CENSUS 6) is not, the gap is representability,
  // not orthonormality.
  {
    Eigen::MatrixXcd Gfine(nbasis,nbasis);
    for(int i=0;i<nbasis;i++){
      for(int j=i;j<nbasis;j++){
        ComplexD g = innerProduct(raw[i],raw[j]);
        Gfine(i,j) = std::complex<double>(g.real(),g.imag());
        Gfine(j,i) = std::conj(Gfine(i,j));
      }
    }
    double GmI = (Gfine - Eigen::MatrixXcd::Identity(nbasis,nbasis)).norm();
    std::cout << GridLogMessage << "CENSUS 2b: fine Gram ||<raw_i|raw_j> - I||_F = " << GmI
              << "   (expect ~0 if fine vectors orthonormal)" << std::endl;
  }

  raw.clear();                    raw.shrink_to_fit();
  AggregatesGCR.subspace.clear(); AggregatesGCR.subspace.shrink_to_fit();

  //////////////////////////////////////////////////////////////////////
  // CENSUS 3: sigma_max^2 = lambda_max( A_c^dag A_c ) by power method
  //////////////////////////////////////////////////////////////////////
  std::cout << GridLogMessage << "=================================================" << std::endl;
  std::cout << GridLogMessage << "CENSUS 3: power method for sigma_max" << std::endl;
  std::cout << GridLogMessage << "=================================================" << std::endl;
  MdagMLinearOperator<LittleDiracOperator,CoarseVector> HermOpAdagA(LittleDiracOpPV);
  random(CRNG,c_x);
  PowerMethod<CoarseVector> PM;
  RealD sigmax2 = PM(HermOpAdagA,c_x);
  std::cout << GridLogMessage << "CENSUS: lambda_max(AdagA) = " << sigmax2
            << "   sigma_max = " << std::sqrt(sigmax2) << std::endl;

  //////////////////////////////////////////////////////////////////////
  // CENSUS 4: low singular values via Chebyshev-filtered IRL on A^dag A
  //   The low end of sigma(A_c) is the pseudospectrum of A_c at z=0.
  //////////////////////////////////////////////////////////////////////
  std::cout << GridLogMessage << "=================================================" << std::endl;
  std::cout << GridLogMessage << "CENSUS 4: Chebyshev-filtered Lanczos, low sigma^2" << std::endl;
  std::cout << GridLogMessage << "=================================================" << std::endl;
  RealD cheby_hi = (ChebyHi > 0.0) ? ChebyHi : 1.1*sigmax2;
  std::cout << GridLogMessage << "Chebyshev filter [" << ChebyLo << "," << cheby_hi << "] order " << ChebyOrder << std::endl;

  // eval/evec/Nconv hoisted out of the block so CENSUS 6 can deflate with them.
  std::vector<RealD>        eval(CensusNm);
  std::vector<CoarseVector> evec(CensusNm,Coarse5d);
  int Nconv=0;
  {
    Chebyshev<CoarseVector>      Cheby(ChebyLo,cheby_hi,ChebyOrder);
    FunctionHermOp<CoarseVector> OpCheby(Cheby,HermOpAdagA);
    PlainHermOp<CoarseVector>    Op     (HermOpAdagA);

    ImplicitlyRestartedLanczos<CoarseVector> IRL(OpCheby,Op,CensusNstop,CensusNk,CensusNm,CensusTol,CensusMaxIt);

    random(CRNG,c_x);
    IRL.calc(eval,evec,c_x,Nconv);

    std::cout << GridLogMessage << "CENSUS: converged " << Nconv << " modes of AdagA" << std::endl;
    for(int i=0;i<Nconv;i++){
      std::cout << GridLogMessage << "CENSUS: sigma[" << i << "]^2 = " << eval[i]
                << "   sigma = " << std::sqrt(std::max(eval[i],0.0)) << std::endl;
    }

    // Optionally persist the low right-singular-vector basis: this IS the
    // deflation basis for the coarse solve (ADEF1 / MultiRHSDeflation).
    // Set CENSUS_EVEC_FILE to enable.
    if ( getenv("CENSUS_EVEC_FILE") && Nconv>0 ) {
#ifdef HAVE_LIME
      std::string evec_file(getenv("CENSUS_EVEC_FILE"));
      std::string eval_file = evec_file + ".evals.xml";
      std::cout << GridLogMessage << "CENSUS: saving " << Nconv << " singular vectors to " << evec_file << std::endl;
      emptyUserRecord record;
      ScidacWriter WR(evec[0].Grid()->IsBoss());
      WR.open(evec_file);
      for(int i=0;i<Nconv;i++) WR.writeScidacFieldRecord(evec[i],record);
      WR.close();
      XmlWriter WRx(eval_file);
      std::vector<RealD> eval_out(eval.begin(),eval.begin()+Nconv); // don't shrink shared eval
      write(WRx,"evals",eval_out);
#endif
    }
  }

  // NB: evec/eval stay sized CensusNm (Lattice has no default ctor, so
  // std::vector<CoarseVector>::resize won't instantiate).  They match in size,
  // which is all DeflatedGuesser asserts; CENSUS 6 only ever indexes [0,Nconv).

  //////////////////////////////////////////////////////////////////////
  // CENSUS 5: half-plane margin from the Hermitian part
  //   lambda_min(H) = min Re W(A_c) > 0 <=> positive-real (EES applies)
  //////////////////////////////////////////////////////////////////////
  std::cout << GridLogMessage << "=================================================" << std::endl;
  std::cout << GridLogMessage << "CENSUS 5: Hermitian part H=(A+Adag)/2, half-plane margin" << std::endl;
  std::cout << GridLogMessage << "=================================================" << std::endl;
  HermitianPartOperator<LittleDiracOperator,CoarseVector> HermPart(LittleDiracOpPV);

  random(CRNG,c_x);
  RealD lamHmax = PM(HermPart,c_x);
  std::cout << GridLogMessage << "CENSUS: lambda_max(H) = " << lamHmax << std::endl;

  // lambda_min(H): most-negative eigenvalue via Chebyshev-filtered IRL on H.
  // A shifted power method cannot separate it from the dense low tail (which is
  // why the earlier -0.006 is suspect); Cheby(lo, hi>=lambda_max) amplifies the
  // most-negative mode hardest so IRL isolates the true bottom of the spectrum.
  RealD hpLo    = getenv("HPLANE_CHEBY_LO")    ? atof(getenv("HPLANE_CHEBY_LO"))    : 0.1;
  RealD hpHi    = getenv("HPLANE_CHEBY_HI")    ? atof(getenv("HPLANE_CHEBY_HI"))    : 1.1*lamHmax;
  int   hpOrder = getenv("HPLANE_CHEBY_ORDER") ? atoi(getenv("HPLANE_CHEBY_ORDER")) : 61;
  // Grid's Chebyshev filter MUST be odd order (positive for x < -1, where the low/
  // negative modes map); an even order flips the sign there and the IRL blows up.
  if(hpOrder%2==0){ hpOrder++;
    std::cout<<GridLogMessage<<"HPLANE_CHEBY_ORDER forced odd -> "<<hpOrder<<std::endl; }
  int   hpNstop = getenv("HPLANE_NSTOP")       ? atoi(getenv("HPLANE_NSTOP"))       : 8;
  int   hpNk    = getenv("HPLANE_NK")          ? atoi(getenv("HPLANE_NK"))          : 24;
  int   hpNm    = getenv("HPLANE_NM")          ? atoi(getenv("HPLANE_NM"))          : 48;
  RealD hpTol   = getenv("HPLANE_TOL")         ? atof(getenv("HPLANE_TOL"))         : 1.0e-4;
  int   hpMaxIt = getenv("HPLANE_MAXIT")       ? atoi(getenv("HPLANE_MAXIT"))       : 20;

  Chebyshev<CoarseVector>      HCheby(hpLo,hpHi,hpOrder);
  FunctionHermOp<CoarseVector> HOpCheby(HCheby,HermPart);
  PlainHermOp<CoarseVector>    HOpPlain(HermPart);
  ImplicitlyRestartedLanczos<CoarseVector> HIRL(HOpCheby,HOpPlain,hpNstop,hpNk,hpNm,hpTol,hpMaxIt);
  std::vector<RealD>        heval(hpNm);
  std::vector<CoarseVector> hevec(hpNm,Coarse5d);
  int hNconv=0;
  random(CRNG,c_x);
  HIRL.calc(heval,hevec,c_x,hNconv);
  RealD lamHmin = (hNconv>0) ? heval[0] : 9.99e99;
  for(int kk=0;kk<hNconv;kk++) lamHmin = std::min(lamHmin, heval[kk]);
  std::cout << GridLogMessage << "CENSUS: IRL H-bottom converged " << hNconv
            << " eigenvalues; most-negative = " << lamHmin << std::endl;
  std::cout << GridLogMessage << "CENSUS: lambda_min(H) = " << lamHmin
            << "  (positive-real / half-plane margin; NEGATIVE => GCR unguaranteed)" << std::endl;

  //////////////////////////////////////////////////////////////////////
  // CENSUS 6: Ritz matrix of the coarse near-null basis + deflated-CG study
  //
  //   C_ij = <psi_c^i | A^dag A | psi_c^j>,   S_ij = <psi_c^i | psi_c^j>.
  //   psi_c are NOT orthonormal (raw near-null projected to coarse), so the
  //   Rayleigh-Ritz problem is the GENERALISED Hermitian one  C v = theta S v.
  //   Its eigenpairs (theta_i, g_i = sum_j V(j,i) psi_c^j) are the best approximate
  //   eigenpairs of A^dag A available from span{psi_c}; Eigen normalises so that
  //   V^dag S V = I, hence <g_i|g_j> = delta_ij and the g_i are an orthonormal
  //   DeflatedGuesser basis.  Compare theta_i to the Lanczos sigma_i^2, then run
  //   three CG solves on A^dag A: [1] no deflation, [2] Lanczos-eigenvector
  //   deflated guess, [3] Ritz global-vector deflated guess (g_i treated as pure
  //   eigenvectors with eigenvalue theta_i).
  //////////////////////////////////////////////////////////////////////
  std::cout << GridLogMessage << "=================================================" << std::endl;
  std::cout << GridLogMessage << "CENSUS 6: Ritz matrix C_ij = <psi_c^i|AdagA|psi_c^j> + deflated CG" << std::endl;
  std::cout << GridLogMessage << "=================================================" << std::endl;

  std::vector<CoarseVector> Apsi(nbasis,Coarse5d);
  for(int j=0;j<nbasis;j++) HermOpAdagA.HermOp(psi_c[j],Apsi[j]);   // A^dag A psi_c^j

  Eigen::MatrixXcd Cmat(nbasis,nbasis);
  Eigen::MatrixXcd Smat(nbasis,nbasis);
  for(int i=0;i<nbasis;i++){
    for(int j=0;j<nbasis;j++){
      ComplexD cij = innerProduct(psi_c[i],Apsi[j]);
      ComplexD sij = innerProduct(psi_c[i],psi_c[j]);
      Cmat(i,j) = std::complex<double>(cij.real(),cij.imag());
      Smat(i,j) = std::complex<double>(sij.real(),sij.imag());
    }
  }

  {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ses(Smat);
    Eigen::MatrixXcd Id = Eigen::MatrixXcd::Identity(nbasis,nbasis);
    double SmI = (Smat - Id).norm();   // ||S - I||_F : ~0 iff psi_c orthonormal
    std::cout << GridLogMessage << "CENSUS 6: Gram S eig range [" << ses.eigenvalues()(0)
              << ", " << ses.eigenvalues()(nbasis-1)
              << "]   ||S - I||_F = " << SmI
              << "   (expect ~0: fine vectors are GlobalOrthonormalise'd => psi_c orthonormal)" << std::endl;
  }
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> ges(Cmat,Smat);
  Eigen::VectorXd  theta = ges.eigenvalues();    // ascending, real
  Eigen::MatrixXcd Vr    = ges.eigenvectors();   // columns; V^dag S V = I

  int ncmp = std::min((int)nbasis,Nconv);
  std::cout << GridLogMessage << "CENSUS 6: Ritz theta vs Lanczos sigma^2 (both ascending):" << std::endl;
  for(int i=0;i<nbasis;i++){
    if(i<ncmp)
      std::cout << GridLogMessage << "CENSUS: theta["<<i<<"] = "<<theta(i)
                <<"   sqrt = "<<std::sqrt(std::max(theta(i),0.0))
                <<"   | sigma^2 = "<<eval[i]<<"   theta/sigma^2 = "<<theta(i)/eval[i] << std::endl;
    else
      std::cout << GridLogMessage << "CENSUS: theta["<<i<<"] = "<<theta(i)
                <<"   sqrt = "<<std::sqrt(std::max(theta(i),0.0)) << std::endl;
  }

  // Ritz global vectors g_i = sum_j V(j,i) psi_c^j  (S-orthonormal), eigenvalue theta_i
  std::vector<CoarseVector> gvec(nbasis,Coarse5d);
  std::vector<RealD>        gval(nbasis);
  for(int i=0;i<nbasis;i++){
    gvec[i] = Zero();
    for(int j=0;j<nbasis;j++){
      ComplexD coeff(Vr(j,i).real(),Vr(j,i).imag());
      axpy(gvec[i],coeff,psi_c[j],gvec[i]);
    }
    gval[i] = theta(i);
  }

  // How good are the diagonalised global vectors as actual eigenvectors of A^dag A?
  {
    CoarseVector Ag(Coarse5d), rr(Coarse5d);
    int nchk = std::min((int)nbasis,16);
    for(int i=0;i<nchk;i++){
      HermOpAdagA.HermOp(gvec[i],Ag);
      axpy(rr,-gval[i],gvec[i],Ag);          // rr = A^dag A g - theta g
      RealD rn = std::sqrt(norm2(rr));
      std::cout << GridLogMessage << "CENSUS 6: Ritz resid ["<<i<<"] ||AdagA g - theta g||/theta = "
                << rn/std::max(gval[i],1.0e-30) << "   (theta="<<gval[i]<<")" << std::endl;
    }
  }

  // --- Three CG solves on A^dag A, common random source ---
  int rankLanc = (DeflRank>0) ? std::min(DeflRank,Nconv)       : Nconv;
  int rankRitz = (DeflRank>0) ? std::min(DeflRank,(int)nbasis) : (int)nbasis;
  std::cout << GridLogMessage << "CENSUS 6: CG tol "<<CGdeflTol<<" maxit "<<CGdeflMaxIt
            << " ; deflation ranks -- Lanczos "<<rankLanc<<", Ritz "<<rankRitz << std::endl;

  CoarseVector cg_src(Coarse5d);  random(CRNG,cg_src);
  CoarseVector cg_x  (Coarse5d);
  ConjugateGradient<CoarseVector> CGdefl(CGdeflTol,CGdeflMaxIt,false);

  cg_x = Zero();
  CGdefl(HermOpAdagA,cg_src,cg_x);
  std::cout << GridLogMessage << "CENSUS 6: [1] no deflation          : iters = "
            << CGdefl.IterationsToComplete << "   true_resid = " << CGdefl.TrueResidual << std::endl;

  if(rankLanc>0){
    DeflatedGuesser<CoarseVector> guessL(evec,eval,rankLanc);
    guessL(cg_src,cg_x);
    CGdefl(HermOpAdagA,cg_src,cg_x);
    std::cout << GridLogMessage << "CENSUS 6: [2] Lanczos-evec deflation : iters = "
              << CGdefl.IterationsToComplete << "   true_resid = " << CGdefl.TrueResidual << std::endl;
  }

  {
    DeflatedGuesser<CoarseVector> guessR(gvec,gval,rankRitz);
    guessR(cg_src,cg_x);
    CGdefl(HermOpAdagA,cg_src,cg_x);
    std::cout << GridLogMessage << "CENSUS 6: [3] Ritz-vector deflation  : iters = "
              << CGdefl.IterationsToComplete << "   true_resid = " << CGdefl.TrueResidual << std::endl;
  }

  //////////////////////////////////////////////////////////////////////
  // Summary
  //////////////////////////////////////////////////////////////////////
  std::cout << GridLogMessage << "=================================================" << std::endl;
  std::cout << GridLogMessage << "CENSUS SUMMARY" << std::endl;
  std::cout << GridLogMessage << "  sigma_max        = " << std::sqrt(sigmax2) << std::endl;
  std::cout << GridLogMessage << "  lambda_max(H)    = " << lamHmax  << std::endl;
  std::cout << GridLogMessage << "  lambda_min(H)    = " << lamHmin  << std::endl;
  std::cout << GridLogMessage << "  low sigma census : see CENSUS 4 table above" << std::endl;
  std::cout << GridLogMessage << "  Compare min sigma with |lambda| from Krylov-Schur (Patrick):" << std::endl;
  std::cout << GridLogMessage << "    sigma_min ~ min|lambda| : effectively normal; deflation rank is the issue" << std::endl;
  std::cout << GridLogMessage << "    sigma_min << min|lambda|: non-normal; need two-sided/singular-vector deflation" << std::endl;
  std::cout << GridLogMessage << "=================================================" << std::endl;

  std::cout << GridLogMessage << "Done" << std::endl;
  Grid_finalize();
  return 0;
}
