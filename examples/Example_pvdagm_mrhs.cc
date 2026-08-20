/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Example_pvdagm_mrhs.cc

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

// MultiRHS (valence) two-level multigrid for PVdagM.
//
// Philosophy: NO block Krylov, NO per-RHS lockstep coefficients. A single
// GCR polynomial for the enlarged block-diagonal system diag(A,...,A):
// all inner products are summed over the RHS index, giving one alpha/beta
// per step shared by every RHS.
//
// Level structure mirrors Example_pvdagm:
//   outer:  MrhsPGCRNonHermitian on PVdagM over std::vector<LatticeFermionD>
//   precon: V-cycle -- per-RHS fine post-smoother (16-step shifted GCR),
//           batched restriction (MultiRHSBlockProject / GEMM),
//           ONE coarse PGCR on the 6D mrhs coarse operator
//           (MultiGeneralCoarsenedMatrix, GEMM mults -- the ~10x win),
//           batched prolongation.
//
// The coarse operator is coarsened once with the standard single-RHS
// machinery (subspace cache reused) and imported via CopyMatrix.
//
// Memory note: outer restart history is 2*mmax*nrhs fine fields
// (49GB/field global at 48^3x96,Ls=24). Default mmax=8, nrhs=12 needs
// ~10TB for history; use OuterMmax / NRHS to fit the partition.
//
// Env vars: MASS, SUBSPACE_FILE, BLOCK (dotted e.g. 4.4.3.4),
//           NRHS (default 12, multiple of Nsimd),
//           FineSmootherShift, FineSmootherOrder,
//           CoarseSolverTol, CoarseSolverOrder,
//           OuterMmax, OuterNstep, OuterTol

#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidualNonHermitian.h>

using namespace std;
using namespace Grid;

RealD FineSmootherShift = 0.1;
int   FineSmootherOrder = 16;
RealD CoarseSolverTol   = 0.03;
int   CoarseSolverOrder = 200;
RealD OuterTol          = 1.0e-8;
int   OuterMmax         = 8;
int   OuterNstep        = 8;
int   Nrhs              = 12;
RealD mass              = 0.00078;

void ParseEnvironment(void)
{
  if(getenv("MASS"))              mass              = atof(getenv("MASS"));
  if(getenv("FineSmootherShift")) FineSmootherShift = atof(getenv("FineSmootherShift"));
  if(getenv("FineSmootherOrder")) FineSmootherOrder = atoi(getenv("FineSmootherOrder"));
  if(getenv("CoarseSolverTol"))   CoarseSolverTol   = atof(getenv("CoarseSolverTol"));
  if(getenv("CoarseSolverOrder")) CoarseSolverOrder = atoi(getenv("CoarseSolverOrder"));
  if(getenv("OuterTol"))          OuterTol          = atof(getenv("OuterTol"));
  if(getenv("OuterMmax"))         OuterMmax         = atoi(getenv("OuterMmax"));
  if(getenv("OuterNstep"))        OuterNstep        = atoi(getenv("OuterNstep"));
  if(getenv("NRHS"))              Nrhs              = atoi(getenv("NRHS"));

  std::cout << GridLogMessage << "PARAM: MASS              " << mass              << std::endl;
  std::cout << GridLogMessage << "PARAM: NRHS              " << Nrhs              << std::endl;
  std::cout << GridLogMessage << "PARAM: FineSmootherShift " << FineSmootherShift << std::endl;
  std::cout << GridLogMessage << "PARAM: FineSmootherOrder " << FineSmootherOrder << std::endl;
  std::cout << GridLogMessage << "PARAM: CoarseSolverTol   " << CoarseSolverTol   << std::endl;
  std::cout << GridLogMessage << "PARAM: CoarseSolverOrder " << CoarseSolverOrder << std::endl;
  std::cout << GridLogMessage << "PARAM: OuterTol          " << OuterTol          << std::endl;
  std::cout << GridLogMessage << "PARAM: OuterMmax         " << OuterMmax         << std::endl;
  std::cout << GridLogMessage << "PARAM: OuterNstep        " << OuterNstep        << std::endl;
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

//////////////////////////////////////////////////////////////////////
// Minimal multi-RHS function interface (preconditioner slot)
//////////////////////////////////////////////////////////////////////
template<class Field>
class MrhsLinearFunction {
public:
  virtual void operator()(std::vector<Field> &in, std::vector<Field> &out) = 0;
};

//////////////////////////////////////////////////////////////////////
// Single-polynomial multi-RHS PGCR (non-Hermitian).
//
// Verbatim adaptation of PrecGeneralisedConjugateResidualNonHermitian
// to std::vector<Field>: every innerProduct / norm2 is SUMMED over the
// RHS index, so one alpha/beta per step is shared by all RHS -- the
// single GCR on the enlarged block-diagonal system.
//////////////////////////////////////////////////////////////////////
template<class Field>
class MrhsPGCRNonHermitian {
public:
  RealD   Tolerance;
  Integer MaxIterations;
  int mmax;
  int nstep;
  int steps;
  int level;
  int ZeroGuess = 0; int FirstCycle = 0;  // caller contract: zero guess => first-cycle r0 = src
  std::string name = "Level 1";
  LinearOperatorBase<Field> &Linop;
  MrhsLinearFunction<Field> &Preconditioner;

  void Level(int lv) { name = "Level " + std::to_string(lv); level=lv; };
  void Name(std::string n) { name = n; };
  void SetZeroGuess(int z) { ZeroGuess=z; };

  MrhsPGCRNonHermitian(RealD tol,Integer maxit,
                       LinearOperatorBase<Field> &_Linop,
                       MrhsLinearFunction<Field> &Prec,
                       int _mmax,int _nstep)
    : Tolerance(tol), MaxIterations(maxit), Linop(_Linop), Preconditioner(Prec),
      mmax(_mmax), nstep(_nstep) { level=1; }

  ///////////////////////////////////////////////////////////////
  // vector-of-fields linear algebra, reductions summed over rhs
  ///////////////////////////////////////////////////////////////
  static RealD vnorm2(std::vector<Field> &x){
    RealD s=0.0; for(auto &f : x) s+=norm2(f); return s;
  }
  static ComplexD vinnerProduct(std::vector<Field> &x, std::vector<Field> &y){
    ComplexD s(0.0); for(int r=0;r<(int)x.size();r++) s+=innerProduct(x[r],y[r]); return s;
  }
  static void vaxpy(std::vector<Field> &z, ComplexD a, std::vector<Field> &x, std::vector<Field> &y){
    for(int r=0;r<(int)z.size();r++) axpy(z[r],a,x[r],y[r]);
  }
  void vOp(std::vector<Field> &in, std::vector<Field> &out){
    for(int r=0;r<(int)in.size();r++) Linop.Op(in[r],out[r]);
  }

  void operator() (std::vector<Field> &src, std::vector<Field> &psi){
    RealD cp, ssq, rsq;
    int nrhs = src.size();
    GridBase *grid = src[0].Grid();

    ssq=vnorm2(src);
    rsq=Tolerance*Tolerance*ssq;

    std::vector<Field> r(nrhs,grid);

    GridStopWatch SolverTimer;
    SolverTimer.Start();

    steps=0;
    FirstCycle=1;
    for(int k=0;k<MaxIterations;k++){

      cp=GCRnStep(src,psi,rsq);

      std::cout<<GridLogMessage<<std::string(level,'\t')<<" "<<name
               <<" MrhsPGCR("<<mmax<<","<<nstep<<") "<<steps<<" steps cp = "<<cp<<" target "<<rsq<<std::endl;

      if(cp<rsq){
        SolverTimer.Stop();
        vOp(psi,r);
        for(int rr=0;rr<nrhs;rr++) axpy(r[rr],-1.0,src[rr],r[rr]);
        RealD tr=vnorm2(r);
        std::cout<<GridLogMessage<<std::string(level,'\t')<<" "<<name
                 <<" MrhsPGCR: Converged on iteration "<<steps
                 <<" computed residual "<<std::sqrt(cp/ssq)
                 <<" true residual "<<std::sqrt(tr/ssq)
                 <<" target "<<Tolerance<<std::endl;
        std::cout<<GridLogMessage<<std::string(level,'\t')<<" "<<name
                 <<" MrhsPGCR Time elapsed: Total "<<SolverTimer.Elapsed()<<std::endl;
        // Per-RHS true residuals: the honest metric under the summed norm
        for(int rr=0;rr<nrhs;rr++){
          RealD rn = std::sqrt(norm2(r[rr])/norm2(src[rr]));
          std::cout<<GridLogMessage<<"MrhsPGCR per-rhs true residual["<<rr<<"] = "<<rn<<std::endl;
        }
        return;
      }
    }
    std::cout<<GridLogMessage<<"MrhsPGCR: did not converge"<<std::endl;
  }

  RealD GCRnStep(std::vector<Field> &src, std::vector<Field> &psi, RealD rsq){

    RealD cp;
    ComplexD a, b, rq;
    RealD zAAz;

    int nrhs = src.size();
    GridBase *grid = src[0].Grid();

    std::vector<Field> r (nrhs,grid);
    std::vector<Field> z (nrhs,grid);
    std::vector<Field> Az(nrhs,grid);

    ////////////////////////////////
    // history for flexible orthog: [mmax][nrhs]
    ////////////////////////////////
    std::vector< std::vector<Field> > q(mmax, std::vector<Field>(nrhs,grid));
    std::vector< std::vector<Field> > p(mmax, std::vector<Field>(nrhs,grid));
    std::vector<RealD> qq(mmax);

    std::cout<<GridLogMessage<<std::string(level,'\t')<<" "<<name<<" MrhsPGCR nStep("<<nstep<<")"<<std::endl;

    if (ZeroGuess && FirstCycle) {
      for(int rr=0;rr<nrhs;rr++){ psi[rr]=Zero(); r[rr]=src[rr]; }
    } else {
      vOp(psi,Az);
      for(int rr=0;rr<nrhs;rr++) r[rr] = src[rr]-Az[rr];
      std::cout<<GridLogMessage<<std::string(level,'\t')<<" "<<name
               <<" MrhsPGCR true residual r = src - A psi   "<<vnorm2(r)<<std::endl;
    }
    FirstCycle=0;

    Preconditioner(r,z);
    vOp(z,Az);
    zAAz=vnorm2(Az);

    p[0]=z;
    q[0]=Az;
    qq[0]=zAAz;

    cp=vnorm2(r);

    for(int k=0;k<nstep;k++){

      steps++;

      int kp     = k+1;
      int peri_k = k %mmax;
      int peri_kp= kp%mmax;

      rq = vinnerProduct(q[peri_k],r);
      a  = rq/qq[peri_k];

      vaxpy(psi,a,p[peri_k],psi);
      vaxpy(r,-a,q[peri_k],r);
      cp = vnorm2(r);

      std::cout<<GridLogMessage<<std::string(level,'\t')<<" "<<name
               <<" MrhsPGCR step["<<steps<<"]  resid "<<cp<<" target "<<rsq<<std::endl;

      if((k==nstep-1)||(cp<rsq)){
        return cp;
      }

      Preconditioner(r,z);
      vOp(z,Az);
      zAAz=vnorm2(Az);

      q[peri_kp]=Az;
      p[peri_kp]=z;

      int northog = ((kp)>(mmax-1))?(mmax-1):(kp);
      for(int back=0;back<northog;back++){
        int peri_back=(k-back)%mmax; GRID_ASSERT((k-back)>=0);
        b = -real(vinnerProduct(q[peri_back],Az))/qq[peri_back];
        vaxpy(p[peri_kp],b,p[peri_back],p[peri_kp]);
        vaxpy(q[peri_kp],b,q[peri_back],q[peri_kp]);
      }
      qq[peri_kp]=vnorm2(q[peri_kp]);
    }
    GRID_ASSERT(0); // never reached
    return cp;
  }
};

//////////////////////////////////////////////////////////////////////
// Trivial multi-RHS preconditioner
//////////////////////////////////////////////////////////////////////
template<class Field>
class TrivialMrhsPrecon : public MrhsLinearFunction<Field> {
public:
  virtual void operator()(std::vector<Field> &in, std::vector<Field> &out){
    for(int r=0;r<(int)in.size();r++) out[r]=in[r];
  }
};

//////////////////////////////////////////////////////////////////////
// MultiRHS two-level V-cycle.
//
// Mirrors MGPreconditioner in Example_pvdagm:
//   out = in (trivial pre)                      [per rhs]
//   r1  = in - A out                            [per rhs]
//   batched blockProject -> pack -> ONE mrhs coarse PGCR -> unpack
//     -> batched blockPromote; out += correction
//   r2  = in - A out                            [per rhs]
//   per-RHS fine post-smoother; out += smooth(r2)
//////////////////////////////////////////////////////////////////////
template<class FineField, class MrhsCoarseVector, class FineSmoother>
class MrhsTwoLevelMG : public MrhsLinearFunction<FineField> {
public:
  typedef MrhsCoarseVector CoarseVector;   // same lattice type on Coarse5d and CoarseMrhs

  LinearOperatorBase<FineField>            &_FineOperator;
  FineSmoother                             &_PostSmoother;   // single-RHS smoother, looped
  MultiRHSBlockProject<FineField>          &_Projector;
  LinearFunction<CoarseVector>             &_CoarseSolve;    // PGCR on the 6D mrhs field
  GridBase                                 *_CoarseGrid;     // Coarse5d (single rhs)
  GridBase                                 *_CoarseGridMrhs; // 6D

  MrhsTwoLevelMG(LinearOperatorBase<FineField> &FineOp,
                 FineSmoother                  &Post,
                 MultiRHSBlockProject<FineField> &Projector,
                 LinearFunction<CoarseVector>  &CoarseSolve,
                 GridBase *CoarseGrid, GridBase *CoarseGridMrhs)
    : _FineOperator(FineOp), _PostSmoother(Post), _Projector(Projector),
      _CoarseSolve(CoarseSolve), _CoarseGrid(CoarseGrid), _CoarseGridMrhs(CoarseGridMrhs) {}

  virtual void operator()(std::vector<FineField> &in, std::vector<FineField> &out){
    int nrhs = in.size();
    GridBase *fgrid = in[0].Grid();
    double t;

    std::vector<FineField> vec1(nrhs,fgrid);
    std::vector<FineField> vec2(nrhs,fgrid);

    // Trivial pre-smoother: out = in (as in Example_pvdagm with simple_fine)
    for(int r=0;r<nrhs;r++) out[r]=in[r];

    // Residual
    for(int r=0;r<nrhs;r++){
      _FineOperator.Op(out[r],vec1[r]);
      sub(vec1[r],in[r],vec1[r]);
    }

    // Batched fine->coarse, pack rhs into 6D field
    std::vector<CoarseVector> Csrc_split(nrhs,_CoarseGrid);
    std::vector<CoarseVector> Csol_split(nrhs,_CoarseGrid);
    CoarseVector CsrcMrhs(_CoarseGridMrhs);
    CoarseVector CsolMrhs(_CoarseGridMrhs);

    t=-usecond();
    _Projector.blockProject(vec1,Csrc_split);
    for(int r=0;r<nrhs;r++) InsertSliceFast(Csrc_split[r],CsrcMrhs,r,0);
    t+=usecond();
    std::cout<<GridLogMessage<<"Mrhs project+pack took "<<t/1000.0<<"ms"<<std::endl;

    // ONE coarse solve for all rhs (GEMM coarse mults)
    t=-usecond();
    CsolMrhs=Zero();
    _CoarseSolve(CsrcMrhs,CsolMrhs);
    t+=usecond();
    std::cout<<GridLogMessage<<"Mrhs coarse solve took "<<t/1000.0<<"ms  ("<<t/1000.0/nrhs<<"ms/rhs)"<<std::endl;

    // Unpack, batched coarse->fine, add correction
    t=-usecond();
    for(int r=0;r<nrhs;r++) ExtractSliceFast(Csol_split[r],CsolMrhs,r,0);
    _Projector.blockPromote(vec1,Csol_split);
    for(int r=0;r<nrhs;r++) add(out[r],out[r],vec1[r]);
    t+=usecond();
    std::cout<<GridLogMessage<<"Mrhs unpack+promote took "<<t/1000.0<<"ms"<<std::endl;

    // Residual
    for(int r=0;r<nrhs;r++){
      _FineOperator.Op(out[r],vec1[r]);
      sub(vec1[r],in[r],vec1[r]);
    }

    // Per-RHS post-smoother (fine level has no batching win; memory-light)
    t=-usecond();
    for(int r=0;r<nrhs;r++){
      vec2[r]=Zero();
      _PostSmoother(vec1[r],vec2[r]);
      add(out[r],out[r],vec2[r]);
    }
    t+=usecond();
    std::cout<<GridLogMessage<<"Mrhs post-smooth took "<<t/1000.0<<"ms  ("<<t/1000.0/nrhs<<"ms/rhs)"<<std::endl;
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
  const int nrhs   = Nrhs;

  GRID_ASSERT(nrhs % vComplex::Nsimd() == 0);

  std::cout << GridLogMessage << "MultiRHS PVdagM MG: mass=" << mass << " Ls=" << Ls
            << " nbasis=" << nbasis << " nrhs=" << nrhs << std::endl;

  std::vector<int> lat_size {48, 48, 48, 96};

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(lat_size, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  // Blocking, env-overridable: BLOCK=2.2.2.2
  Coordinate clatt = lat_size;
  Coordinate Block({4,4,3,4});
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

  ////////////////////////////////////////////////////////////
  // 6D multi-RHS coarse grid: rhs is dim 0, SIMD across rhs
  // (pattern: tests/debug/Test_general_coarse_hdcg_phys48.cc)
  ////////////////////////////////////////////////////////////
  Coordinate mpi=GridDefaultMpi();
  Coordinate rhMpi ({1,1,mpi[0],mpi[1],mpi[2],mpi[3]});
  Coordinate rhLatt({nrhs,1,clatt[0],clatt[1],clatt[2],clatt[3]});
  Coordinate rhSimd({vComplex::Nsimd(),1, 1,1,1,1});
  GridCartesian *CoarseMrhs = new GridCartesian(rhLatt,rhSimd,rhMpi);

  GridParallelRNG RNG5(FGrid);  RNG5.SeedFixedIntegers({5,6,7,8});
  GridParallelRNG RNG4(UGrid);  RNG4.SeedFixedIntegers({1,2,3,4});

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
  typedef MultiGeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> MrhsLittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector                           CoarseVector;
  typedef Aggregation<vSpinColourVector,vTComplex,nbasis>             Subspace;

  PVdagM_t        PVdagM(Ddwf,Dpv);
  ShiftedPVdagM_t ShiftedPVdagM(FineSmootherShift,Ddwf,Dpv);

  NextToNearestStencilGeometry5D geom(Coarse5d);

  ////////////////////////////////////////////////////////////
  // Subspace: load from cache or generate
  ////////////////////////////////////////////////////////////
  std::string subspace_file = "/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/subspace_nb"
                              + std::to_string(nbasis) + ".scidac";
  if ( getenv("SUBSPACE_FILE") ) subspace_file = std::string(getenv("SUBSPACE_FILE"));

  uint64_t file_exists = 0;
  if ( UGrid->IsBoss() ) { std::ifstream f(subspace_file); file_exists = f.good() ? 1 : 0; }
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

  ////////////////////////////////////////////////////////////
  // Coarsen once (single-RHS machinery), import into mrhs op.
  // NB CoarsenOperator block-orthogonalises subspace in place;
  // ImportBasis AFTER so the projector matches the coarse op.
  ////////////////////////////////////////////////////////////
  MrhsLittleDiracOperator mrhsLittleDiracOpPV(geom,CoarseMrhs);
  {
    // Scope the single-RHS operator so its padded _A is freed after import.
    // At small local volumes the depth-2 padded cell inflates ~8x
    // (e.g. 2^4 blocking on 432 ranks: local 8x4x4x12 -> padded 12x8x8x16,
    // ~23GB/GCD for _A alone -> OOM if kept alive).
    LittleDiracOperator LittleDiracOpPV(geom,FGrid,Coarse5d);
    LittleDiracOpPV.CoarsenOperator(PVdagM, AggregatesGCR);
    mrhsLittleDiracOpPV.CopyMatrix(LittleDiracOpPV);
  }

  MultiRHSBlockProject<LatticeFermionD> MrhsProjector;
  MrhsProjector.Allocate(nbasis,FGrid,Coarse5d);
  MrhsProjector.ImportBasis(AggregatesGCR.subspace);

  ////////////////////////////////////////////////////////////
  // Solvers
  ////////////////////////////////////////////////////////////
  NonHermitianLinearOperator<MrhsLittleDiracOperator,CoarseVector> mrhsLinOpCoarse(mrhsLittleDiracOpPV);

  TrivialPrecon<CoarseVector> simpleC;
  PrecGeneralisedConjugateResidualNonHermitian<CoarseVector>
    L2PGCRmrhs(CoarseSolverTol, CoarseSolverOrder/20, mrhsLinOpCoarse, simpleC, 20, 20);
  L2PGCRmrhs.Level(2);
  L2PGCRmrhs.Name("Couter");
  L2PGCRmrhs.SetZeroGuess(1);          // caller zeroes CsolMrhs

  TrivialPrecon<LatticeFermionD> simple_fine;
  PrecGeneralisedConjugateResidualNonHermitian<LatticeFermionD>
    SmootherGCR(0.0, 1, ShiftedPVdagM, simple_fine, FineSmootherOrder, FineSmootherOrder);
  SmootherGCR.Level(1);
  SmootherGCR.Name("Fsmoother");
  SmootherGCR.SetZeroGuess(1);         // caller zeroes vec2[r]

  typedef PrecGeneralisedConjugateResidualNonHermitian<LatticeFermionD> FineSmoother_t;
  MrhsTwoLevelMG<LatticeFermionD,CoarseVector,FineSmoother_t>
    TwoLevelPrecon(PVdagM, SmootherGCR, MrhsProjector, L2PGCRmrhs, Coarse5d, CoarseMrhs);

  MrhsPGCRNonHermitian<LatticeFermionD>
    L1PGCRmrhs(OuterTol, 1000, PVdagM, TwoLevelPrecon, OuterMmax, OuterNstep);
  L1PGCRmrhs.Level(1);
  L1PGCRmrhs.Name("Fouter");
  L1PGCRmrhs.SetZeroGuess(1);          // sol[r]=Zero() at source setup

  ////////////////////////////////////////////////////////////
  // Sources and solve
  ////////////////////////////////////////////////////////////
  std::vector<LatticeFermionD> src(nrhs,FGrid);
  std::vector<LatticeFermionD> sol(nrhs,FGrid);
  for(int r=0;r<nrhs;r++){
    gaussian(RNG5,src[r]);
    sol[r]=Zero();
  }

  std::cout << GridLogMessage << "**********************************************" << std::endl;
  std::cout << GridLogMessage << " MultiRHS two-level solve: " << nrhs << " RHS " << std::endl;
  std::cout << GridLogMessage << "**********************************************" << std::endl;

  GridStopWatch w; w.Start();
  L1PGCRmrhs(src,sol);
  w.Stop();

  std::cout << GridLogMessage << "MultiRHS solve total " << w.Elapsed()
            << "  (per RHS: " << w.useconds()/1.0e6/nrhs << " s)" << std::endl;

  // Independent final verification, per RHS
  {
    LatticeFermionD Ax(FGrid);
    RealD worst=0.0;
    for(int r=0;r<nrhs;r++){
      PVdagM.Op(sol[r],Ax);
      Ax = Ax - src[r];
      RealD rn = std::sqrt(norm2(Ax)/norm2(src[r]));
      std::cout << GridLogMessage << "FINAL: rhs["<<r<<"] true residual = " << rn << std::endl;
      worst = std::max(worst,rn);
    }
    std::cout << GridLogMessage << "FINAL: worst-case residual = " << worst << std::endl;
  }

  std::cout << GridLogMessage << "Done" << std::endl;
  Grid_finalize();
  return 0;
}
