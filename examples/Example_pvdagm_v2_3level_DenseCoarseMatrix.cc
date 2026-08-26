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
// STAGES ONE AND TWO: grids, types, subspace, and the L1 and L2 coarsenings.
// The dense bottom and the solves are not here yet.
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
// Env: LATT LS MASS NBASIS(compile time) NRHS BLOCK BLOCK2 COARSEN_BATCH
//      HOT_START CONFIG SUBSPACE_FILE V1_CHECK MRHS_COARSEN
//

#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidualNonHermitian.h>
#include <Grid/algorithms/multigrid/DenseCoarseMatrix.h>

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

// Solver tuning.  PRINCIPLE (PB, 2026-08-24): the defaults ARE the current
// optimum, so an unset environment reproduces the best banked result; they
// are updated as and when a better point is found, and every change is
// dated here.  Environment variables of the same names override for sweeps.
//
// Current optimum: 2026-08-24, slurm-5335492 F4, 48^3x96 Ls=24 on 288 GCDs,
//   17.2 s/RHS at Nrhs=4, 32.2 s at Nrhs=1 (exact-halo FINAL ~1e-8 pending
//   the exact-outer rerun).  Smoother mmax == order (full GCR history);
//   PB's mmax=1 trial gave 72 vs ~60 outer iterations and was slower.
RealD FineSmootherShift    = 0.1;
int   FineSmootherOrder    = 6;
int   FineSmootherMmax     = 6;
RealD CoarseSmootherShift  = 0.1;
int   PowerIterations      = 0;     // >0: power-iterate the smoother operators before the solves (spectral edge)
// Smoother implementation per level (Smoothers.h):
//   gcr    : the adaptive PGCR (default, as always)
//   replay : run the PGCR with a coefficient recorder for the first
//            PolyRecordIters outer steps, then switch to GCRReplaySmoother
//            (same polynomial, no inner products) -- 1402.2585 p.13 revisited
//   cheb   : ChebyshevNonHermitianSmoother, 1/x on [ChebLo,ChebHi], order
//            = the GCR step count of that level
std::string FineSmootherMode   = "gcr";
std::string CoarseSmootherMode = "gcr";
int   PolyRecordIters      = 4;
RealD FineChebLo   = 3.0,  FineChebHi   = 137.0;   // from the harvested polynomial and PowerIterations edge
RealD CoarseChebLo = 8.0,  CoarseChebHi = 45.0;
int   CoarseSmootherNstep  = 2;
int   CoarseSmootherMmax   = 2;
RealD CoarseSolverTol      = 0.05;
int   CoarseSolverOrder    = 200;
int   CoarseSolverMmax     = 16;
RealD OuterTol             = 1.0e-8;
int   OuterMmax            = 4;
int   OuterNstep           = 8;

// "It's legal to get the same answer faster, not to get a less correct
// answer."  (PB, 2026-08-24)
//
// Halo-precision POLICY: reduced-precision (fp32 wire)
// halos belong in the PRECONDITIONER -- the smoother, the V-cycle's own
// residuals, and the coarsening -- and NEVER in the outer Krylov.  The
// outer operator's applications define what "converged" means; making them
// sloppy turns the stopping criterion into a statement about the wrong
// operator (measured: solver stops at computed 9.8e-9 while the true
// residual is 3.4e-8).  Exactness costs one exact fine matvec per outer
// iteration, a few percent of the solve.
//
// Stencil::SloppyComms is a free runtime setter (Stencil.h:303), so the
// policy is implemented by SCOPED toggling: SetFineSloppy(1) on entering
// the preconditioner / coarsening, SetFineSloppy(0) on leaving.  The
// operators default to EXACT.  FineSloppyComms therefore now means
// "sloppy inside the preconditioner"; =0 makes everything exact.
int   FineSloppyComms      = 1;
// SmootherCoeffLog=1 : print the GCR step lengths a_k and orthogonalisation
// coefficients b_kj of BOTH smoothers every call -- the harvest for a fixed
// polynomial smoother (stable coefficients => stationary p(A), no reductions).
int   SmootherCoeffLog     = 0;
std::function<void(int)> SetFineSloppy = [](int){};

void ParseEnvironment(void)
{
  if(getenv("MASS"))          mass        = atof(getenv("MASS"));
  if(getenv("NRHS"))          Nrhs        = atoi(getenv("NRHS"));
  if(getenv("LS"))            Ls          = atoi(getenv("LS"));
  if(getenv("COARSEN_BATCH")) CoarsenBatch= atoi(getenv("COARSEN_BATCH"));
  if(getenv("FineSmootherShift"))  FineSmootherShift  = atof(getenv("FineSmootherShift"));
  if(getenv("FineSmootherOrder"))  FineSmootherOrder  = atoi(getenv("FineSmootherOrder"));
  if(getenv("FineSmootherMmax"))   FineSmootherMmax   = atoi(getenv("FineSmootherMmax"));
  if(getenv("CoarseSmootherShift"))CoarseSmootherShift= atof(getenv("CoarseSmootherShift"));
  if(getenv("PowerIterations"))    PowerIterations    = atoi(getenv("PowerIterations"));
  if(getenv("FineSmootherMode"))   FineSmootherMode   = getenv("FineSmootherMode");
  if(getenv("CoarseSmootherMode")) CoarseSmootherMode = getenv("CoarseSmootherMode");
  if(getenv("PolyRecordIters"))    PolyRecordIters    = atoi(getenv("PolyRecordIters"));
  if(getenv("FineChebLo"))         FineChebLo         = atof(getenv("FineChebLo"));
  if(getenv("FineChebHi"))         FineChebHi         = atof(getenv("FineChebHi"));
  if(getenv("CoarseChebLo"))       CoarseChebLo       = atof(getenv("CoarseChebLo"));
  if(getenv("CoarseChebHi"))       CoarseChebHi       = atof(getenv("CoarseChebHi"));
  if(getenv("CoarseSmootherNstep"))CoarseSmootherNstep= atoi(getenv("CoarseSmootherNstep"));
  if(getenv("CoarseSmootherMmax")) CoarseSmootherMmax = atoi(getenv("CoarseSmootherMmax"));
  if(getenv("CoarseSolverTol"))    CoarseSolverTol    = atof(getenv("CoarseSolverTol"));
  if(getenv("CoarseSolverOrder"))  CoarseSolverOrder  = atoi(getenv("CoarseSolverOrder"));
  if(getenv("CoarseSolverMmax"))   CoarseSolverMmax   = atoi(getenv("CoarseSolverMmax"));
  if(getenv("OuterTol"))           OuterTol           = atof(getenv("OuterTol"));
  if(getenv("OuterMmax"))          OuterMmax          = atoi(getenv("OuterMmax"));
  if(getenv("FineSloppyComms"))    FineSloppyComms    = atoi(getenv("FineSloppyComms"));
  if(getenv("SmootherCoeffLog"))   SmootherCoeffLog   = atoi(getenv("SmootherCoeffLog"));
  if(getenv("OuterNstep"))         OuterNstep         = atoi(getenv("OuterNstep"));
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
// ||<v|v> - I||_F over a set of coarse vectors. Small means the raw near null
// content survived the projection; see GramGuard for where a leak lands.
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

// On a leak every image collapses to the block unit e_k, the Gram becomes
// N*I, and the defect lands at (N-1)*sqrt(nbasis) -- orders above the ~0.2
// of a content preserving projection. Trip well below that so a mis-set
// threshold costs a log line rather than the run.
template<class CoarseField>
void GramGuard(const std::string &name,std::vector<CoarseField> &v,GridBase *grid)
{
  RealD defect = GramDefect(v);
  RealD N      = (RealD)grid->gSites();
  RealD leak   = (N-1.0)*std::sqrt((RealD)v.size());
  RealD trip   = std::sqrt(N);
  std::cout << GridLogMessage << "GUARD: ||<"<<name<<"|"<<name<<"> - I||_F = " << defect
	    << "   (e_k leak would be " << leak << ", trip at " << trip << ")" << std::endl;
  GRID_ASSERT( defect < trip );
}


//////////////////////////////////////////////////////////////////////
// Shifted variants for the smoothers
//////////////////////////////////////////////////////////////////////
template<class Matrix,class Field>
class ShiftedPVdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat; Matrix &_PV;
public:
  RealD shift;
  ShiftedPVdagMLinearOperator(RealD _shift,Matrix &Mat,Matrix &PV): shift(_shift),_Mat(Mat),_PV(PV){};
  void OpDiag (const Field &in, Field &out) { assert(0); }
  void OpDir  (const Field &in, Field &out,int dir,int disp) { assert(0); }
  void OpDirAll  (const Field &in, std::vector<Field> &out){ assert(0); };
  void Op     (const Field &in, Field &out){ Field tmp(in.Grid()); _Mat.M(in,tmp); _PV.Mdag(tmp,out); out = out + shift*in; }
  void AdjOp  (const Field &in, Field &out){ Field tmp(in.Grid()); _PV.M(tmp,out); _Mat.Mdag(in,tmp); out = out + shift*in; }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){ assert(0); }
  void HermOp(const Field &in, Field &out){ Field tmp(in.Grid()); Op(in,tmp); AdjOp(tmp,out); }
};

template<class Field>
class ShiftedLinearOperator : public LinearOperatorBase<Field> {
  LinearOperatorBase<Field> &_Op; RealD shift;
public:
  ShiftedLinearOperator(RealD _shift, LinearOperatorBase<Field> &Op) : _Op(Op), shift(_shift) {}
  void OpDiag  (const Field &in, Field &out) { assert(0); }
  void OpDir   (const Field &in, Field &out,int dir,int disp) { assert(0); }
  void OpDirAll (const Field &in, std::vector<Field> &out) { assert(0); }
  void Op      (const Field &in, Field &out) { _Op.Op(in,out);    out = out + shift*in; }
  void AdjOp   (const Field &in, Field &out) { _Op.AdjOp(in,out); out = out + shift*in; }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){ assert(0); }
  void HermOp  (const Field &in, Field &out) { Field tmp(in.Grid()); Op(in,tmp); AdjOp(tmp,out); }
};

//////////////////////////////////////////////////////////////////////
// Power iteration on a (non-Hermitian) operator: the spectral edge the
// smoother polynomial must not exceed.  Reports
//   step 0 : |A v|/|v| on a RANDOM unit v -- a one-sample lower bound on
//            sigma_max(A).  If this and the converged value agree, the
//            operator is near-normal and the spectral picture (R_m(lambda)
//            on the spectrum) is trustworthy; if not, the field of values
//            sets the safe interval and the spectrum understates it.
//   step k : |A v_k|/|v_k| -> |lambda_max| as v_k -> the dominant
//            eigenvector; the complex Rayleigh quotient <v,Av> gives its
//            phase (real => on the axis).  A non-converging oscillation
//            means a complex-conjugate pair of equal modulus at the top.
// Uses Op(), not HermOp(): this is the operator the smoother sees.
//////////////////////////////////////////////////////////////////////
template<class Field>
void PowerIteration(const std::string &name, LinearOperatorBase<Field> &Op, GridBase *grid, int iters)
{
  GRID_TRACE("PowerIteration");
  GridParallelRNG RNG(grid); RNG.SeedFixedIntegers(std::vector<int>({7,11,13,17}));
  Field v(grid), Av(grid);
  gaussian(RNG,v);
  RealD nv = std::sqrt(norm2(v)); v = v*(1.0/nv);
  RealD ratio=0.0, ratio0=0.0; ComplexD rq(0.0);
  for(int i=0;i<iters;i++){
    Op.Op(v,Av);
    RealD nAv = std::sqrt(norm2(Av));
    rq = innerProduct(v,Av);                 // v is unit
    ratio = nAv;
    if ( i==0 ) {
      ratio0 = ratio;
      std::cout << GridLogMessage << "PowerIteration " << name << " step 0 (random v): |Av|/|v| = " << ratio
                << "   [lower bound on sigma_max]" << std::endl;
    }
    if ( (i%10==0) || (i==iters-1) )
      std::cout << GridLogMessage << "PowerIteration " << name << " step " << i << " |Av|/|v| = " << ratio
                << "  Rayleigh = (" << real(rq) << "," << imag(rq) << ")" << std::endl;
    v = Av*(1.0/nAv);
  }
  std::cout << GridLogMessage << "PowerIteration " << name << " SUMMARY: |lambda_max| ~ " << ratio
            << "  Rayleigh (" << real(rq) << "," << imag(rq) << ")"
            << "  phase " << std::atan2(imag(rq),real(rq)) << " rad"
            << "  step-0 ratio / converged = " << ratio0/ratio
            << (ratio0/ratio > 1.2 ? "   ** non-normal: sigma_max well above |lambda_max| **" : "   (near-normal)")
            << std::endl;
}

//////////////////////////////////////////////////////////////////////
// Dense L3 solve on the packed D+1 coarse-coarse field
//////////////////////////////////////////////////////////////////////
template<class DenseType, class CoarseCoarseField>
class MrhsDenseCCSolve : public LinearFunction<CoarseCoarseField> {
public:
  DenseType &_Dense;
  int _nrhs;
  MrhsDenseCCSolve(DenseType &D, int nrhs) : _Dense(D), _nrhs(nrhs) {}
  using LinearFunction<CoarseCoarseField>::operator();
  virtual void operator()(const CoarseCoarseField &in, CoarseCoarseField &out){
    _Dense.ApplyBatch6D(in, out, _nrhs);
  }
};

//////////////////////////////////////////////////////////////////////
// mrhs interfaces + single-polynomial mrhs PGCR
//////////////////////////////////////////////////////////////////////
template<class Field>
class MrhsLinearFunction {
public:
  virtual void operator()(std::vector<Field> &in, std::vector<Field> &out) = 0;
};

template<class Field>
class MrhsPGCRNonHermitian {
public:
  RealD Tolerance; Integer MaxIterations; int mmax,nstep,steps,level;
  int ZeroGuess = 0; int FirstCycle = 0;
  std::string name = "Level 1";
  LinearOperatorBase<Field> &Linop;
  MrhsLinearFunction<Field> &Preconditioner;
  std::function<void(int)> OnStep;      // called with the outer step count after every step (smoother switching)
  void Level(int lv){ name = "Level " + std::to_string(lv); level=lv; }
  void Name(std::string n){ name = n; }
  void SetZeroGuess(int z){ ZeroGuess=z; }
  MrhsPGCRNonHermitian(RealD tol,Integer maxit,LinearOperatorBase<Field> &_Linop,MrhsLinearFunction<Field> &Prec,int _mmax,int _nstep)
    : Tolerance(tol),MaxIterations(maxit),Linop(_Linop),Preconditioner(Prec),mmax(_mmax),nstep(_nstep){ level=1; }
  static RealD vnorm2(std::vector<Field> &x){ RealD s=0; for(auto &f:x) s+=norm2(f); return s; }
  static ComplexD vinnerProduct(std::vector<Field> &x,std::vector<Field> &y){ ComplexD s(0); for(int r=0;r<(int)x.size();r++) s+=innerProduct(x[r],y[r]); return s; }
  static void vaxpy(std::vector<Field> &z,ComplexD a,std::vector<Field> &x,std::vector<Field> &y){ for(int r=0;r<(int)z.size();r++) axpy(z[r],a,x[r],y[r]); }
  void vOp(std::vector<Field> &in,std::vector<Field> &out){ GRID_TRACE("MrhsPGCR::vOp"); for(int r=0;r<(int)in.size();r++) Linop.Op(in[r],out[r]); }
  void operator()(std::vector<Field> &src,std::vector<Field> &psi){
    RealD cp,ssq,rsq; int nrhs=src.size(); GridBase *grid=src[0].Grid();
    ssq=vnorm2(src); rsq=Tolerance*Tolerance*ssq;
    std::vector<Field> r(nrhs,grid);
    GridStopWatch T; T.Start(); steps=0; FirstCycle=1;
    for(int k=0;k<MaxIterations;k++){
      cp=GCRnStep(src,psi,rsq);
      std::cout<<GridLogMessage<<std::string(level,'\t')<<" "<<name<<" MrhsPGCR("<<mmax<<","<<nstep<<") "<<steps<<" steps cp = "<<cp<<" target "<<rsq<<std::endl;
      if(cp<rsq){
        T.Stop(); vOp(psi,r); for(int rr=0;rr<nrhs;rr++) axpy(r[rr],-1.0,src[rr],r[rr]);
        RealD tr=vnorm2(r);
        std::cout<<GridLogMessage<<std::string(level,'\t')<<" "<<name<<" MrhsPGCR: Converged on iteration "<<steps
                 <<" computed residual "<<std::sqrt(cp/ssq)<<" true residual "<<std::sqrt(tr/ssq)<<" target "<<Tolerance<<std::endl;
        std::cout<<GridLogMessage<<std::string(level,'\t')<<" "<<name<<" MrhsPGCR Time elapsed: Total "<<T.Elapsed()<<std::endl;
        return;
      }
    }
    std::cout<<GridLogMessage<<"MrhsPGCR: did not converge"<<std::endl;
  }
  RealD GCRnStep(std::vector<Field> &src,std::vector<Field> &psi,RealD rsq){
    RealD cp; ComplexD a,b,rq; int nrhs=src.size(); GridBase *grid=src[0].Grid();
    std::vector<Field> r(nrhs,grid),Az(nrhs,grid);   // Az: restart residual scratch only
    std::vector< std::vector<Field> > q(mmax,std::vector<Field>(nrhs,grid));
    std::vector< std::vector<Field> > p(mmax,std::vector<Field>(nrhs,grid));
    std::vector<RealD> qq(mmax);
    if (ZeroGuess && FirstCycle) { for(int rr=0;rr<nrhs;rr++){ psi[rr]=Zero(); r[rr]=src[rr]; } }
    else                         { vOp(psi,Az); for(int rr=0;rr<nrhs;rr++) r[rr]=src[rr]-Az[rr]; }
    FirstCycle=0;
    // p[0]=Prec(r), q[0]=A p[0], produced directly in the history slots (no copies)
    Preconditioner(r,p[0]); vOp(p[0],q[0]); qq[0]=vnorm2(q[0]); cp=vnorm2(r);
    for(int k=0;k<nstep;k++){
      steps++; int kp=k+1, peri_k=k%mmax, peri_kp=kp%mmax;
      if ( OnStep ) OnStep(steps);
      rq=vinnerProduct(q[peri_k],r); a=rq/qq[peri_k];
      vaxpy(psi,a,p[peri_k],psi); vaxpy(r,-a,q[peri_k],r); cp=vnorm2(r);
      std::cout<<GridLogMessage<<std::string(level,'\t')<<" "<<name<<" MrhsPGCR step["<<steps<<"]  resid "<<cp<<" target "<<rsq<<std::endl;
      if((k==nstep-1)||(cp<rsq)) return cp;
      // New direction straight into its history slot: p=Prec(r), q=A p.
      Preconditioner(r,p[peri_kp]);
      vOp(p[peri_kp],q[peri_kp]);
      int northog=((kp)>(mmax-1))?(mmax-1):(kp);
      {
	GRID_TRACE("MrhsPGCR orthog");
        // Classical Gram-Schmidt: all coefficients against the UN-updated new q
        // (independent, batchable), then apply.  Complex coefficient: the
        // operator is non-Hermitian, real(<q_j,Aq>) alone left q's non-orthogonal.
        // Batched per rhs (one fused kernel + one reduction each), the shared
        // coefficient summed over rhs on the host, ONE GlobalSumVector.
        std::vector<ComplexD> bcoef(northog,ComplexD(0.0)), part;
        for(int rr=0;rr<nrhs;rr++){
          std::vector<const Field*> qwin(northog);
          for(int back=0;back<northog;back++){ int peri_back=(k-back)%mmax; GRID_ASSERT((k-back)>=0); qwin[back]=&q[peri_back][rr]; }
          rankInnerProductMulti(part,qwin,q[peri_kp][rr]);
          for(int back=0;back<northog;back++) bcoef[back]+=part[back];
        }
        if(northog) grid->GlobalSumVector(&bcoef[0],northog);
        for(int back=0;back<northog;back++){ int peri_back=(k-back)%mmax; bcoef[back]=-bcoef[back]/qq[peri_back]; }
        for(int rr=0;rr<nrhs;rr++){
          std::vector<const Field*> qwin(northog), pwin(northog);
          for(int back=0;back<northog;back++){ int peri_back=(k-back)%mmax; qwin[back]=&q[peri_back][rr]; pwin[back]=&p[peri_back][rr]; }
          axpyMulti(p[peri_kp][rr],bcoef,pwin);
          axpyMulti(q[peri_kp][rr],bcoef,qwin);
        }
      }
      qq[peri_kp]=vnorm2(q[peri_kp]);
    }
    GRID_ASSERT(0); return cp;
  }
};

//////////////////////////////////////////////////////////////////////
// L2->L3 mrhs V-cycle on the D+1 coarse field
//////////////////////////////////////////////////////////////////////
template<class CoarseField, class CoarseCoarseField>
class MrhsCoarseThreeLevelPrec : public LinearFunction<CoarseField> {
public:
  LinearOperatorBase<CoarseField>          &_CoarseOp;
  LinearFunction<CoarseField>              &_CoarseSmoother;
  MultiRHSBlockProject<CoarseField>        &_Projector;
  LinearFunction<CoarseCoarseField>        &_CoarseCoarseSolve;
  GridBase *_Coarse5d, *_CoarseCoarse5d, *_CoarseCoarseMrhs;
  int _nrhs;
  MrhsCoarseThreeLevelPrec(LinearOperatorBase<CoarseField> &CoarseOp,
                           LinearFunction<CoarseField> &CoarseSmoother,
                           MultiRHSBlockProject<CoarseField> &Projector,
                           LinearFunction<CoarseCoarseField> &CoarseCoarseSolve,
                           GridBase *Coarse5d, GridBase *CoarseCoarse5d, GridBase *CoarseCoarseMrhs, int nrhs)
    : _CoarseOp(CoarseOp), _CoarseSmoother(CoarseSmoother), _Projector(Projector),
      _CoarseCoarseSolve(CoarseCoarseSolve),
      _Coarse5d(Coarse5d), _CoarseCoarse5d(CoarseCoarse5d), _CoarseCoarseMrhs(CoarseCoarseMrhs), _nrhs(nrhs) {}
  using LinearFunction<CoarseField>::operator();
  virtual void operator()(const CoarseField &in, CoarseField &out) {
    int nrhs=_nrhs;
    CoarseField vec1(in.Grid());
    CoarseField vec2(in.Grid());
    out = in;
    _CoarseOp.Op(out,vec1);  sub(vec1,in,vec1);

    // restrict, through the mixed blockProject: D+1 coarse in, D+1 cc out
    CoarseCoarseField CCsrc(_CoarseCoarseMrhs);
    CoarseCoarseField CCsol(_CoarseCoarseMrhs);
    _Projector.blockProject(vec1,CCsrc);

    //    CCsol=Zero(); // Is this necessary?
    _CoarseCoarseSolve(CCsrc,CCsol);

    _Projector.blockPromote(vec1,CCsol);
    add(out,out,vec1);

    _CoarseOp.Op(out,vec1);  sub(vec1,in,vec1);
    //    vec2=Zero(); // Zero guess
    _CoarseSmoother(vec1,vec2);
    add(out,out,vec2);
  }
};

//////////////////////////////////////////////////////////////////////
// L1->L2 mrhs V-cycle
//////////////////////////////////////////////////////////////////////
template<class FineField, class MrhsCoarseVector, class FineSmoother>
class MrhsTwoLevelMG : public MrhsLinearFunction<FineField> {
public:
  typedef MrhsCoarseVector CoarseVector;
  LinearOperatorBase<FineField>   &_FineOperator;
  FineSmoother                    &_PostSmoother;
  MultiRHSBlockProject<FineField> &_Projector;
  LinearFunction<CoarseVector>    &_CoarseSolve;
  GridBase *_CoarseGrid, *_CoarseGridMrhs;
  MrhsTwoLevelMG(LinearOperatorBase<FineField> &FineOp, FineSmoother &Post,
                 MultiRHSBlockProject<FineField> &Projector, LinearFunction<CoarseVector> &CoarseSolve,
                 GridBase *CoarseGrid, GridBase *CoarseGridMrhs)
    : _FineOperator(FineOp),_PostSmoother(Post),_Projector(Projector),_CoarseSolve(CoarseSolve),
      _CoarseGrid(CoarseGrid),_CoarseGridMrhs(CoarseGridMrhs){}
  virtual void operator()(std::vector<FineField> &in, std::vector<FineField> &out){
    // The whole V-cycle is preconditioner: its fine residuals and the
    // smoother run with sloppy halos; the caller (the outer Krylov) gets
    // the exact operator back on exit.
    GRID_TRACE("MGVcycle");
    SetFineSloppy(FineSloppyComms);
    int nrhs=in.size(); GridBase *fgrid=in[0].Grid();
    std::vector<FineField> vec1(nrhs,fgrid),vec2(nrhs,fgrid);
    for(int r=0;r<nrhs;r++) out[r]=in[r];
    { GRID_TRACE("MGFineResidual");
      for(int r=0;r<nrhs;r++){ _FineOperator.Op(out[r],vec1[r]); sub(vec1[r],in[r],vec1[r]); }
    }
    // fine vector -> D+1 coarse, via the mixed blockProject
    CoarseVector CsrcMrhs(_CoarseGridMrhs), CsolMrhs(_CoarseGridMrhs);
    { GRID_TRACE("MGProject");
      _Projector.blockProject(vec1,CsrcMrhs);
    }
    CsolMrhs=Zero();
    { GRID_TRACE("MGCoarseSolve");
      _CoarseSolve(CsrcMrhs,CsolMrhs);
    }
    { GRID_TRACE("MGPromote");
      _Projector.blockPromote(vec1,CsolMrhs);
      for(int r=0;r<nrhs;r++) add(out[r],out[r],vec1[r]);
    }
    { GRID_TRACE("MGFineResidual2");
      for(int r=0;r<nrhs;r++){ _FineOperator.Op(out[r],vec1[r]); sub(vec1[r],in[r],vec1[r]); }
    }
    { GRID_TRACE("MGPostSmooth");
      for(int r=0;r<nrhs;r++){
	//	vec2[r]=Zero();
	_PostSmoother(vec1[r],vec2[r]); add(out[r],out[r],vec2[r]);
      }
    }
    SetFineSloppy(0);
  }
};

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
  Coordinate Block({2,2,3,3});         // L1->L2 blocking; banked optimum 2026-08-24 (env BLOCK overrides)
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

  // PVdagM and ShiftedPVdagM are thin wrappers over these same two objects,
  // so this one callback controls every fine halo in the program.  Default
  // EXACT; the preconditioner and the coarsening turn sloppiness on for
  // their own scope only (policy note at FineSloppyComms).
  SetFineSloppy = [&Ddwf,&Dpv](int sloppy){
    Ddwf.SloppyComms(sloppy);
    Dpv .SloppyComms(sloppy);
  };
  SetFineSloppy(0);
  std::cout << GridLogMessage << "Fine halo policy: preconditioner+coarsening "
            << (FineSloppyComms ? "SLOPPY (fp32 wire)" : "exact")
            << ", outer Krylov EXACT" << std::endl;

  typedef PVdagMLinearOperator<MobiusFermionD,LatticeFermionD>        PVdagM_t;
  typedef ShiftedPVdagMLinearOperator<MobiusFermionD,LatticeFermionD> ShiftedPVdagM_t;
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
  SetFineSloppy(FineSloppyComms);   // coarsening builds the PRECONDITIONER
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
  SetFineSloppy(0);

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

  GramGuard("psi_coarse",psi_coarse,Coarse5d);

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
    SetFineSloppy(FineSloppyComms);
    LittleDiracOpPV.CoarsenOperator(PVdagM,AggV);
    SetFineSloppy(0);

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
  Coordinate Block2({4,4,2,4});        // L2->L3 blocking; banked optimum 2026-08-24 (env BLOCK2 overrides)
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

    GramGuard("psi_cc",psi_cc,CoarseCoarse5d);
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

  //////////////////////////////////////////////////////////////////////
  // STAGE THREE (part one): the dense bottom on L2.
  //
  // DenseCoarseMatrix is bilingual: it takes the elements through
  // Geometry()/ExtractMatrix(), so the V2 operator serves directly. It does
  // detect that a multiRHS op cannot apply on the D dimensional grid and
  // skips its own certificate and VERIFY, so the equivalent check is done
  // here instead, driving the L2 operator at Nrhs 1 through a slice.
  //////////////////////////////////////////////////////////////////////
  typedef DenseCoarseMatrix<CComplexS2,nbasis> DenseCC_t;
  std::unique_ptr<DenseCC_t> DenseCC;

  if ( getenv("DENSE_CC")==nullptr || atoi(getenv("DENSE_CC")) ) {

    std::cout << GridLogMessage << "*** L3 dense bottom: import from the V2 L2 operator ***" << std::endl;
    DenseCC.reset(new DenseCC_t(CoarseCoarse5d));
    DenseCC->Import(CoarseOpL2);

    ////////////////////////////////////////////////////////////////////
    // ||A Ainv x - x|| / ||x||, the check Import could not run itself
    ////////////////////////////////////////////////////////////////////
    Coordinate cc1latt({1,1,cclatt[0],cclatt[1],cclatt[2],cclatt[3]});
    GridCartesian *CoarseCoarseOne = new GridCartesian(cc1latt,cmsimd,cmmpi);

    CoarseOpL2.SetGrid(CoarseCoarseOne);

    CoarseCoarseVector x(CoarseCoarse5d),y(CoarseCoarse5d),z(CoarseCoarse5d);
    GridParallelRNG dRNG(CoarseCoarse5d); dRNG.SeedFixedIntegers({11,12,13,14});
    random(dRNG,x);

    (*DenseCC)(x,y);                    // y = Ainv x

    CoarseCoarseVector y1(CoarseCoarseOne),z1(CoarseCoarseOne);
    InsertSliceFast(y,y1,0,0);
    CoarseOpL2.M(y1,z1);                // z = A y
    ExtractSliceFast(z,z1,0,0);

    z = z - x;
    RealD rel = std::sqrt(norm2(z)/norm2(x));
    std::cout << GridLogMessage << "L3 dense: ||A Ainv x - x||/||x|| = " << rel << std::endl;
    GRID_ASSERT( rel < 1.0e-2 );

    CoarseOpL2.SetGrid(CoarseCoarseMrhs);
    delete CoarseCoarseOne;
  }

  //////////////////////////////////////////////////////////////////////
  // STAGE THREE (part two): the solves.
  //
  // Both operators are driven from the SAME objects at whatever Nrhs is
  // asked for -- the matrix elements were built once and survive SetGrid --
  // so single RHS and multiRHS are the same code path with a different grid.
  //////////////////////////////////////////////////////////////////////
  GRID_ASSERT(DenseCC != nullptr);   // the PGCR bottom is not ported yet

  typedef PrecGeneralisedConjugateResidualNonHermitian<LatticeFermionD> FineSmoother_t;

  ShiftedPVdagM_t ShiftedPVdagM(FineSmootherShift,Ddwf,Dpv);
  TrivialPrecon<LatticeFermionD> simple_fine;
  TrivialPrecon<CoarseVector>    simpleC;

  auto RunSolve = [&](int nr)
  {
    std::cout << GridLogMessage << "**********************************************" << std::endl;
    std::cout << GridLogMessage << " V2 THREE-level solve, Nrhs = " << nr << std::endl;
    std::cout << GridLogMessage << "**********************************************" << std::endl;

    Coordinate cml({nr,1,clatt[0],clatt[1],clatt[2],clatt[3]});
    Coordinate ccml({nr,1,cclatt[0],cclatt[1],cclatt[2],cclatt[3]});
    GridCartesian *CMrhs  = new GridCartesian(cml, cmsimd,cmmpi);
    GridCartesian *CCMrhs = new GridCartesian(ccml,cmsimd,cmmpi);

    CoarseOpPV.SetGrid(CMrhs);
    CoarseOpL2.SetGrid(CCMrhs);

    NonHermitianLinearOperator<CoarseOperator,CoarseVector>            LinOpC (CoarseOpPV);
    NonHermitianLinearOperator<CoarseCoarseOperator,CoarseCoarseVector> LinOpCC(CoarseOpL2);

    MrhsDenseCCSolve<DenseCC_t,CoarseCoarseVector> ccSolve(*DenseCC,nr);

    ShiftedLinearOperator<CoarseVector> ShiftedC(CoarseSmootherShift, LinOpC);
    if ( PowerIterations > 0 ) {
      // Spectral edges the smoother polynomials must respect (see
      // scripts/gcr_polynomial.py: |R_m|>1 beyond the edge = amplification).
      PowerIteration<CoarseVector>   ("CoarseSmootherOp(shift="+std::to_string(CoarseSmootherShift)+")", ShiftedC, CMrhs, PowerIterations);
      PowerIteration<CoarseVector>   ("CoarseOp(unshifted)", LinOpC, CMrhs, PowerIterations);
      PowerIteration<LatticeFermionD>("FineSmootherOp(shift="+std::to_string(FineSmootherShift)+")", ShiftedPVdagM, Ddwf.FermionGrid(), PowerIterations);
    }
    PrecGeneralisedConjugateResidualNonHermitian<CoarseVector>
      CoarseSmootherGCR(0.01,1,ShiftedC,simpleC,CoarseSmootherMmax,CoarseSmootherNstep);
    CoarseSmootherGCR.Level(2); CoarseSmootherGCR.Name("Csmoother"); CoarseSmootherGCR.SetZeroGuess(1);

    SwitchableSmoother<CoarseVector> CoarseSmootherSlot(CoarseSmootherGCR,"Csmoother GCR");
    MrhsCoarseThreeLevelPrec<CoarseVector,CoarseCoarseVector>
      L2to3Precon(LinOpC, CoarseSmootherSlot, MrhsProjectorL2, ccSolve,
                  Coarse5d, CoarseCoarse5d, CCMrhs, nr);

    PrecGeneralisedConjugateResidualNonHermitian<CoarseVector>
      L2PGCR(CoarseSolverTol, CoarseSolverOrder/16, LinOpC, L2to3Precon, CoarseSolverMmax, 16);
    L2PGCR.Level(2); L2PGCR.Name("Couter"); L2PGCR.SetZeroGuess(1);

    FineSmoother_t SmootherGCR(0.0,1,ShiftedPVdagM,simple_fine,FineSmootherMmax,FineSmootherOrder);
    SmootherGCR.Level(1); SmootherGCR.Name("Fsmoother"); SmootherGCR.SetZeroGuess(1);
    SmootherGCR.LogCoefficients(SmootherCoeffLog);
    CoarseSmootherGCR.LogCoefficients(SmootherCoeffLog);

    SwitchableSmoother<LatticeFermionD> FineSmootherSlot(SmootherGCR,"Fsmoother GCR");
    MrhsTwoLevelMG<LatticeFermionD,CoarseVector,SwitchableSmoother<LatticeFermionD> >
      ThreeLevelPrecon(PVdagM, FineSmootherSlot, MrhsProjector, L2PGCR, Coarse5d, CMrhs);

    MrhsPGCRNonHermitian<LatticeFermionD>
      L1PGCR(OuterTol,1000,PVdagM,ThreeLevelPrecon,OuterMmax,OuterNstep);
    L1PGCR.Level(1); L1PGCR.Name("Fouter"); L1PGCR.SetZeroGuess(1);

    //////////////////////////////////////////////////////////////////////
    // Smoother modes.  Objects live for the duration of this RunSolve.
    //////////////////////////////////////////////////////////////////////
    std::unique_ptr<ChebyshevNonHermitianSmoother<LatticeFermionD> > FineCheb;
    std::unique_ptr<ChebyshevNonHermitianSmoother<CoarseVector> >    CoarseCheb;
    std::unique_ptr<GCRReplaySmoother<LatticeFermionD> > FineReplay;
    std::unique_ptr<GCRReplaySmoother<CoarseVector> >    CoarseReplay;
    GCRCoefficients recF, recC;
    if ( FineSmootherMode == "cheb" ) {
      FineCheb.reset(new ChebyshevNonHermitianSmoother<LatticeFermionD>(FineChebLo,FineChebHi,FineSmootherOrder,ShiftedPVdagM));
      FineSmootherSlot.Set(*FineCheb,"Fsmoother Chebyshev");
    }
    if ( CoarseSmootherMode == "cheb" ) {
      CoarseCheb.reset(new ChebyshevNonHermitianSmoother<CoarseVector>(CoarseChebLo,CoarseChebHi,CoarseSmootherNstep,ShiftedC));
      CoarseSmootherSlot.Set(*CoarseCheb,"Csmoother Chebyshev");
    }
    if ( FineSmootherMode   == "replay" ) SmootherGCR.SetCoefficientRecorder(&recF);
    if ( CoarseSmootherMode == "replay" ) CoarseSmootherGCR.SetCoefficientRecorder(&recC);
    L1PGCR.OnStep = [&](int step){
      if ( step != PolyRecordIters ) return;
      if ( FineSmootherMode == "replay" ) {
        SmootherGCR.SetCoefficientRecorder(nullptr);
        recF.Report("Fsmoother");
        FineReplay.reset(new GCRReplaySmoother<LatticeFermionD>(ShiftedPVdagM,recF));
        FineSmootherSlot.Set(*FineReplay,"Fsmoother replay");
      }
      if ( CoarseSmootherMode == "replay" ) {
        CoarseSmootherGCR.SetCoefficientRecorder(nullptr);
        recC.Report("Csmoother");
        CoarseReplay.reset(new GCRReplaySmoother<CoarseVector>(ShiftedC,recC));
        CoarseSmootherSlot.Set(*CoarseReplay,"Csmoother replay");
      }
    };
    std::cout << GridLogMessage << "Smoother modes: fine " << FineSmootherMode << "  coarse " << CoarseSmootherMode
              << (FineSmootherMode=="replay"||CoarseSmootherMode=="replay" ? "  (record for "+std::to_string(PolyRecordIters)+" outer steps)" : "")
              << std::endl;

    std::vector<LatticeFermionD> src(nr,FGrid), sol(nr,FGrid);
    for(int r=0;r<nr;r++){ gaussian(RNG5,src[r]); sol[r]=Zero(); }

    GridStopWatch w; w.Start();
    L1PGCR(src,sol);
    w.Stop();
    std::cout << GridLogMessage << "V2 3-level solve Nrhs "<<nr<<" total " << w.Elapsed()
              << "  (per RHS: " << w.useconds()/1.0e6/nr << " s)" << std::endl;

    // The outer operator is exact by policy; assert the state rather than
    // trust it -- a preconditioner that failed to restore would surface here.
    SetFineSloppy(0);
    { LatticeFermionD Ax(FGrid); RealD worst=0.0;
      for(int r=0;r<nr;r++){ PVdagM.Op(sol[r],Ax); Ax=Ax-src[r];
        RealD rn=std::sqrt(norm2(Ax)/norm2(src[r]));
        std::cout << GridLogMessage << "FINAL Nrhs "<<nr<<": rhs["<<r<<"] true residual = " << rn << std::endl;
        worst=std::max(worst,rn); }
      std::cout << GridLogMessage << "FINAL Nrhs "<<nr<<": worst-case residual = " << worst
                << "   (exact-halo verification)" << std::endl;
    }


    // The operators borrow these grids and build a PaddedCell on them, so
    // they must let go before the grids are destroyed.
    CoarseOpPV.ReleaseGrid();
    CoarseOpL2.ReleaseGrid();
    delete CMrhs; delete CCMrhs;
  };

  RunSolve(nrhs);
  if ( getenv("SOLVE_SRHS")==nullptr || atoi(getenv("SOLVE_SRHS")) ) RunSolve(1);

  std::cout << GridLogMessage << "*** stage three complete: solves done ***" << std::endl;

  Grid_finalize();
}
