/*
 * Example_pvdagm_halfplane.cc
 *
 * Standalone fine-operator diagnostic: the EES half-plane margin of the
 * (non-Hermitian) PV-preconditioned Mobius DWF operator
 *
 *     A(m_adj) = D_adj^dag D_light      (D_adj plays the Pauli-Villars role)
 *
 * as a function of the adjoint mass m_adj, dialled from the light quark mass
 * up to the Pauli-Villars mass (=1).  No coarse grid, no subspace, no Lanczos
 * -- pure power-method spectral tests on the fine grid.
 *
 * Purpose: A is the LEFT preconditioner for inverting the light operator.
 * To solve  D_light X = B  we iterate the preconditioned system
 *     (D_adj^dag D_light) X = D_adj^dag B ,
 * whose solution X is independent of m_adj -- only the conditioning and the
 * iterative convergence change.  m_adj = m_light is the usual CGNR (symmetric
 * normal equations); m_adj = 1 is the Pauli-Villars preconditioned system.
 * The sweep asks which m_adj keeps the preconditioned operator well-behaved
 * (positive-real / EES-guaranteed) while buying the wider spectral range.
 *
 * For the Hermitian part  H(A) = (A + A^dag)/2  we measure, per m_adj:
 *
 *   lambda_max(H)   -- power method on H
 *   lambda_min(H)   -- power method on (sI - H) => min Re W(A), the half-plane
 *                      margin.  EES (Eisenstat-Elman-Schultz 1983, Thm 3.3)
 *                      GUARANTEES GCR convergence with rate
 *                          [ 1 - lambda_min(H)^2 / sigma_max^2 ]^{1/2}
 *                      ONLY when lambda_min(H) > 0 (positive-real / A's field
 *                      of values in the open right half-plane).  A negative
 *                      value means the guarantee is lost (not that GCR
 *                      diverges); the magnitude is then the distance-to-
 *                      positive-realness, i.e. the shift/deflation needed to
 *                      recover it.
 *   sigma_max       -- power method on A^dag A (= A.HermOp)
 *
 * Endpoints:
 *   m_adj = m_light  =>  A = M^dag M, Hermitian PD, positive-real by
 *                        construction, lambda_min(H) = sigma_min^2 > 0 (the
 *                        squared / CGNR operator).
 *   m_adj = 1        =>  A = PV^dag M, the standard PVdagM operator.
 *
 * Env: MASS, M5, MOBIUS_B, MOBIUS_C, LS, CONFIG,
 *      MADJ_LIST (comma separated) OR MADJ_MIN / MADJ_MAX / MADJ_N (geometric).
 *
 * Caveat: lambda_min(H) via a shifted power method can be soft when it sits
 * near zero over a dense low spectrum.  The SIGN and the TREND across m_adj
 * are the robust signal; confirm an individual near-zero value with a proper
 * shifted Lanczos if it is load-bearing.
 */

#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

//////////////////////////////////////////////////////////////////////
// A = PV^dag M : Op = _PV.Mdag . _Mat.M ,  AdjOp = _Mat.Mdag . _PV.M
//////////////////////////////////////////////////////////////////////
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
  void HermOp(const Field &in, Field &out){ // A^dag A
    Field tmp(in.Grid());
    Op(in,tmp);
    AdjOp(tmp,out);
  }
};

//////////////////////////////////////////////////////////////////////
// H = (A + A^dag)/2 for a general non-Hermitian LinearOperator A.
//////////////////////////////////////////////////////////////////////
template<class Field>
class HermitianPartLinOp : public LinearOperatorBase<Field> {
  LinearOperatorBase<Field> &_A;
public:
  HermitianPartLinOp(LinearOperatorBase<Field> &A): _A(A) {};
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
    _A.Op(in,out);      // A in
    _A.AdjOp(in,tmp);   // A^dag in
    out = 0.5*(out + tmp);
  }
};

//////////////////////////////////////////////////////////////////////
// s*I - Op : power method on this gives s - lambda_min(Op) for Hermitian Op.
//////////////////////////////////////////////////////////////////////
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

  RealD mass = 0.00078;
  RealD M5   = 1.8;
  RealD b    = 1.5;
  RealD c    = 0.5;
  int   Ls   = 24;
  std::string config("ckpoint_lat.1000");

  if(getenv("MASS"))     mass = atof(getenv("MASS"));
  if(getenv("M5"))       M5   = atof(getenv("M5"));
  if(getenv("MOBIUS_B")) b    = atof(getenv("MOBIUS_B"));
  if(getenv("MOBIUS_C")) c    = atof(getenv("MOBIUS_C"));
  if(getenv("LS"))       Ls   = atoi(getenv("LS"));
  if(getenv("CONFIG"))   config = std::string(getenv("CONFIG"));

  // Adjoint-mass sweep: explicit list, or geometric MADJ_MIN..MADJ_MAX in MADJ_N steps.
  std::vector<RealD> madj_list;
  if(getenv("MADJ_LIST")){
    std::stringstream ss(getenv("MADJ_LIST"));
    std::string tok;
    while(std::getline(ss,tok,',')) if(tok.size()) madj_list.push_back(std::stod(tok));
  } else {
    int   N  = getenv("MADJ_N")   ? atoi(getenv("MADJ_N"))   : 6;
    RealD lo = getenv("MADJ_MIN") ? atof(getenv("MADJ_MIN")) : mass;
    RealD hi = getenv("MADJ_MAX") ? atof(getenv("MADJ_MAX")) : 1.0;
    GRID_ASSERT(N>=1);
    for(int i=0;i<N;i++)
      madj_list.push_back( (N==1) ? lo : lo*std::pow(hi/lo, double(i)/double(N-1)) );
  }

  // lambda_min(H) is the most-negative eigenvalue; resolved by Chebyshev-filtered
  // Lanczos on H (a shifted power method cannot separate it from the dense low tail).
  RealD HalfChebyLo    = getenv("HALF_CHEBY_LO")    ? atof(getenv("HALF_CHEBY_LO"))    : 0.1;
  RealD HalfChebyHi    = getenv("HALF_CHEBY_HI")    ? atof(getenv("HALF_CHEBY_HI"))    : 0.0; // 0 => auto
  int   HalfChebyOrder = getenv("HALF_CHEBY_ORDER") ? atoi(getenv("HALF_CHEBY_ORDER")) : 61;
  // Grid's Chebyshev filter MUST be odd order: only then is the polynomial positive
  // for x < -1, the region the low/negative modes map to.  An even order flips the
  // sign there, the filtered operator explodes negative, and the IRL never converges.
  if(HalfChebyOrder%2==0){ HalfChebyOrder++;
    std::cout<<GridLogMessage<<"HALF_CHEBY_ORDER forced odd -> "<<HalfChebyOrder<<std::endl; }
  int   HalfNstop      = getenv("HALF_NSTOP")       ? atoi(getenv("HALF_NSTOP"))       : 8;
  int   HalfNk         = getenv("HALF_NK")          ? atoi(getenv("HALF_NK"))          : 24;
  int   HalfNm         = getenv("HALF_NM")          ? atoi(getenv("HALF_NM"))          : 48;
  RealD HalfTol        = getenv("HALF_TOL")         ? atof(getenv("HALF_TOL"))         : 1.0e-4;
  int   HalfMaxIt      = getenv("HALF_MAXIT")       ? atoi(getenv("HALF_MAXIT"))       : 20;

  std::vector<int> lat = {48,48,48,96};

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(lat, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  GridParallelRNG RNG5(FGrid); RNG5.SeedFixedIntegers({5,6,7,8});

  std::cout << GridLogMessage << "PARAM: MASS(light) " << mass << "  M5 " << M5
            << "  b " << b << "  c " << c << "  Ls " << Ls << std::endl;
  std::cout << GridLogMessage << "PARAM: CONFIG " << config << std::endl;

  LatticeGaugeField Umu(UGrid);
  FieldMetaData header;
  std::cout << GridLogMessage << "Reading gauge field " << config << std::endl;
  NerscIO::readConfiguration(Umu,header,config);

  // Fixed light operator (never changes across the sweep).
  MobiusFermionD Dlight(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid, mass, M5, b, c);

  LatticeFermionD x(FGrid);

  std::cout << GridLogMessage << "=================================================" << std::endl;
  std::cout << GridLogMessage << "FINE HALF-PLANE SWEEP  A(m_adj) = D_adj^dag D_light" << std::endl;
  std::cout << GridLogMessage << "  m_adj = " << mass << " => M^dag M (positive-real); m_adj = 1 => PVdagM" << std::endl;
  std::cout << GridLogMessage << "=================================================" << std::endl;

  for(auto madj : madj_list){

    MobiusFermionD Dadj(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid, madj, M5, b, c);

    PVdagMLinearOperator<MobiusFermionD,LatticeFermionD> A(Dlight,Dadj); // A = Dadj^dag Dlight
    HermitianPartLinOp<LatticeFermionD>                  H(A);

    PowerMethod<LatticeFermionD> PM;

    random(RNG5,x); RealD lamHmax = PM(H,x);

    // lambda_min(H): most-negative eigenvalue via Chebyshev-filtered IRL on H.
    // Cheby(lo,hi) amplifies eigenvalues below lo; with hi>=lambda_max(H) the most
    // negative mode is amplified hardest, so IRL isolates the true bottom of the
    // (possibly indefinite) spectrum where the shifted power method could not.
    RealD fhi = (HalfChebyHi>0.0)? HalfChebyHi : 1.1*lamHmax;
    Chebyshev<LatticeFermionD>      Cheby(HalfChebyLo,fhi,HalfChebyOrder);
    FunctionHermOp<LatticeFermionD> OpCheby(Cheby,H);
    PlainHermOp<LatticeFermionD>    OpPlain(H);
    ImplicitlyRestartedLanczos<LatticeFermionD> IRL(OpCheby,OpPlain,HalfNstop,HalfNk,HalfNm,HalfTol,HalfMaxIt);
    std::vector<RealD>           heval(HalfNm);
    std::vector<LatticeFermionD> hevec(HalfNm,FGrid);
    int hNconv=0;
    random(RNG5,x);
    IRL.calc(heval,hevec,x,hNconv);
    RealD lamHmin = (hNconv>0) ? heval[0] : 9.99e99;
    for(int kk=0;kk<hNconv;kk++) lamHmin = std::min(lamHmin, heval[kk]);
    std::cout << GridLogMessage << "  (IRL H-bottom: " << hNconv
              << " converged, most-negative eval " << lamHmin << ")" << std::endl;

    random(RNG5,x); RealD sigmax2 = PM(A,x);   // A.HermOp = A^dag A
    RealD sigmax  = std::sqrt(sigmax2);

    bool  posreal = (lamHmin > 0.0);
    RealD ratefac = posreal ? std::sqrt(1.0 - lamHmin*lamHmin/sigmax2) : 0.0; // EES per-iter
    RealD iters8  = (posreal && ratefac < 1.0) ? std::log(1.0e-8)/std::log(ratefac) : 0.0;

    std::cout << GridLogMessage << "HALFPLANE: m_adj " << madj
              << "  lambda_min(H) " << lamHmin
              << "  lambda_max(H) " << lamHmax
              << "  sigma_max " << sigmax
              << "  positive_real " << (posreal ? "YES" : "NO ")
              << (posreal
                    ? ("  EES_rate " + std::to_string(ratefac) + "  EES_iters(1e-8) " + std::to_string(iters8))
                    : ("  margin_below_zero " + std::to_string(-lamHmin) + " (EES guarantee lost)"))
              << std::endl;
  }

  std::cout << GridLogMessage << "=================================================" << std::endl;
  std::cout << GridLogMessage << "Reading: lambda_min(H) > 0 => EES guarantees GCR at the quoted rate." << std::endl;
  std::cout << GridLogMessage << "         crossing to < 0 as m_adj -> 1 marks loss of positive-realness." << std::endl;
  std::cout << GridLogMessage << "         (non-normality: eigenvalues may still be right-half-plane.)" << std::endl;
  std::cout << GridLogMessage << "Done" << std::endl;

  Grid_finalize();
  return 0;
}
