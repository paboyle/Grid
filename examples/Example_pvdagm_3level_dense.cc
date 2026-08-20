/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Example_pvdagm_3level_dense.cc

    Copyright (C) 2023

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
//////////////////////////////////////////////////////////////////////////////////////
// SUPERCOARSE 3-LEVEL WITH DENSE NON-ITERATIVE COARSE-COARSE INVERSE.
//
// Derived from Example_pvdagm_3level.cc with two changes:
//  1. Maximal second coarsening: Coarse [24,24,24,48] --block {8,4,3,6}--> CC [3,6,8,8]
//     (the geometric floor on --mpi 3.6.4.4; validated to HOLD outer count = 34).
//  2. The CC solve is replaced by an EXPLICIT ROW-DISTRIBUTED DENSE INVERSE
//     (DistributedDenseInverse below): rank = gSites*nbasis = 1152*60 = 69120.
//     Exact, non-iterative: plays the role deflation plays in the Hermitian case.
//
// First-cut design (bulletproof over fast; every step exact or verified):
//   Setup (once per configuration):
//     - PROBE assembly: column j of A_cc = LinOpCC.Op(e_j), unit vectors via global
//       pokeSite; each rank extracts ITS OWN sites' values locally -> assembles its
//       own ROWS of A directly. ~69120 matvecs x ~2.6ms ~ 3-4 min. (Per-trajectory
//       HMC upgrade later: harvest the 33-point stencil blocks directly (~2.2GB) or
//       distance-3 colored probing (81 colors x 60 = 4860 matvecs ~ 13s).)
//     - Gather full A (fp32, 38GB) to boss host by CHUNKED zero-fill+GlobalSum.
//       Exact: each matrix element has exactly ONE nonzero contributor and adding
//       0.0f is exact in IEEE, so the allreduce is bitwise assembly, not arithmetic.
//     - Invert on the boss GCD: rocSOLVER cgetrf+cgetri in fp32 (~tens of seconds).
//       Accuracy: kappa(A_cc) ~ sigma_max/sigma_min ~ 5e3 (census) => fp32 relative
//       error ~ kappa*eps ~ 3e-4, three orders inside the 0.2-tolerance duty this
//       solve replaces. Layout note: we assemble ROW-major and hand the buffer to
//       COLUMN-major LAPACK unchanged; getrf/getri then invert A^T and the col-major
//       result read back row-major is exactly A^{-1} -- the transposes cancel.
//     - Scatter: boss Broadcasts the inverse in chunks; each rank keeps only ITS
//       OWN sites' rows (row ownership == geometric ownership: 4 sites x 60 = 240
//       rows/rank, perfectly uniform on this grid).
//     - VERIFY: || A * (Ainv * x) - x || / ||x|| printed and asserted < 1e-2.
//   Apply (per CC solve; replaces a measured mean 72s iterative solve):
//     - Assemble full x on every rank: zero-fill + GlobalSumVector (1.1 MB, exact).
//     - Local slab GEMV (host, thread_for over my 240 rows, fp32 slab x fp64 x).
//     - pokeLocalSite my own sites. Ownership-aligned rows => NO second global sum.
//     Comms per apply: ONE allreduce of a CC vector. Predicted ~1ms total vs 72s.
//
// Build notes (Frontier/HIP):
//   - Link needs rocSOLVER:  LIBS="-lrocsolver -lrocblas"  (HIP path only).
//   - The boss GCD needs ~38.3GB free HBM for the in-place inversion: run with
//     --device-mem <= 20000 so Grid's pool leaves room (64GB GCD).
//   - Non-HIP builds fall back to Eigen (fine for small local tests; hours at 69k).
//
// Env:  DENSE_CC=1 (default) dense bottom | DENSE_CC=0 iterative L3PGCR + Luscher
//       guesser (the previous supercoarse configuration, for A/B in one binary).
//       CC_TOL (iterative branch tolerance, default 0.05), coarse_smoother_shift,
//       coarse_smoother_nstep, MASS, SUBSPACE_FILE, CC_NEXTRA as before.
//////////////////////////////////////////////////////////////////////////////////////
#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>

#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidual.h>
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidualNonHermitian.h>
#include <Grid/algorithms/iterative/BiCGSTAB.h>

#include <unordered_map>
#include <memory>

#ifdef GRID_HIP
#include <rocsolver/rocsolver.h>
#endif

using namespace std;
using namespace Grid;

template <class T> void readFile(T& out, std::string const fname){
  #ifdef HAVE_LIME
    std::cout << Grid::GridLogMessage << "Reading: " << fname << std::endl;
    Grid::emptyUserRecord record;
    Grid::ScidacReader SR;
    SR.open(fname);
    SR.readScidacFieldRecord(out, record);
    SR.close();
  #endif
}

template <class T> void writeFile(T& in, std::string const fname){
  #ifdef HAVE_LIME
    std::cout << Grid::GridLogMessage << "Writing: " << fname << std::endl;
    Grid::emptyUserRecord record;
    Grid::ScidacWriter SW(in.Grid()->IsBoss());
    SW.open(fname);
    SW.writeScidacFieldRecord(in, record);
    SW.close();
  #endif
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
  int nApp;
  int nAppDag;
public:
  PVdagMLinearOperator(Matrix &Mat,Matrix &PV): _Mat(Mat),_PV(PV), nApp(0), nAppDag(0) {};

  void OpDiag (const Field &in, Field &out) {    assert(0);  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {    assert(0);  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){    assert(0);  };
  void Op     (const Field &in, Field &out){
    Field tmp(in.Grid());
    _Mat.M(in,tmp);
    _PV.Mdag(tmp,out);
    nApp++;
  }
  void AdjOp     (const Field &in, Field &out){
    Field tmp(in.Grid());
    _PV.M(in,tmp);
    _Mat.Mdag(tmp,out);
    nAppDag++;
  }
  void clear() { nApp = 0; nAppDag = 0; }
  void getApplications() {
    std::cout << GridLogMessage << "# applications of PVdagM: " << nApp << std::endl;
    std::cout << GridLogMessage << "# applications of PVdagM^dag: " << nAppDag << std::endl;
    std::cout << GridLogMessage << "# applications total: " << nApp + nAppDag << std::endl;
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
  RealD shift;
public:
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
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){ assert(0); }
  void HermOp(const Field &in, Field &out){
    Field tmp(in.Grid());
    Op(in,tmp);
    AdjOp(tmp,out);
  }
};

// Lüscher deflated guesser (arXiv:0706.2298 Sec A.3) for a non-Hermitian solve.
// C_{st} = <psi[s] | LinOp | psi[t]>;  guess = sum_s c_s psi[s]  where c = C^{-1} psi† src.
template<class Field>
class LuscherGuesser : public LinearFunction<Field> {
  const std::vector<Field> &psi;
  Eigen::MatrixXcd          C_inv;
public:
  using LinearFunction<Field>::operator();
  LuscherGuesser(const std::vector<Field> &psi_, const Eigen::MatrixXcd &Cinv_)
    : psi(psi_), C_inv(Cinv_) {}
  virtual void operator()(const Field &src, Field &guess) {
    int N = psi.size();
    Eigen::VectorXcd b(N);
    for (int t = 0; t < N; t++)
      b(t) = TensorRemove(innerProduct(psi[t], src));
    Eigen::VectorXcd c = C_inv * b;
    guess = Zero();
    for (int s = 0; s < N; s++)
      guess += ComplexD(c(s)) * psi[s];
  }
};

//////////////////////////////////////////////////////////////////////////////////////
// DistributedDenseInverse: exact solve of a (small) coarse operator by explicit,
// row-distributed dense inverse.  See file header for the full design rationale.
//
//   Row ownership == geometric ownership: rank owns rows (site,b) for its own local
//   sites => after the local GEMV the output already lives where it is needed and
//   the apply needs only ONE GlobalSum (assembly of x), not two.
//
//   All tensor-nest depth issues are dodged by treating the site scalar_object as a
//   contiguous array of nbasis ComplexD (static-asserted): works at any MG depth.
//////////////////////////////////////////////////////////////////////////////////////
template<class Field>
class DistributedDenseInverse : public LinearFunction<Field> {
public:
  using LinearFunction<Field>::operator();
  typedef typename Field::vector_object vobj;
  typedef typename vobj::scalar_object  sobj;

  GridBase *grid;
  LinearOperatorBase<Field> &_Op;   // kept for env-gated defect checks (DENSE_CC_CHECK)
  int      nbasis;
  int64_t  N;               // dense rank = gSites * nbasis
  int      nd;
  int      lsites;          // my local sites
  int      nrows;           // my rows = lsites * nbasis
  std::vector<Coordinate> myLcoor;      // local coordinate of my site ss
  std::vector<int64_t>    myGsite;      // global lex site index of my site ss
  std::vector<ComplexF>   slab;         // nrows x N row-major; A rows during probe,
                                        // A^{-1} rows after setup
  static const int64_t CHUNKROWS = 256; // gather/broadcast chunk (256 x N x 8B ~ 141MB)

  DistributedDenseInverse(LinearOperatorBase<Field> &Op, GridBase *g, int nbasis_)
    : grid(g), _Op(Op), nbasis(nbasis_)
  {
    GRID_ASSERT( sizeof(sobj) == nbasis*sizeof(ComplexD) ); // site object == nbasis ComplexD
    nd     = grid->_ndimension;
    N      = grid->gSites() * nbasis;
    lsites = grid->lSites();
    nrows  = lsites * nbasis;

    Coordinate ldims = grid->LocalDimensions();
    Coordinate gdims = grid->GlobalDimensions();

    std::cout << GridLogMessage << "DistributedDenseInverse: N = " << N
              << " (" << grid->gSites() << " sites x " << nbasis << ")"
              << "  rows/rank = " << nrows
              << "  slab = " << (double)nrows*N*sizeof(ComplexF)/1024./1024. << " MB/rank"
              << std::endl;

    ////////////////////////////////////////////////////////////////////
    // Enumerate my sites: local coords and global lexicographic indices
    ////////////////////////////////////////////////////////////////////
    myLcoor.resize(lsites);
    myGsite.resize(lsites);
    for(int ss=0; ss<lsites; ss++){
      Coordinate lcoor(nd);
      Lexicographic::CoorFromIndex(lcoor, ss, ldims);
      Coordinate gcoor(nd);
      for(int d=0; d<nd; d++) gcoor[d] = grid->_lstart[d] + lcoor[d];
      int64_t gsite;
      Lexicographic::IndexFromCoor(gcoor, gsite, gdims);
      myLcoor[ss] = lcoor;
      myGsite[ss] = gsite;
    }

    slab.resize((uint64_t)nrows * N);

    ////////////////////////////////////////////////////////////////////
    // 1. PROBE assembly of my rows of A: column (jsite,b) = Op(e_{jsite,b})
    //    (bulletproof: uses only the public operator interface; upgrade
    //     path for HMC = stencil harvest or colored probing, see header)
    ////////////////////////////////////////////////////////////////////
    double t0 = usecond();
    Field e(grid);
    Field Ae(grid);
    Coordinate gcoorj(nd);
    int64_t gsitesN = grid->gSites();
    for(int64_t jsite=0; jsite<gsitesN; jsite++){
      Lexicographic::CoorFromIndex(gcoorj, jsite, gdims);
      for(int b=0; b<nbasis; b++){
        int64_t col = jsite*nbasis + b;
        sobj s = Zero();
        ((ComplexD *)&s)[b] = ComplexD(1.0,0.0);
        e = Zero();
        pokeSite(s, e, gcoorj);
        Op.Op(e, Ae);
        for(int ss=0; ss<lsites; ss++){
          sobj t;
          peekLocalSite(t, Ae, myLcoor[ss]);
          for(int a=0; a<nbasis; a++){
            slab[ (uint64_t)(ss*nbasis+a)*N + col ] = ComplexF( ((ComplexD *)&t)[a] );
          }
        }
      }
      if ( (jsite % (gsitesN/20+1)) == 0 ) {
        std::cout << GridLogMessage << "DistributedDenseInverse: probed site "
                  << jsite << " / " << gsitesN << std::endl;
      }
    }
    double t1 = usecond();
    std::cout << GridLogMessage << "DistributedDenseInverse: probe assembly took "
              << (t1-t0)/1.0e6 << " s" << std::endl;

    ////////////////////////////////////////////////////////////////////
    // 2. Gather full A (fp32) by chunked zero-fill+GlobalSum.
    //    Exact: one contributor per element; adding 0.0f is exact.
    //    HIP path: chunks stream STRAIGHT INTO THE DEVICE inversion buffer --
    //    the boss host never holds the 38GB matrix (a host Afull on top of the
    //    ~100GB boss-rank RSS OOMed the boss node on the first attempt).
    //    Non-HIP (Eigen) path keeps the host Afull.
    ////////////////////////////////////////////////////////////////////
    int boss = grid->IsBoss();
    std::vector<ComplexF> Afull;
#ifdef GRID_HIP
    rocblas_float_complex *dA = nullptr;
    rocblas_float_complex *dB = nullptr;   // getrs_64 RHS block (CHUNKROWS identity columns)
    rocblas_handle rochandle;
    int64_t *dIpiv = nullptr;              // ILP64 pivots, live from factor to last getrs
    uint64_t Abytes = (uint64_t)N * N * sizeof(ComplexF);
    // Make HBM space for the naked hipMalloc: flush BOTH retention layers, in
    // this order -- EvictAll() pushes device copies into the allocator's
    // deferred-free pool, FreePool() then returns that pool to the driver.
    // (PB's MemoryManager additions; keeps --device-mem and cache semantics
    // intact for normal operation. Boss-only would suffice; all-ranks is cheap.)
    MemoryManager::EvictAll();
    // MemoryManager::FreePool();  // not needed at NB_CC=30 (8.9GB buffer; EvictAll
    // suffices); re-enable with the type-dispatched FreePool once the corrected
    // MemoryManager lands (blanket acceleratorFreeDevice poisoned HIP runtime).
    if (boss) {
      auto aerr = hipMalloc((void **)&dA, Abytes);
      if (aerr != hipSuccess) {
        std::cout << GridLogMessage << "DistributedDenseInverse: hipMalloc of "
                  << Abytes/1024./1024./1024. << " GB FAILED -- reduce --device-mem "
                  << "(need ~36GB free on the boss GCD; try --device-mem 16000)" << std::endl;
        GRID_ASSERT(aerr == hipSuccess);
      }
      std::cout << GridLogMessage << "DistributedDenseInverse: device inversion buffer allocated ("
                << Abytes/1024./1024./1024. << " GB)" << std::endl;
    }
#else
    if (boss) Afull.resize((uint64_t)N * N);
#endif
    {
      std::unordered_map<int64_t,int> rowmap;   // global row -> my slab row
      for(int ss=0; ss<lsites; ss++)
        for(int a=0; a<nbasis; a++)
          rowmap[ myGsite[ss]*nbasis + a ] = ss*nbasis + a;

      std::vector<ComplexF> chunk((uint64_t)CHUNKROWS * N);
      for(int64_t row0=0; row0<N; row0+=CHUNKROWS){
        int64_t nrow = std::min(CHUNKROWS, N-row0);
        uint64_t nelem = (uint64_t)nrow * N;
        for(uint64_t i=0;i<nelem;i++) chunk[i]=ComplexF(0.0,0.0);
        for(int64_t r=row0; r<row0+nrow; r++){
          auto it = rowmap.find(r);
          if (it != rowmap.end()) {
            uint64_t src = (uint64_t)(it->second) * N;
            uint64_t dst = (uint64_t)(r-row0) * N;
            for(int64_t j=0;j<N;j++) chunk[dst+j] = slab[src+j];
          }
        }
        grid->GlobalSumVector(&chunk[0], (int)nelem);
        if (boss) {
#ifdef GRID_HIP
          GRID_ASSERT( hipMemcpy((char *)dA + (uint64_t)row0*N*sizeof(ComplexF),
                                 &chunk[0], nelem*sizeof(ComplexF),
                                 hipMemcpyHostToDevice) == hipSuccess );
#else
          uint64_t dst = (uint64_t)row0 * N;
          for(uint64_t i=0;i<nelem;i++) Afull[dst+i] = chunk[i];
#endif
        }
      }
    }
    double t2 = usecond();
    std::cout << GridLogMessage << "DistributedDenseInverse: gather to boss took "
              << (t2-t1)/1.0e6 << " s" << std::endl;

    ////////////////////////////////////////////////////////////////////
    // 3. Invert in place on the boss (fp32).
    //    Row-major buffer handed to column-major LAPACK => it inverts A^T,
    //    whose col-major result read back row-major is exactly A^{-1}.
    ////////////////////////////////////////////////////////////////////
    if (boss) {
#ifdef GRID_HIP
      std::cout << GridLogMessage << "DistributedDenseInverse: rocSOLVER cgetrf_64 (ILP64 LU) N=" << N
                << " in place on resident device buffer" << std::endl;
      auto hst = rocblas_create_handle(&rochandle);
      std::cout << GridLogMessage << "DistributedDenseInverse: rocblas handle status " << (int)hst << std::endl;
      int64_t *dInfo;
      GRID_ASSERT( hipMalloc((void **)&dIpiv, N*sizeof(int64_t)) == hipSuccess );
      GRID_ASSERT( hipMalloc((void **)&dInfo, sizeof(int64_t))   == hipSuccess );
      auto st1 = rocsolver_cgetrf_64(rochandle, (int64_t)N, (int64_t)N, dA, (int64_t)N, dIpiv, dInfo);
      hipDeviceSynchronize();
      int64_t info_h = -1;
      hipMemcpy(&info_h, dInfo, sizeof(int64_t), hipMemcpyDeviceToHost);
      std::cout << GridLogMessage << "DistributedDenseInverse: cgetrf_64 status " << (int)st1
                << " info = " << (int)info_h << std::endl;
      GRID_ASSERT(st1 == rocblas_status_success);
      GRID_ASSERT(info_h == 0);
      hipFree(dInfo);
      GRID_ASSERT( hipMalloc((void **)&dB, (uint64_t)CHUNKROWS*N*sizeof(ComplexF)) == hipSuccess );
      // dA holds the LU of A^T (row-major assembly == col-major A^T); rows of
      // A^{-1} are produced blockwise in the scatter loop via cgetrs_64 on
      // identity-column blocks: A^T X = E  =>  X columns = rows of A^{-1},
      // in exactly the linear layout the existing harvest expects.
#else
      // Eigen fallback: intended for small local CPU tests; hours at N ~ 69k.
      std::cout << GridLogMessage << "DistributedDenseInverse: Eigen fallback inversion N=" << N
                << (N > 10000 ? "  (WARNING: this will be SLOW; use the HIP/rocSOLVER path)" : "")
                << std::endl;
      typedef Eigen::Matrix<std::complex<float>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> MatF;
      Eigen::Map<MatF> A(reinterpret_cast<std::complex<float>*>(&Afull[0]), N, N);
      MatF Ainv = A.inverse();
      A = Ainv;
#endif
    }
    double t3 = usecond();
    std::cout << GridLogMessage << "DistributedDenseInverse: inversion took "
              << (t3-t2)/1.0e6 << " s" << std::endl;

    ////////////////////////////////////////////////////////////////////
    // 4. Scatter: chunked Broadcast of A^{-1}; each rank keeps its rows.
    ////////////////////////////////////////////////////////////////////
    {
      std::unordered_map<int64_t,int> rowmap;
      for(int ss=0; ss<lsites; ss++)
        for(int a=0; a<nbasis; a++)
          rowmap[ myGsite[ss]*nbasis + a ] = ss*nbasis + a;

      std::vector<ComplexF> chunk((uint64_t)CHUNKROWS * N);
      for(int64_t row0=0; row0<N; row0+=CHUNKROWS){
        int64_t nrow = std::min(CHUNKROWS, N-row0);
        uint64_t nelem = (uint64_t)nrow * N;
        if (boss) {
#ifdef GRID_HIP
          // Identity block E: column j = e_{row0+j}; solve A^T X = E so X's
          // columns are rows [row0,row0+nrow) of A^{-1} (see factor comment).
          for(uint64_t i=0;i<nelem;i++) chunk[i] = ComplexF(0.0,0.0);
          for(int64_t j=0;j<nrow;j++) chunk[(uint64_t)j*N + (uint64_t)(row0+j)] = ComplexF(1.0,0.0);
          GRID_ASSERT( hipMemcpy(dB, &chunk[0], nelem*sizeof(ComplexF), hipMemcpyHostToDevice) == hipSuccess );
          auto strs = rocsolver_cgetrs_64(rochandle, rocblas_operation_none,
                                          (int64_t)N, (int64_t)nrow,
                                          dA, (int64_t)N, dIpiv, dB, (int64_t)N);
          GRID_ASSERT(strs == rocblas_status_success);
          hipDeviceSynchronize();
          GRID_ASSERT( hipMemcpy(&chunk[0], dB, nelem*sizeof(ComplexF), hipMemcpyDeviceToHost) == hipSuccess );
#else
          uint64_t src = (uint64_t)row0 * N;
          for(uint64_t i=0;i<nelem;i++) chunk[i] = Afull[src+i];
#endif
        }
        grid->Broadcast(0, &chunk[0], nelem*sizeof(ComplexF));
        for(int64_t r=row0; r<row0+nrow; r++){
          auto it = rowmap.find(r);
          if (it != rowmap.end()) {
            uint64_t dst = (uint64_t)(it->second) * N;
            uint64_t src = (uint64_t)(r-row0) * N;
            for(int64_t j=0;j<N;j++) slab[dst+j] = chunk[src+j];
          }
        }
      }
    }
#ifdef GRID_HIP
    if (boss) {
      if (dA)    hipFree(dA);
      if (dB)    hipFree(dB);
      if (dIpiv) hipFree(dIpiv);
      rocblas_destroy_handle(rochandle);
    }
#endif
    double t4 = usecond();
    std::cout << GridLogMessage << "DistributedDenseInverse: blocked getrs solve+scatter of inverse took "
              << (t4-t3)/1.0e6 << " s" << std::endl;

    ////////////////////////////////////////////////////////////////////
    // 5. VERIFY: || A (Ainv x) - x || / ||x|| on a deterministic test vector,
    //    plus one timed apply. Catches indexing, layout, transpose and
    //    precision errors in a single number. Expect ~1e-4 (fp32, kappa~5e3).
    ////////////////////////////////////////////////////////////////////
    {
      Field x(grid); Field y(grid); Field z(grid);
      x = ComplexD(1.0,0.0);
      double ta = usecond();
      (*this)(x, y);
      double tb = usecond();
      Op.Op(y, z);
      z = z - x;
      RealD rel = std::sqrt(norm2(z)/norm2(x));
      std::cout << GridLogMessage << "DistributedDenseInverse: VERIFY ||A Ainv x - x||/||x|| = "
                << rel << "   (one apply took " << (tb-ta)/1000.0 << " ms)" << std::endl;
      GRID_ASSERT(rel < 1.0e-2);
    }
    std::cout << GridLogMessage << "DistributedDenseInverse: setup complete, total "
              << (usecond()-t0)/1.0e6 << " s" << std::endl;
  }

  ////////////////////////////////////////////////////////////////////
  // Apply: psi = A^{-1} src.
  //   ONE GlobalSum (assemble x everywhere; exact), local slab GEMV,
  //   poke my own sites (ownership-aligned rows: no output comms).
  ////////////////////////////////////////////////////////////////////
  virtual void operator()(const Field &src, Field &psi)
  {
    std::vector<ComplexD> x((uint64_t)N, ComplexD(0.0,0.0));
    for(int ss=0; ss<lsites; ss++){
      sobj s;
      peekLocalSite(s, src, myLcoor[ss]);
      for(int b=0; b<nbasis; b++)
        x[ myGsite[ss]*nbasis + b ] = ((ComplexD *)&s)[b];
    }
    grid->GlobalSumVector(&x[0], (int)N);

    std::vector<ComplexD> y(nrows);
    thread_for(r, nrows, {
      ComplexD acc(0.0,0.0);
      const ComplexF *row = &slab[(uint64_t)r * N];
      for(int64_t j=0; j<N; j++) acc += ComplexD(row[j]) * x[j];
      y[r] = acc;
    });

    for(int ss=0; ss<lsites; ss++){
      sobj s;
      for(int b=0; b<nbasis; b++)
        ((ComplexD *)&s)[b] = y[ss*nbasis + b];
      pokeLocalSite(s, psi, myLcoor[ss]);
    }

    // Env-gated per-apply defect: one CC matvec (~few ms) restores the
    // true-residual visibility an iterative bottom would have printed.
    if ( getenv("DENSE_CC_CHECK") ) {
      Field tmp(grid);
      _Op.Op(psi, tmp);
      tmp = tmp - src;
      std::cout << GridLogMessage << "DistributedDenseInverse: apply defect ||A x - b||/||b|| = "
                << std::sqrt(norm2(tmp)/norm2(src)) << std::endl;
    }
  }
};

template<class Fobj,class CComplex,int nbasis>
class MGPreconditioner : public LinearFunction< Lattice<Fobj> > {
public:
  using LinearFunction<Lattice<Fobj> >::operator();

  typedef Aggregation<Fobj,CComplex,nbasis> Aggregates;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::FineField    FineField;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::CoarseVector CoarseVector;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::CoarseMatrix CoarseMatrix;
  typedef LinearOperatorBase<FineField>                            FineOperator;
  typedef LinearFunction    <FineField>                            FineSmoother;
  typedef LinearOperatorBase<CoarseVector>                         CoarseOperator;
  typedef LinearFunction    <CoarseVector>                         CoarseSolver;

  Aggregates     & _Aggregates;
  FineOperator   & _FineOperator;
  FineSmoother   & _PreSmoother;
  FineSmoother   & _PostSmoother;
  CoarseOperator & _CoarseOperator;
  CoarseSolver   & _CoarseSolve;
  CoarseSolver   & _CoarseGuesser;

  int    level;  void Level(int lv) {level = lv; };

  MGPreconditioner(Aggregates &Agg,
		   FineOperator &Fine,
		   FineSmoother &PreSmoother,
		   FineSmoother &PostSmoother,
		   CoarseOperator &CoarseOperator_,
		   CoarseSolver &CoarseSolve_,
		   CoarseSolver &CoarseGuesser_)
    : _Aggregates(Agg),
      _FineOperator(Fine),
      _PreSmoother(PreSmoother),
      _PostSmoother(PostSmoother),
      _CoarseOperator(CoarseOperator_),
      _CoarseSolve(CoarseSolve_),
      _CoarseGuesser(CoarseGuesser_),
      level(1)  {  }

  virtual void operator()(const FineField &in, FineField & out)
  {
    GridBase *CoarseGrid = _Aggregates.CoarseGrid;
    CoarseVector Csrc(CoarseGrid);
    CoarseVector Csol(CoarseGrid);
    FineField vec1(in.Grid());
    FineField vec2(in.Grid());

    double t;
    out = Zero();
    t=-usecond();
    _PreSmoother(in,out);
    t+=usecond();
    std::cout<<GridLogMessage << "PreSmoother took "<< t/1000.0<< "ms" <<std::endl;

    _FineOperator.Op(out,vec1);  sub(vec1, in ,vec1);

    t=-usecond();
    _Aggregates.ProjectToSubspace(Csrc,vec1);
    t+=usecond();
    std::cout<<GridLogMessage << "Project to coarse took "<< t/1000.0<< "ms" <<std::endl;

    t=-usecond();
    _CoarseGuesser(Csrc,Csol);
    _CoarseSolve(Csrc,Csol);
    t+=usecond();
    std::cout<<GridLogMessage << "Coarse solve took "<< t/1000.0<< "ms" <<std::endl;

    t=-usecond();
    _Aggregates.PromoteFromSubspace(Csol,vec1);
    add(out,out,vec1);
    t+=usecond();
    std::cout<<GridLogMessage << "Promote to this level took "<< t/1000.0<< "ms" <<std::endl;

    _FineOperator.Op(out,vec1);  sub(vec1 ,in , vec1);

    t=-usecond();
    vec2=Zero();
    _PostSmoother(vec1,vec2);
    t+=usecond();
    std::cout<<GridLogMessage << "PostSmoother took "<< t/1000.0<< "ms" <<std::endl;

    add(out,out,vec2);
  }
};

// Generic shifted linear operator: wraps any LinearOperatorBase and adds shift*I.
// Used to condition the coarse-level GCR smoother, analogous to ShiftedPVdagMLinearOperator
// at the fine level.
template<class Field>
class ShiftedLinearOperator : public LinearOperatorBase<Field> {
  LinearOperatorBase<Field> &_Op;
  RealD shift;
public:
  ShiftedLinearOperator(RealD _shift, LinearOperatorBase<Field> &Op) : shift(_shift), _Op(Op) {}
  void OpDiag   (const Field &in, Field &out)                         { assert(0); }
  void OpDir    (const Field &in, Field &out, int dir, int disp)      { assert(0); }
  void OpDirAll (const Field &in, std::vector<Field> &out)            { assert(0); }
  void Op       (const Field &in, Field &out) { _Op.Op(in, out);    out = out + shift * in; }
  void AdjOp    (const Field &in, Field &out) { _Op.AdjOp(in, out); out = out + shift * in; }
  void HermOpAndNorm(const Field &in, Field &out, RealD &n1, RealD &n2) { assert(0); }
  void HermOp   (const Field &in, Field &out) { Field tmp(in.Grid()); Op(in,tmp); AdjOp(tmp,out); }
};

template<int NB, class PVdagM_t, class ShiftedPVdagM_t, class Subspace, class LittleDiracOperator, class CoarseVector, class TwoLevelMG>
void runMG(
  GridCartesian *FGrid,
  GridCartesian *Coarse5d,
  GridCartesian *CoarseCoarse5d,
  NextToNearestStencilGeometry5D geom,
  PVdagM_t &PVdagM,
  ShiftedPVdagM_t &ShiftedPVdagM,
  Subspace &AggregatesPD
) {
  std::vector<LatticeFermion> subspace = AggregatesPD.subspace;
  assert((int)subspace.size() == NB);
  const int nbasis = NB;
  const int cb = 0;

  CoarseVector c_src(Coarse5d);
  CoarseVector c_res(Coarse5d);
  Complex one(1.0);

  LatticeFermionD f_src(FGrid);
  LatticeFermionD f_res(FGrid);

  TrivialPrecon<CoarseVector>      simpleC;
  TrivialPrecon<LatticeFermionD>   simple_fine;

  //////////////////////////////////////////////////////////////////////
  // Level 0→1: coarsen PVdagM, build LinOpCoarse
  //////////////////////////////////////////////////////////////////////
  LittleDiracOperator LittleDiracOpPV(geom, FGrid, Coarse5d);
  LittleDiracOpPV.CoarsenOperator(PVdagM, AggregatesPD);

  NonHermitianLinearOperator<LittleDiracOperator,CoarseVector> LinOpCoarse(LittleDiracOpPV);

  //////////////////////////////////////////////////////////////////////
  // Baseline: plain PGCR on LinOpCoarse (reference for comparison)
  //////////////////////////////////////////////////////////////////////
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout<<GridLogMessage<<" Level 1 solve: plain PGCR baseline"<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;

  PrecGeneralisedConjugateResidualNonHermitian<CoarseVector> L2PGCR_baseline(3.0e-2,1100,LinOpCoarse,simpleC,10,10);
  L2PGCR_baseline.Level(2);
  L2PGCR_baseline.Name("Cbaseline");
  c_src = one;
  c_res = Zero();
  L2PGCR_baseline(c_src,c_res);

  //////////////////////////////////////////////////////////////////////
  // psi_coarse: coarse projections of pre-GS fine null vectors.
  //////////////////////////////////////////////////////////////////////
  std::vector<CoarseVector> psi_coarse(nbasis, Coarse5d);
  for (int k = 0; k < nbasis; k++)
    AggregatesPD.ProjectToSubspace(psi_coarse[k], subspace[k]);

  //////////////////////////////////////////////////////////////////////
  // Optional sigma-ordering of psi_coarse (SVD_REORDER set): replace the crude
  // first-NB_CC slice with the NB_CC most-null directions of span(psi_coarse)
  // under LinOpCoarse.  Nullness measure = singular values (eig of the
  // Gram-whitened Psi†A†APsi), NOT the numerical range Q†AQ which
  // non-normality contaminates.  The printed sigma spectrum shows where the
  // truncation cliff sits.  Unset => raw first-30 (crude GS-ordered slice).
  //////////////////////////////////////////////////////////////////////
  if ( getenv("SVD_REORDER") ) {
    std::cout << GridLogMessage << "SVD_REORDER: sigma-ordering psi_coarse under LinOpCoarse" << std::endl;

    Eigen::MatrixXcd G(nbasis,nbasis);                    // Gram  = Psi^dag Psi
    for (int i=0;i<nbasis;i++)
    for (int j=0;j<nbasis;j++)
      G(i,j) = TensorRemove(innerProduct(psi_coarse[i],psi_coarse[j]));

    std::vector<CoarseVector> Apsi(nbasis, Coarse5d);
    for (int j=0;j<nbasis;j++) LinOpCoarse.Op(psi_coarse[j], Apsi[j]);

    Eigen::MatrixXcd M(nbasis,nbasis);                    // Psi^dag A^dag A Psi
    for (int i=0;i<nbasis;i++)
    for (int j=0;j<nbasis;j++)
      M(i,j) = TensorRemove(innerProduct(Apsi[i],Apsi[j]));

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> esG(G);
    Eigen::VectorXd g = esG.eigenvalues();               // ascending, real
    RealD gmax = g(nbasis-1);
    RealD gtol = 1.0e-9 * gmax;
    int keep = 0; for (int i=0;i<nbasis;i++) if (g(i) > gtol) keep++;
    std::cout << GridLogMessage << "  Gram spectrum: min=" << g(0) << " max=" << gmax
              << " cond=" << gmax/std::max(g(0),1.0e-300) << " keep=" << keep << "/" << nbasis << std::endl;

    Eigen::MatrixXcd T(nbasis, keep);                    // whitening (largest-g first)
    { int c=0;
      for (int i=nbasis-1;i>=0;i--) if (g(i) > gtol) { T.col(c) = esG.eigenvectors().col(i)/std::sqrt(g(i)); c++; }
    }
    Eigen::MatrixXcd Mw = T.adjoint() * M * T;           // whitened A^dagA (Hermitian)
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> esM(Mw);
    Eigen::VectorXd s2 = esM.eigenvalues();              // ascending sigma^2: most-null first
    std::cout << GridLogMessage << "  Singular spectrum sigma_k (most-null first):" << std::endl;
    for (int k=0;k<keep;k++)
      std::cout << GridLogMessage << "    sigma[" << k << "] = " << std::sqrt(std::max(s2(k),0.0)) << std::endl;

    Eigen::MatrixXcd R = T * esM.eigenvectors();         // coeffs over Psi, sigma-ordered
    std::vector<CoarseVector> phi(keep, Coarse5d);
    for (int k=0;k<keep;k++) {
      phi[k] = Zero();
      for (int j=0;j<nbasis;j++)
        phi[k] = phi[k] + ComplexD(R(j,k)) * psi_coarse[j];
    }
    for (int k=0;k<keep;k++) psi_coarse[k] = phi[k];     // [0..NB_CC-1] = most-null dirs
    std::cout << GridLogMessage << "SVD_REORDER: psi_coarse replaced by sigma-ordered directions" << std::endl;
  }

  //////////////////////////////////////////////////////////////////////
  // Level 1→2: SUPERCOARSE aggregation using psi_coarse as subspace.
  // Maximal block {8,4,3,6}: CC = [3,6,8,8] = the dense-invertible floor.
  //
  // UNBLOCKING TRUNCATION NB_CC = 30 (first-30 slice of psi_coarse):
  //   N = 1152*30 = 34560, N^2 = 1.19e9 < 2^31 => stock rocSOLVER cgetrf/cgetri
  //   is int32-safe.  At NB=60 (N=69120, N^2=4.78e9) rocSOLVER faults on
  //   32-bit internal indexing (GPU memory access fault inside cgetrf; ROCm
  //   7.2.0 has rocsolver_cgetrf_64 but NO cgetri_64).  Restore NB_CC=60 when
  //   the getrf_64 + blocked-getrs (or tiled) inversion lands.
  //   NOTE: 30-basis CC operator = weaker coarse correction; outer count may
  //   drift from 34 -- this run validates the dense machinery end-to-end, and
  //   doubles as the bottom-thinning experiment (5-level held at 30).
  //////////////////////////////////////////////////////////////////////
  // NB_CC = 60: the sigma spectrum (denseCC.log.5110401) is FLAT (0.0255..0.0327,
  // factor 1.29 over all 60) => every promoted direction is equally near-null;
  // truncation to 30 costs ~51->91 L2 iterations and NO ordering can recover it.
  // N = 69120 => rocSOLVER int32 getri faults; use cgetrf_64 + blocked cgetrs_64
  // (ILP64, present in ROCm 7.2) producing rows of A^{-1} blockwise.
  constexpr int NB_CC = 60;
  std::cout << GridLogMessage << "PARAM NB_CC (CC basis) = " << NB_CC
            << "  (dense rank N = " << CoarseCoarse5d->gSites()*NB_CC << ")" << std::endl;

  typedef typename CoarseVector::vector_object                            CoarseSiteObj;
  typedef iScalar<vTComplex>                                              vTTComplex;
  typedef GeneralCoarsenedMatrix<CoarseSiteObj,vTTComplex,NB_CC>         LittleDiracOperatorL2;
  typedef typename LittleDiracOperatorL2::CoarseVector                   CoarseCoarseVector;
  typedef Aggregation<CoarseSiteObj,vTTComplex,NB_CC>                    SubspaceL2;
  typedef MGPreconditioner<CoarseSiteObj,vTTComplex,NB_CC>               L1to2MG;

  SubspaceL2 AggregatesL2(CoarseCoarse5d, Coarse5d, cb);
  for (int k = 0; k < NB_CC; k++)
    AggregatesL2.subspace[k] = psi_coarse[k];      // first NB_CC raw promoted vectors

  NextToNearestStencilGeometry5D geom2(CoarseCoarse5d);
  LittleDiracOperatorL2 LittleDiracOpL2(geom2, Coarse5d, CoarseCoarse5d);
  LittleDiracOpL2.CoarsenOperator(LinOpCoarse, AggregatesL2);

  NonHermitianLinearOperator<LittleDiracOperatorL2,CoarseCoarseVector> LinOpCC(LittleDiracOpL2);

  TrivialPrecon<CoarseCoarseVector> simpleCC;

  //////////////////////////////////////////////////////////////////////
  // CC solve: DENSE (exact, non-iterative) by default; DENSE_CC=0 gives
  // the previous iterative L3PGCR + Lüscher-guesser path for A/B.
  //////////////////////////////////////////////////////////////////////
  int use_dense = 1;
  if (getenv("DENSE_CC")) use_dense = atoi(getenv("DENSE_CC"));
  RealD cc_tol = 0.05;
  if (getenv("CC_TOL")) cc_tol = atof(getenv("CC_TOL"));
  std::cout << GridLogMessage << "PARAM DENSE_CC = " << use_dense
            << (use_dense ? "  (dense exact CC inverse)" : "  (iterative CC solve)") << std::endl;

  std::unique_ptr<DistributedDenseInverse<CoarseCoarseVector>> DenseCC;
  std::unique_ptr<LuscherGuesser<CoarseCoarseVector>>          CCDeflGuesser;
  std::vector<CoarseCoarseVector> psi_cc;   // kept alive for the guesser branch

  PrecGeneralisedConjugateResidualNonHermitian<CoarseCoarseVector> L3PGCR(cc_tol,200,LinOpCC,simpleCC,16,16);
  L3PGCR.Level(3);
  L3PGCR.Name("CCouter");

  LinearFunction<CoarseCoarseVector> *ccSolve;
  LinearFunction<CoarseCoarseVector> *ccGuess;

  if (use_dense) {
    std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
    std::cout<<GridLogMessage<<" Dense CC inverse setup"<<std::endl;
    std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
    DenseCC.reset(new DistributedDenseInverse<CoarseCoarseVector>(LinOpCC, CoarseCoarse5d, NB_CC));
    ccSolve = DenseCC.get();
    ccGuess = &simpleCC;      // exact solve ignores/overwrites any guess
  } else {
    ////////////////////////////////////////////////////////////////////
    // Lüscher deflation guesser (arXiv:0706.2298 A.3) for the iterative CC
    // solve, as in the earlier supercoarse configuration.
    ////////////////////////////////////////////////////////////////////
    psi_cc.resize(nbasis, CoarseCoarse5d);
    for (int k = 0; k < nbasis; k++)
      AggregatesL2.ProjectToSubspace(psi_cc[k], psi_coarse[k]);
    {
      int Nextra = nbasis;
      if ( getenv("CC_NEXTRA") ) Nextra = atoi(getenv("CC_NEXTRA"));
      GridParallelRNG RNG_CC(CoarseCoarse5d);
      RNG_CC.SeedFixedIntegers({11,13,17,19});
      PrecGeneralisedConjugateResidualNonHermitian<CoarseCoarseVector>
        nullGCR(1e-2, 200, LinOpCC, simpleCC, 32, 32);
      CoarseCoarseVector tmp(CoarseCoarse5d);
      for (int k = 0; k < Nextra; k++) {
        CoarseCoarseVector src(CoarseCoarse5d);
        gaussian(RNG_CC, src);
        tmp = Zero();
        nullGCR(src, tmp);
        psi_cc.push_back(tmp);
      }
      std::cout << GridLogMessage << "LinOpCC deflation basis: " << nbasis
                << " promoted + " << Nextra << " bred = " << psi_cc.size() << " total" << std::endl;
    }
    const int Naug = psi_cc.size();
    Eigen::MatrixXcd Ccc = Eigen::MatrixXcd::Zero(Naug, Naug);
    {
      CoarseCoarseVector Acc(CoarseCoarse5d);
      for (int l = 0; l < Naug; l++) {
        LinOpCC.Op(psi_cc[l], Acc);
        for (int k = 0; k < Naug; k++)
          Ccc(k,l) = TensorRemove(innerProduct(psi_cc[k], Acc));
      }
    }
    Eigen::MatrixXcd Ccc_inv = Ccc.inverse();
    CCDeflGuesser.reset(new LuscherGuesser<CoarseCoarseVector>(psi_cc, Ccc_inv));
    ccSolve = &L3PGCR;
    ccGuess = CCDeflGuesser.get();
  }

  //////////////////////////////////////////////////////////////////////
  // Coarse-level GCR smoother for Level 1→2 V-cycle.
  //////////////////////////////////////////////////////////////////////
  RealD coarse_smoother_shift = 0.1;
  int   coarse_smoother_nstep = 2;
  if(getenv("coarse_smoother_shift")) coarse_smoother_shift = atof(getenv("coarse_smoother_shift"));
  if(getenv("coarse_smoother_nstep")) coarse_smoother_nstep = atoi(getenv("coarse_smoother_nstep"));

  ShiftedLinearOperator<CoarseVector> ShiftedLinOpCoarse(coarse_smoother_shift, LinOpCoarse);
  PrecGeneralisedConjugateResidualNonHermitian<CoarseVector> CoarseSmootherGCR(0.01,1,ShiftedLinOpCoarse,simpleC,coarse_smoother_nstep,coarse_smoother_nstep);
  CoarseSmootherGCR.Level(2);
  CoarseSmootherGCR.Name("Csmoother");
  CoarseSmootherGCR.SetZeroGuess(1);  // post-smoother slot: caller zeroes vec2 (NOT L2MGsolver: it takes the Luscher guess)

  //////////////////////////////////////////////////////////////////////
  // Level 1→2 V-cycle preconditioner.
  //////////////////////////////////////////////////////////////////////
  L1to2MG L1to2Precon(AggregatesL2,
                       LinOpCoarse,
                       simpleC,            // no pre-smoother
                       CoarseSmootherGCR,  // post-smoother: shallow shifted GCR
                       LinOpCC,
                       *ccSolve,
                       *ccGuess);

  //////////////////////////////////////////////////////////////////////
  // Standalone Level 1 two-level solve test.
  //////////////////////////////////////////////////////////////////////
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout<<GridLogMessage<<" Level 1 solve: two-level MG preconditioned PGCR"<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;

  PrecGeneralisedConjugateResidualNonHermitian<CoarseVector> L2MGsolver(3.0e-2,200,LinOpCoarse,L1to2Precon,16,16);
  L2MGsolver.Level(2);
  L2MGsolver.Name("Couter");
  c_res = Zero();
  L2MGsolver(c_src,c_res);

  std::cout << GridLogMessage << "Level 1 two-level test: PVdagM operator uses:" << std::endl;
  PVdagM.getApplications();
  PVdagM.clear();

  //////////////////////////////////////////////////////////////////////
  // Full three-level outer solve
  //////////////////////////////////////////////////////////////////////
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout<<GridLogMessage<<" Three-level outer solve (dense CC bottom)"<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;

  PrecGeneralisedConjugateResidualNonHermitian<LatticeFermionD> SmootherGCR(0.01,1,ShiftedPVdagM,simple_fine,8,8);
  SmootherGCR.SetZeroGuess(1);        // pre+post smoother slots both zero their guess: saves 2 fine mults/outer
  SmootherGCR.Level(1);
  SmootherGCR.Name("Fsmoother");

  f_src = one;

  TwoLevelMG ThreeLevelPrecon(AggregatesPD,
                               PVdagM,
                               simple_fine,
                               SmootherGCR,
                               LinOpCoarse,
                               L2MGsolver,
                               simpleC);

  PrecGeneralisedConjugateResidualNonHermitian<LatticeFermion> L1PGCR(1.0e-8,1000,PVdagM,ThreeLevelPrecon,16,16);
  L1PGCR.Level(1);
  L1PGCR.Name("Fouter");

  f_res = Zero();
  L1PGCR(f_src,f_res);

  std::cout << GridLogMessage << "Three-level outer solve: PVdagM operator uses:" << std::endl;
  PVdagM.getApplications();
  PVdagM.clear();
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls   = 24;
  RealD     M5   = 1.8;
  RealD     b    = 1.5;
  RealD     c    = 0.5;
  RealD     mass = 0.00078;
  if ( getenv("MASS") ) mass = atof(getenv("MASS"));

  const int nbasis = 60;

  std::cout << GridLogMessage << "Mass: " << mass << ", Ls: " << Ls << ", b=" << b << ", c=" << c << std::endl;
  std::cout << GridLogMessage << "nbasis: " << nbasis << std::endl;

  std::vector<int> lat_size {48, 48, 48, 96};

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(lat_size, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  // Level 1 coarse grid: block 2^4 from fine (48×48×48×96 → 24×24×24×48, Ls=1)
  Coordinate clatt = lat_size;
  //  Coordinate Block1({2,2,2,2});
  //  Coordinate Block2({8,4,3,6});
  //  Coordinate Block1({2,2,3,3});
  //  Coordinate Block2({8,4,2,4});
  //  Coordinate Block1({4,2,3,3}); // 144s with Fsmoother 12
  //  Coordinate Block2({4,4,2,4});
  Coordinate Block1({4,4,3,3}); 
  Coordinate Block2({4,2,2,4});
  for (int d = 0; d < 4; d++) clatt[d] /= Block1[d];
  std::cout << GridLogMessage << "Level 1 coarse lattice: " << clatt << std::endl;

  GridCartesian *Coarse4d  = SpaceTimeGrid::makeFourDimGrid(clatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridCartesian *Coarse5d  = SpaceTimeGrid::makeFiveDimGrid(1,Coarse4d);

  // Level 2 SUPERCOARSE grid: maximal block {8,4,3,6} from Level 1:
  //   [24,24,24,48] -> [3,6,8,8] = 1152 sites -- the geometric floor on
  //   --mpi 3.6.4.4: x>=3 (mpi 3), y>=6 (mpi 6), z,t>=8 (mpi 4 x SIMD-even).
  //   CC local = [1,1,2,2] = 4 sites/rank: z,t-local = 2 even (SIMD {1,1,2,2}) OK.
  //   Rows/rank for the dense inverse = 4 x 60 = 240, perfectly uniform.
  //   Validated (3level.supercoarse logs): outer count HELD at 34 under this
  //   blocking; the iterative CC solve was the sole remaining cost -- which the
  //   dense inverse removes.
  Coordinate clatt2 = clatt;
  for (int d = 0; d < 4; d++) clatt2[d] /= Block2[d];
  std::cout << GridLogMessage << "Level 2 supercoarse lattice: " << clatt2 << std::endl;

  GridCartesian *CoarseCoarse4d = SpaceTimeGrid::makeFourDimGrid(clatt2, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridCartesian *CoarseCoarse5d = SpaceTimeGrid::makeFiveDimGrid(1,CoarseCoarse4d);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG  RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG  RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeGaugeField Umu(UGrid);
  std::cout << GridLogMessage << "Reading gauge field" << std::endl;
  FieldMetaData header;
  std::string file("/ccs/home/poare/ckpoint_lat.1000");
  NerscIO::readConfiguration(Umu,header,file);

  RealD b_ = 1.5;
  RealD c_ = 0.5;
  MobiusFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b_,c_);
  MobiusFermionD Dpv (Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,1.0, M5,b_,c_);

  typedef PVdagMLinearOperator<MobiusFermionD,LatticeFermionD>        PVdagM_t;
  typedef ShiftedPVdagMLinearOperator<MobiusFermionD,LatticeFermionD> ShiftedPVdagM_t;
  typedef GeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis>  LittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector                           CoarseVector;
  typedef Aggregation<vSpinColourVector,vTComplex,nbasis>             Subspace;
  typedef MGPreconditioner<vSpinColourVector,vTComplex,nbasis>        TwoLevelMG;

  PVdagM_t        PVdagM(Ddwf,Dpv);
  ShiftedPVdagM_t ShiftedPVdagM(0.2,Ddwf,Dpv);

  NextToNearestStencilGeometry5D geom(Coarse5d);

  // Subspace cache: save after generation, reload on subsequent runs to skip expensive setup.
  // Set SUBSPACE_FILE to override the default path.
  std::string subspace_file = "/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/subspace_nb"
                              + std::to_string(nbasis) + ".scidac";
  if ( getenv("SUBSPACE_FILE") ) subspace_file = std::string(getenv("SUBSPACE_FILE"));

  // Check if subspace file exists (boss rank checks, result broadcast via GlobalSum).
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
    // Insurance: GLOBAL (whole-lattice) orthonormalise, matching CreateSubspaceGCR
    // (Aggregates.h:196), in case the cached file predates it.  Span-preserving
    // and globally orthonormal -- NOT the block Orthogonalise() below, which would
    // defeat the raw-null discipline (runMG promotes the RAW subspace to build L2;
    // block-GS here -> psi_coarse = e_k).  The raw copy in runMG happens AFTER this.
    AggregatesGCR.GlobalOrthonormalise();
    //    AggregatesGCR.Orthogonalise();
    std::cout << GridLogMessage << "Subspace loaded, globally orthonormalised (raw block basis preserved)." << std::endl;
  } else {
    std::cout << GridLogMessage << "*** GCR subspace generation ***" << std::endl;
    AggregatesGCR.CreateSubspaceGCR(RNG5,PVdagM,nbasis);
    std::cout << GridLogMessage << "Subspace generation: PVdagM operator uses:" << std::endl;
    PVdagM.getApplications();
    PVdagM.clear();
    saveSubspace(AggregatesGCR.subspace, subspace_file);
    std::cout << GridLogMessage << "Subspace saved to: " << subspace_file << std::endl;
  }

  runMG<nbasis,PVdagM_t,ShiftedPVdagM_t,Subspace,LittleDiracOperator,CoarseVector,TwoLevelMG>(
    FGrid,
    Coarse5d,
    CoarseCoarse5d,
    geom,
    PVdagM,
    ShiftedPVdagM,
    AggregatesGCR
  );

  std::cout << GridLogMessage << "Done" << std::endl;
  Grid_finalize();
  return 0;
}
