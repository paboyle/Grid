/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: Test_schur_dense_coarse.cc

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
// T6 of the RecursiveSchurInverse regression chain
// (schur_recursive_inverse_plan.txt 4B.5): the DenseCoarseMatrix GLUE,
// on a real (tiny) lattice coarse operator, CPU laptop build.
//
// Builds a genuine GeneralCoarsenedMatrix (DWF MdagM + 0.5 shift for a
// guaranteed-invertible Galerkin coarse op, random aggregation basis,
// nbasis=8, 4^4 x Ls/1 blocking) and constructs DenseCoarseMatrix in
// DENSE_SCHUR=2 AUDIT mode with small DENSE_PANEL_BYTES (multi-panel
// gathers exercised through the glue).  The constructor then runs, in
// order, all the certificates this stage exists to check:
//   - fresh ImportDense (no SLAB_FILE) + IMPORT CERTIFICATE vs Op.M
//   - InvertDenseSingle (the oracle)
//   - InvertDenseSchur: self-certifying rank-major map, fp64 diagonal
//     import certificate vs the fp32 slab, distributed recursion,
//     growth telemetry
//   - AUDIT: max|Ainv_schur - Ainv_single| over the full slab
//   - VERIFY ||A Ainv x - x||/||x|| through the SCHUR result
// This program adds asserts on the audit number and a random-vector
// round trip.
//
// Uniform local volume 12.12.12.12 (fine), per-dim blocks {4,4,3,3},
// coarse 3.3.4.4/rank, nbasis 4 (N = 576n):
//   mpirun -n 1 ./Test_schur_dense_coarse --grid 12.12.12.12 --mpi 1.1.1.1
//   mpirun -n 2 ./Test_schur_dense_coarse --grid 12.12.12.24 --mpi 1.1.1.2
//   mpirun -n 3 ./Test_schur_dense_coarse --grid 12.12.12.36 --mpi 1.1.1.3
//   mpirun -n 4 ./Test_schur_dense_coarse --grid 12.12.12.48 --mpi 1.1.1.4
//
#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>
#include <Grid/algorithms/multigrid/DenseCoarseMatrix.h>

using namespace std;
using namespace Grid;

///////////////////////////////////////////////////////////////////////
// MdagM + shift: Galerkin projection of a PD operator plus sigma I is
// safely invertible whatever the (random) subspace quality.
///////////////////////////////////////////////////////////////////////
template<class Field>
class ShiftedHermOpAdaptor : public LinearOperatorBase<Field>
{
  LinearOperatorBase<Field> &wrapped;
  RealD                      shift;
public:
  ShiftedHermOpAdaptor(LinearOperatorBase<Field> &wrapme, RealD s)
    : wrapped(wrapme), shift(s) {};

  void Op(const Field &in, Field &out)
  {
    wrapped.HermOp(in, out);
    out = out + shift*in;
  }
  void AdjOp(const Field &in, Field &out)
  {
    Op(in, out);
  }
  void HermOp(const Field &in, Field &out)
  {
    Op(in, out);
  }
  void OpDiag(const Field &in, Field &out)                        { GRID_ASSERT(0); }
  void OpDir (const Field &in, Field &out, int dir, int disp)     { GRID_ASSERT(0); }
  void OpDirAll(const Field &in, std::vector<Field> &out)         { GRID_ASSERT(0); }
  void HermOpAndNorm(const Field &in, Field &out, RealD &n1, RealD &n2) { GRID_ASSERT(0); }
};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls     = 4;
  const int nbasis = 4;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                                   GridDefaultSimd(Nd,vComplex::Nsimd()),
                                                                   GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  // Per-dimension blocking {4,4,3,3}: fine 12.12.12.12 -> coarse
  // 3.3.4.4.  Two constraints meet here (both MEASURED today):
  //  - coarse dims of 2 hit the probing pathology (health probe below)
  //  - GEN-simd lanes {1,1,2,2} must land on even coarse dims, so the
  //    odd production-like 3s go on the lane-free x,y axes (exactly the
  //    production [3,6,8,8] trick).
  Coordinate blocks({4,4,3,3});
  Coordinate clatt = GridDefaultLatt();
  for(int d=0; d<clatt.size(); d++)
  {
    GRID_ASSERT( (clatt[d] % blocks[d]) == 0 );
    clatt[d] = clatt[d]/blocks[d];
  }
  GridCartesian *Coarse4d = SpaceTimeGrid::makeFourDimGrid(clatt,
                                                           GridDefaultSimd(Nd,vComplex::Nsimd()),
                                                           GridDefaultMpi());
  GridCartesian *Coarse5d = SpaceTimeGrid::makeFiveDimGrid(1,Coarse4d);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  std::vector<int> cseeds({9,10,11,12});
  GridParallelRNG RNG4(UGrid);    RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG RNG5(FGrid);    RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG CRNG(Coarse5d); CRNG.SeedFixedIntegers(cseeds);

  LatticeGaugeField Umu(UGrid);
  SU<Nc>::HotConfiguration(RNG4,Umu);

  RealD mass = 0.1;
  RealD M5   = 1.8;
  DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);

  MdagMLinearOperator<DomainWallFermionD,LatticeFermion> HermDefOp(Ddwf);
  ShiftedHermOpAdaptor<LatticeFermionD> HOA(HermDefOp, 0.5);

  std::cout << GridLogMessage << "Building random aggregation space, nbasis " << nbasis << std::endl;
  typedef Aggregation<vSpinColourVector,vTComplex,nbasis> Subspace;
  Subspace Aggregates(Coarse5d,FGrid,0);
  Aggregates.CreateSubspaceRandom(RNG5);

  std::cout << GridLogMessage << "Coarsening shifted MdagM" << std::endl;
  typedef GeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> LittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector CoarseVector;
  NextToNextToNextToNearestStencilGeometry5D geom(Coarse5d);
  LittleDiracOperator LittleDiracOp(geom,FGrid,Coarse5d);
  LittleDiracOp.CoarsenOperator(HOA,Aggregates);

  ///////////////////////////////////////////////////////////////////////
  // Operator health probes (independent of DenseCoarseMatrix).
  //
  // MEASURED PATHOLOGY, banked 2026-08-14: on coarse dims of 2 (fine
  // 8.8.8.8, block 4 -> coarse 2.2.2.2) the coarsened operator is
  // rank 16/128 with 112 zero ROWS (output support = 2 of 16 sites)
  // and Hermiticity violation 0.17 -- the probing construction breaks
  // on the size-2 torus.  The import certificate cannot see this
  // (dense and M share _A).  Out of scope here; coarse dims >= 3.
  //
  // Cheap any-size probes: output support + Hermiticity via inner
  // products on random vectors.
  ///////////////////////////////////////////////////////////////////////
  {
    CoarseVector px(Coarse5d);
    CoarseVector py(Coarse5d);
    CoarseVector Mx(Coarse5d);
    CoarseVector My(Coarse5d);
    random(CRNG, px);
    random(CRNG, py);
    LittleDiracOp.M(px, Mx);
    LittleDiracOp.M(py, My);
    ComplexD ip1 = innerProduct(py, Mx);        // <y, M x>
    ComplexD ip3 = innerProduct(px, My);        // <x, M y>
    RealD hermdev = abs(ip1 - conj(ip3)) / std::sqrt(norm2(Mx)*norm2(py));
    RealD support = norm2(Mx) / norm2(px);
    std::cout << GridLogMessage << "Operator health: ||Mx||^2/||x||^2 = " << support
              << "  herm-dev " << hermdev << std::endl;
    // Hermitian fine op => exactly Hermitian Galerkin coarse op.
    // (A measured herm-dev of 1.3e-4 here was the CPU SIMT-lane
    // CoarsenOperator bug -- fixed 2026-08-14, now 4e-15.  A loud
    // failure here means _A population is broken again.)
    GRID_ASSERT( support > 1.0e-3 );
    GRID_ASSERT( hermdev < 1.0e-10 );
  }

  ///////////////////////////////////////////////////////////////////////
  // Full-matrix conditioning probe at small N: dense columns by
  // applying M to unit vectors, fp64 Eigen SVD.
  ///////////////////////////////////////////////////////////////////////
  {
    int64_t Nprobe = Coarse5d->gSites() * nbasis;
    if ( Nprobe <= 700 )
    {
      Eigen::MatrixXcd eA(Nprobe, Nprobe);
      CoarseVector e(Coarse5d);
      CoarseVector Me(Coarse5d);
      for(int64_t j=0; j<Nprobe; j++)
      {
        int64_t gsite = j / nbasis;
        int     b     = j % nbasis;
        e = Zero();
        Coordinate gcoor(Coarse5d->_ndimension);
        Lexicographic::CoorFromIndex(gcoor, gsite, Coarse5d->GlobalDimensions());
        typedef typename CoarseVector::vector_object::scalar_object csobj;
        csobj s;
        s = Zero();
        ((ComplexD *)&s)[b] = ComplexD(1.0,0.0);
        pokeSite(s, e, gcoor);
        LittleDiracOp.M(e, Me);
        for(int64_t i=0; i<Nprobe; i++)
        {
          int64_t gsi = i / nbasis;
          int     bi  = i % nbasis;
          Coordinate gci(Coarse5d->_ndimension);
          Lexicographic::CoorFromIndex(gci, gsi, Coarse5d->GlobalDimensions());
          csobj si;
          peekSite(si, Me, gci);
          // Explicit re/im at the thrust/std boundary (HIP builds)
          ComplexD zz = ((ComplexD *)&si)[bi];
          eA(i,j) = std::complex<double>(zz.real(), zz.imag());
        }
      }
      Eigen::JacobiSVD<Eigen::MatrixXcd> svd(eA);
      double smax = svd.singularValues()(0);
      double smin = svd.singularValues()(Nprobe-1);
      int64_t rank = 0;
      for(int64_t i=0; i<Nprobe; i++)
      {
        if ( svd.singularValues()(i) > 1.0e-10*smax ) rank++;
      }
      double herm = (eA - eA.adjoint()).cwiseAbs().maxCoeff();
      int64_t zrows = 0;
      int64_t zcols = 0;
      for(int64_t i=0; i<Nprobe; i++)
      {
        if ( eA.row(i).cwiseAbs().maxCoeff() < 1.0e-12 ) zrows++;
        if ( eA.col(i).cwiseAbs().maxCoeff() < 1.0e-12 ) zcols++;
      }
      std::cout << GridLogMessage << "Operator probe: N=" << Nprobe
                << "  sigma_max " << smax
                << "  sigma_min " << smin
                << "  rank " << rank << "/" << Nprobe
                << "  herm-dev " << herm
                << "  zero rows/cols " << zrows << "/" << zcols
                << std::endl;
    }
  }

  ///////////////////////////////////////////////////////////////////////
  // T6: AUDIT mode, fresh import, multi-panel gathers.  The constructor
  // runs every certificate in the chain (see banner).
  ///////////////////////////////////////////////////////////////////////
  setenv("DENSE_SCHUR","2",1);
  setenv("DENSE_PANEL_BYTES","65536",1);
  unsetenv("SLAB_FILE");

  // Restructured interface: DenseCoarseMatrix<CComplex,nbasis>, constructed
  // on the grid and fed by Import (which runs the certificate chain).
  typedef DenseCoarseMatrix<vTComplex,nbasis> DenseCC;
  DenseCC dcm(Coarse5d);
  dcm.Import(LittleDiracOp);

  std::cout << GridLogMessage << "T6  audit relative slab difference (schur vs single) = "
            << dcm.schurAuditRel << std::endl;
  GRID_ASSERT( dcm.schurAuditRel >= 0.0 );      // audit actually ran
  GRID_ASSERT( dcm.schurAuditRel <  1.0e-3 );

  ///////////////////////////////////////////////////////////////////////
  // Random-vector round trip through the SCHUR inverse
  ///////////////////////////////////////////////////////////////////////
  CoarseVector x(Coarse5d);
  CoarseVector y(Coarse5d);
  CoarseVector z(Coarse5d);
  random(CRNG, x);
  dcm(x, y);
  LittleDiracOp.M(y, z);
  z = z - x;
  RealD rel = std::sqrt(norm2(z)/norm2(x));
  std::cout << GridLogMessage << "T6  round trip ||A Ainv x - x||/||x|| (random x) = "
            << rel << std::endl;
  GRID_ASSERT( rel < 1.0e-2 );

  std::cout << GridLogMessage << "Test_schur_dense_coarse: T6 ALL PASS" << std::endl;

  Grid_finalize();
}
