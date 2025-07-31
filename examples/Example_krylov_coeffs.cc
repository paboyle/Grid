/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_padded_cell.cc

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

// Tests code written to read off the Krylov coefficients

#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>

#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidual.h>
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidualNonHermitian.h>
#include <Grid/algorithms/iterative/BiCGSTAB.h>
#include <Grid/algorithms/iterative/ConjugateGradient.h>

using namespace std;
using namespace Grid;

// Hermitize a DWF operator by squaring it
template<class Matrix,class Field>
class SquaredLinearOperator : public LinearOperatorBase<Field> {

  public:
  Matrix &_Mat;

  public:
    SquaredLinearOperator(Matrix &Mat): _Mat(Mat) {};

    void OpDiag (const Field &in, Field &out) {    assert(0);  }
    void OpDir  (const Field &in, Field &out,int dir,int disp) {    assert(0);  }
    void OpDirAll  (const Field &in, std::vector<Field> &out){    assert(0);  };
    void Op     (const Field &in, Field &out){
      // std::cout << "Op is overloaded as HermOp" << std::endl;
      HermOp(in, out);
    }
    void AdjOp     (const Field &in, Field &out){
      HermOp(in, out);
    }
    void _Op     (const Field &in, Field &out){
      // std::cout << "Op: M "<<std::endl;
      _Mat.M(in, out);
    }
    void _AdjOp     (const Field &in, Field &out){
      // std::cout << "AdjOp: Mdag "<<std::endl;
      _Mat.Mdag(in, out);
    }
    void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){    assert(0);  }
    void HermOp(const Field &in, Field &out){
      // std::cout << "HermOp: Mdag M Mdag M"<<std::endl;
      Field tmp(in.Grid());
      _Op(in,tmp);
      _AdjOp(tmp,out);
    }
};

/**
 * Computes the coefficients in the Krylov expansion for 1/D ~ \sum_{i=0}^N c_i D^i. 
 * 
 * Parameters
 * ----------
 * std::vector<double> &coeffs
 *    Polynomial coeffients to return, with indexing order (c_0, c_1, c_2, ..., c_n). 
 * LinearOperatorBase<FineField> &DiracOp
 *    Dirac operator D. 
 * FineField src
 *    Source field b. 
 * FineField psiStar
 *    Output approximation for D^{-1} b coming from a Krylov method. 
 * int N
 *    Dimension of the polynomial approximation (Krylov space K_{N-1} = {b, Db, D^2 b, ..., D^{N-1} b}).
 */
void poly_coeffs(std::vector<std::complex<double>> &coeffs, LinearOperatorBase<LatticeFermion> &DiracOp, LatticeFermion src, 
                  LatticeFermion psiStar, GridCartesian* FGrid, int N, bool use_herm = false)
{
  // stdBasis = {b, Db, D^2 b, ..., D^N b}, kryBasis = {k0, k1, ..., kN}
  std::vector<LatticeFermion> kryBasis;
  Eigen::VectorXcd psiStarCoeffs (N);

  // Normalize by 1 / ||src||; does not change the polynomial coefficients
  double srcNorm   = 1 / std::sqrt(norm2(src));
  kryBasis.push_back(srcNorm * src);                // normalized source
  psiStar          = srcNorm * psiStar;
  psiStarCoeffs(0) = innerProduct(kryBasis[0], psiStar);

  // orthonormalize canonical Krylov basis {b, Db, D^2 b, ..., D^{N-1} b} <--> {k_i} and compute components <k_i | psi*>
  LatticeFermion tmp (FGrid);
  for (int i = 0; i < N - 1; i++) {               // construct ONB for {b, Db, ..., D^{i+1} b}
    if (use_herm) {
      DiracOp.HermOp(kryBasis.back(), tmp);         // tmp \in span{(D^\dag D)^{i+1} b} \oplus span{(D^\dag D)^i b, ..., D^\dag D b, b}
    } else {
      DiracOp.Op(kryBasis.back(), tmp);             // tmp \in span{D^{i+1} b} \oplus span{D^i b, ..., Db, b}
    }

    for (int j = 0; j < i+1; j++) {               // orthogonalize tmp with previous basis vectors
      ComplexD coeff = innerProduct(kryBasis[j], tmp);      // <k_j | tmp>
      tmp -= coeff * kryBasis[j];                           // subtract off |k_j><k_j | tmp>; now tmp is perp to |k_j>
    }
    double tmpNorm = 1 / std::sqrt(norm2(tmp));
    kryBasis.push_back(
      tmpNorm * tmp
    );                                                      // normalize |k_i> and add to kryBasis
    psiStarCoeffs(i+1) = innerProduct(kryBasis[i+1], psiStar);  // compute < k_i | psi* >
  }

  // To verify the basis is ONB
  // for (int i = 0; i < N; i++) {
  //   for (int j = 0; j < N; j++) {
  //     std::cout << "<ki|kj> for (i,j) = (" << i << ", " << j << ") = "  << innerProduct(kryBasis[i], kryBasis[j]) << std::endl;
  //   }
  // }

  // Compute change of basis matrix
  LatticeFermion tmp2 (FGrid);
  Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(N, N);
  tmp = kryBasis[0];       // current Krylov vector; starts with tmp = src (normalized)
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < i + 1; j++) {    // fill column with components of kryVec. Only need j <= i to get orthonormal components
      M(j, i) = innerProduct(kryBasis[j], tmp);
    }    
    if (use_herm) {     // tmp --> D^\dag D(tmp)
      DiracOp.HermOp(tmp, tmp2);
      tmp = tmp2;
    } else {      // tmp --> D(tmp). Note that DiracOp.Op(tmp, tmp) will cause a bug
      DiracOp.Op(tmp, tmp2);
      tmp = tmp2;
    }
  }
  
  // Compute M^{-1} @ psiStarCoeffs and copy to coeffs
  Eigen::VectorXcd res (N);
  res = M.inverse() * psiStarCoeffs;
  for (int i = 0; i < N; i++) {
    coeffs[i] = res(i);
  }

}

// out file for poly coefficients (should it be complex?)
// class PolynomialFile: Serializable {
// public:
//   GRID_SERIALIZABLE_CLASS_MEMBERS(OutputFile, std::vector< Real >, data);
// };

std::complex<double> poly_approx(std::complex<double> x, std::vector<std::complex<double>> coeffs) {
  std::complex<double> px;
  for (int i = 0; i < coeffs.size(); i++) {
    px += coeffs[i] * std::pow(x, i);
  }
  return px;
}

/**
 * Returns the approximation psi = \sum_i c_i D^i b resulting from a Krylov solver.
 * 
 * Parameters
 * ----------
 * LatticeFermion &psi
 *    Approximation field, returned psi = \sum_i c_i D^i b.
 * LatticeFermion src
 *    Source b used to generate the Krylov space K_n(D, b).
 * LinearOperatorBase<LatticeFermion> &Linop
 *    Dirac operator used to generate the Krylov space K_n(D, b).
 * std::vector<std::complex<double>> coeffs
 *    Polynomial coefficients returned from the solver. 
 */
void krylovApprox(LatticeFermion &psi, LatticeFermion src, LinearOperatorBase<LatticeFermion> &Linop, std::vector<std::complex<double>> coeffs) {
  psi = Zero();
  LatticeFermion tmp (psi.Grid());
  tmp = src;
  LatticeFermion tmp2 (psi.Grid());
  for (int i = 0; i < coeffs.size(); i++) {
      psi = psi + coeffs[i] * tmp;
      Linop.Op(tmp, tmp2);              // tmp = D*tmp
      tmp = tmp2;
  }
}

int main (int argc, char ** argv)
{
  Grid_init(&argc, &argv);
  const int Ls = 8;
  std::vector<int> lat_size {16, 16, 16, 32};
  std::cout << "Lattice size: " << lat_size << std::endl;
  GridCartesian * UGrid = SpaceTimeGrid::makeFourDimGrid(lat_size, 
								          GridDefaultSimd(Nd,vComplex::Nsimd()),
								          GridDefaultMpi());

  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  //////////////////////////////////////////////////////////////////////
  // You can manage seeds however you like.
  // Recommend SeedUniqueString.
  //////////////////////////////////////////////////////////////////////
  // std::vector<int> seeds4({1, 2, 3, 4}); 
  // GridParallelRNG RNG4(UGrid);
  // RNG4.SeedFixedIntegers(seeds4);

  // std::vector<int> seeds5({1, 2, 3, 4, 5}); 
  // GridParallelRNG RNG5(FGrid);
  // RNG5.SeedFixedIntegers(seeds5);

  // std::string outStrStem = "/Users/patrickoare/Dropbox (MIT)/research/multigrid/grid_out/";

  LatticeGaugeField Umu(UGrid);
  FieldMetaData header;
  std::string file("/Users/patrickoare/libraries/PETSc-Grid/ckpoint_lat.4000");
  NerscIO::readConfiguration(Umu, header, file);
  
  RealD mass=0.01;
  RealD M5=1.8;
  // RealD M5=1.0;
  RealD b=1.5;// Scale factor b+c=2, b-c=1
  RealD c=0.5;

  // load in Dirac operators that we'll use; square it to Hermitize
  // Dsq just needs to be a Hermitian operator so we can use CG on it
  DomainWallFermionD Ddwf(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5);
  SquaredLinearOperator<DomainWallFermionD, LatticeFermionD> Dsq (Ddwf);
  NonHermitianLinearOperator<DomainWallFermionD, LatticeFermionD> DLinOp (Ddwf);

  LatticeFermion src (FGrid); src = 1.0;                              // Source to use 
  LatticeFermion psiCG (FGrid); psiCG = Zero();                       // Field to solve with for CG
  LatticeFermion psiGCR (FGrid); psiGCR = Zero();                     // Field to solve with for GCR

  std::cout << GridLogMessage << "*******************************************" << std::endl;
  std::cout << GridLogMessage << "********** TESTING CG POLY COEFFS *********" << std::endl;
  std::cout << GridLogMessage << "*******************************************" << std::endl << std::endl;

  double tol = 1.0e-8;
  int N = 5;           // max iterations (size of Krylov basis)

  // GCR variables
  int outer_iters = 1;                  // num restarts for GCR
  TrivialPrecon<LatticeFermionD> prec;  // trivial preconditioner

  ConjugateGradientPolynomial<LatticeFermion> CGP(tol, N, false);
  CGP(Dsq, src, psiCG);

  // Compute Krylov coeffs directly and compare
  std::vector<std::complex<double>> cg_coeffs (N);
  poly_coeffs(cg_coeffs, Dsq, src, psiCG, FGrid, N, true);

  PolynomialFile PF;

  // Use GCR solver, also get poly coeffs
  std::vector<std::complex<double>> gcr_sym_coeffs (N);     // Can try N --> N + 3 to test to see if the last 3 comps are 0
  PGCRPolynomial<LatticeFermionD> GCRPolySym(tol, outer_iters, Dsq, prec, N+1, N, PF);    // mmax sets the memory, note the last beta doesn't really matter for updating the polynomial
  GCRPolySym(src, psiGCR);
  // poly_coeffs(gcr_sym_coeffs, Dsq, src, psi, FGrid, N, true);
  poly_coeffs(gcr_sym_coeffs, Dsq, src, psiGCR, FGrid, N, true);

  std::cout << GridLogMessage << std::endl << "******** CG POLYNOMIAL COEFFICIENTS *******" << std::endl;
  std::cout << GridLogMessage << CGP.polynomial << std::endl << std::endl;

  std::cout << GridLogMessage << "****** DIRECT POLYNOMIAL COEFFICIENTS *****" << std::endl;
  std::cout << GridLogMessage << cg_coeffs << std::endl << std::endl;

  // TODO: try GCR with a Hermitian operator (Dsq)
  std::cout << GridLogMessage << "****** GCR COEFFICIENTS *****" << std::endl;
  std::cout << GridLogMessage << GCRPolySym.polynomial << std::endl << std::endl;

  std::cout << GridLogMessage << "****** DIRECT GCR COEFFICIENTS *****" << std::endl;
  std::cout << GridLogMessage << gcr_sym_coeffs << std::endl << std::endl;

  // test how good the decomposition is
  std::cout << "Testing fidelity of decomposition by computing ||psi* - sum_i c_i D^i b||^2!" << std::endl;
  LatticeFermion psiPrime (FGrid);

  // for CG
  krylovApprox(psiPrime, src, Dsq, cg_coeffs);
  std::cout << "CG with Dsq, ||psi - psiPrime||^2 = " << norm2(psiCG - psiPrime) << std::endl;

  // for GCR with alpha / beta computation
  krylovApprox(psiPrime, src, Dsq, GCRPolySym.polynomial);
  std::cout << "GCR with Dsq, ||psi - psiPrime||^2 = " << norm2(psiGCR - psiPrime) << std::endl;

  // for GCR with alpha / beta computation
  krylovApprox(psiPrime, src, Dsq, gcr_sym_coeffs);
  std::cout << "GCR direct with Dsq, ||psi - psiPrime||^2 = " << norm2(psiGCR - psiPrime) << std::endl;

  
  // std::vector<double> real_cg_diff (N);
  // for (int i = 0; i < N; i++) { real_cg_diff[i] = std::abs(cg_coeffs[i].real() - CGP.polynomial[i]); }
  // std::cout << GridLogMessage << "************* COEFF DIFFERENCE ************" << std::endl;
  // std::cout << GridLogMessage << real_cg_diff << std::endl << std::endl;

  // GCR polynomial reconstruction with Ddwf!
  std::cout << GridLogMessage << "*******************************************" << std::endl;
  std::cout << GridLogMessage << "********* TESTING GCR POLY COEFFS *********" << std::endl;
  std::cout << GridLogMessage << "*******************************************" << std::endl << std::endl;

  // re-init variables
  src = 1.0;
  src = (1 / std::sqrt(norm2(src))) * src;
  psiGCR = Zero(); psiPrime = Zero();

  // test GCR poly
  PGCRPolynomial<LatticeFermionD> GCRPoly(tol, outer_iters, DLinOp, prec, N+1, N, PF);    // mmax sets the memory, note the last beta doesn't really matter for updating the polynomial
  GCRPoly(src, psiGCR);

  // Compute Krylov coeffs directly and compare
  // N = 1;    // compare the N > 1 decomposition with the psi* resulting from N = 1
  std::vector<std::complex<double>> gcr_coeffs (N);   // note N --> N + k should just give k coeffs that are 0; this works as intended
  poly_coeffs(gcr_coeffs, DLinOp, src, psiGCR, FGrid, N, false);

  std::cout << GridLogMessage << "******* GCR POLYNOMIAL COEFFICIENTS *******" << std::endl;
  std::cout << GridLogMessage << GCRPoly.polynomial << std::endl << std::endl;

  std::cout << GridLogMessage << "****** DIRECT POLYNOMIAL COEFFICIENTS *****" << std::endl;
  std::cout << GridLogMessage << gcr_coeffs << std::endl << std::endl;

  // test how good the decomposition is
  std::cout << "Testing fidelity of decomposition by computing ||psi* - sum_i c_i D^i b||^2!" << std::endl;

  // for GCR with alpha / beta computation
  krylovApprox(psiPrime, src, DLinOp, GCRPoly.polynomial);
  std::cout << "GCR with Dsq, ||psi - psiPrime||^2 = " << norm2(psiGCR - psiPrime) << std::endl;

  // for GCR with alpha / beta computation
  krylovApprox(psiPrime, src, DLinOp, gcr_coeffs);
  std::cout << "GCR direct with Dsq, ||psi - psiPrime||^2 = " << norm2(psiGCR - psiPrime) << std::endl;

  // TESTS TO DO THE N = 2 CASE DIRECTLY
  /*
  std::vector<std::complex<double>> alphas {
    std::complex(0.244300601, 0.00013007545), 
    std::complex(0.285370971, -0.000160704481)
  };
  std::complex<double> beta00 (-0.184661284, -6.52153945e-05);
  LatticeFermion psi2 (FGrid);
  LatticeFermion Dsrc (FGrid);
  DLinOp.Op(src, Dsrc);
  std::complex<double> c1 = alphas[0] + alphas[1] * (1. + beta00);
  std::complex<double> c2 = -alphas[0] * alphas[1];
  psi2 = c1 * src + c2 * Dsrc;

  std::cout << "||b|| = " << norm2(src) << std::endl;
  std::cout << "||Db|| = " << norm2(Dsrc) << std::endl;
  // fail; so far this is giving something different than what's being computed in krylovApprox (idk how?)

  std::cout << "c1 and c2 are: " << c1 << " and " << c2 << std::endl;
  std::cout << "GCRPoly polynomial coeffs are (should equal c1 and c2): " << GCRPoly.polynomial << std::endl;
  std::cout << "||GCRpsi - psi2||_2^2 = " << norm2(psiGCR - psi2) << std::endl;
  // pass
  
  LatticeFermion src2 (FGrid);
  src2 = 1.0;
  src2 = (1 / std::sqrt(norm2(src2))) * src2;
  std::cout << "||ones - src|| (to verify that src is the same throughout, should be 0) = " << norm2(src2 - src) << std::endl;
  // pass

  krylovApprox(psiPrime, src, DLinOp, GCRPoly.polynomial);
  std::cout << "GCR with Dsq, ||psi2 - psiPrime||^2 = " << norm2(psi2 - psiPrime) << std::endl;

  std::vector<std::complex<double>> psi2_coeffs (N);   // note N --> N + k should just give k coeffs that are 0; this works as intended
  poly_coeffs(psi2_coeffs, DLinOp, src, psi2, FGrid, N, false);
  krylovApprox(psiPrime, src, DLinOp, psi2_coeffs);
  std::cout << "GCR direct with Dsq, ||psi - psiPrime||^2 = " << norm2(psi2 - psiPrime) << std::endl;
  */

  // std::complex z (10.0, 0.0);     // z = 10
  // std::cout << GridLogMessage << "************* GCR POLY(z = 10) *************" << std::endl;
  // std::cout << GridLogMessage << poly_approx(z, GCRPoly.polynomial) << std::endl;
  // std::cout << GridLogMessage << "************ DIRECT POLY(z = 10) ***********" << std::endl;
  // std::cout << GridLogMessage << poly_approx(z, gcr_coeffs) << std::endl;

  // std::vector<std::complex<double>> gcr_diff (N);
  // for (int i = 0; i < N; i++) { gcr_diff[i] = gcr_coeffs[i] - GCRPoly.polynomial[i]; }
  // std::cout << GridLogMessage << "*********** GCR COEFF DIFFERENCE **********" << std::endl;
  // std::cout << GridLogMessage << gcr_diff << std::endl << std::endl;

  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage << "Done "<< std::endl;

  Grid_finalize();
  return 0;
}
