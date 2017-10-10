/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/debug/Test_reweight_dwf_eofa.cc

Copyright (C) 2017

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: David Murphy <dmurphy@phys.columbia.edu>

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
using namespace Grid::QCD;

typedef typename GparityDomainWallFermionR::FermionField FermionField;

// parameters for test
const std::vector<int> grid_dim = { 8, 8, 8, 8 };
const int Ls = 8;
const int Nhits = 10;
const int max_iter = 5000;
const RealD b  = 2.5;
const RealD c  = 1.5;
const RealD mf = 0.1;
const RealD mb = 0.11;
const RealD M5 = 1.8;
const RealD stop_tol = 1.0e-12;

RealD mean(const std::vector<RealD>& data)
{
    int N = data.size();
    RealD mean(0.0);
    for(int i=0; i<N; ++i){ mean += data[i]; }
    return mean/RealD(N);
}

RealD jack_mean(const std::vector<RealD>& data, int sample)
{
    int N = data.size();
    RealD mean(0.0);
    for(int i=0; i<N; ++i){ if(i != sample){ mean += data[i]; } }
    return mean/RealD(N-1);
}

RealD jack_std(const std::vector<RealD>& jacks, RealD mean)
{
    int N = jacks.size();
    RealD std(0.0);
    for(int i=0; i<N; ++i){ std += std::pow(jacks[i]-mean, 2.0); }
    return std::sqrt(RealD(N-1)/RealD(N)*std);
}

std::vector<RealD> jack_stats(const std::vector<RealD>& data)
{
    int N = data.size();
    std::vector<RealD> jack_samples(N);
    std::vector<RealD> jack_stats(2);

    jack_stats[0] = mean(data);
    for(int i=0; i<N; i++){ jack_samples[i] = jack_mean(data,i); }
    jack_stats[1] = jack_std(jack_samples, jack_stats[0]);
    return jack_stats;
}

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  // Initialize spacetime grid
  std::cout << GridLogMessage << "Lattice dimensions: "
    << grid_dim << "   Ls: " << Ls << std::endl;
  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(grid_dim,
      GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian* FGrid = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);
  GridRedBlackCartesian* FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid);

  // Set up RNGs
  std::vector<int> seeds4({1, 2, 3, 4});
  std::vector<int> seeds5({5, 6, 7, 8});
  GridParallelRNG RNG5(FGrid);
  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG RNG4(UGrid);
  RNG4.SeedFixedIntegers(seeds4);

  // Random gauge field
  LatticeGaugeField Umu(UGrid);
  SU3::HotConfiguration(RNG4, Umu);

  // Initialize RHMC fermion operators
  GparityDomainWallFermionR::ImplParams params;
  GparityMobiusFermionR Ddwf_f(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mf, M5, b, c, params);
  GparityMobiusFermionR Ddwf_b(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mb, M5, b, c, params);
  SchurDiagMooeeOperator<GparityMobiusFermionR, FermionField> MdagM(Ddwf_f);
  SchurDiagMooeeOperator<GparityMobiusFermionR, FermionField> VdagV(Ddwf_b);

  // Degree 12 rational approximations to x^(1/4) and x^(-1/4)
  double     lo = 0.0001;
  double     hi = 95.0;
  int precision = 64;
  int    degree = 12;
  AlgRemez remez(lo, hi, precision);
  std::cout << GridLogMessage << "Generating degree " << degree << " for x^(1/4)" << std::endl;
  remez.generateApprox(degree, 1, 4);
  MultiShiftFunction PowerQuarter(remez, stop_tol, false);
  MultiShiftFunction PowerNegQuarter(remez, stop_tol, true);

  // Stochastically estimate reweighting factor via RHMC
  RealD scale = std::sqrt(0.5);
  std::vector<RealD> rw_rhmc(Nhits);
  ConjugateGradientMultiShift<FermionField> msCG_V(max_iter, PowerQuarter);
  ConjugateGradientMultiShift<FermionField> msCG_M(max_iter, PowerNegQuarter);
  std::cout.precision(12);

  for(int hit=0; hit<Nhits; hit++){

    // Gaussian source
    FermionField Phi    (Ddwf_f.FermionGrid());
    FermionField PhiOdd (Ddwf_f.FermionRedBlackGrid());
    std::vector<FermionField> tmp(2, Ddwf_f.FermionRedBlackGrid());
    gaussian(RNG5, Phi);
    Phi = Phi*scale;

    pickCheckerboard(Odd, PhiOdd, Phi);

    // evaluate -log(rw)
    msCG_V(VdagV, PhiOdd, tmp[0]);
    msCG_M(MdagM, tmp[0], tmp[1]);
    rw_rhmc[hit] = norm2(tmp[1]) - norm2(PhiOdd);
    std::cout << std::endl << "==================================================" << std::endl;
    std::cout << " --- RHMC: Hit " << hit << ": rw = " << rw_rhmc[hit];
    std::cout << std::endl << "==================================================" << std::endl << std::endl;

  }

  // Initialize EOFA fermion operators
  RealD shift_L = 0.0;
  RealD shift_R = -1.0;
  int pm = 1;
  GparityMobiusEOFAFermionR Deofa_L(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mf, mf, mb, shift_L, pm, M5, b, c, params);
  GparityMobiusEOFAFermionR Deofa_R(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mb, mf, mb, shift_R, pm, M5, b, c, params);
  MdagMLinearOperator<GparityMobiusEOFAFermionR, FermionField> LdagL(Deofa_L);
  MdagMLinearOperator<GparityMobiusEOFAFermionR, FermionField> RdagR(Deofa_R);

  // Stochastically estimate reweighting factor via EOFA
  RealD k = Deofa_L.k;
  std::vector<RealD> rw_eofa(Nhits);
  ConjugateGradient<FermionField> CG(stop_tol, max_iter);
  SchurRedBlackDiagMooeeSolve<FermionField> SchurSolver(CG);

  // Compute -log(Z), where: ( RHMC det ratio ) = Z * ( EOFA det ratio )
  RealD Z = std::pow(b+c+1.0,Ls) + mf*std::pow(b+c-1.0,Ls);
  Z /= std::pow(b+c+1.0,Ls) + mb*std::pow(b+c-1.0,Ls);
  Z = -12.0*grid_dim[0]*grid_dim[1]*grid_dim[2]*grid_dim[3]*std::log(Z);

  for(int hit=0; hit<Nhits; hit++){

    // Gaussian source
    FermionField Phi       (Deofa_L.FermionGrid());
    FermionField spProj_Phi(Deofa_L.FermionGrid());
    std::vector<FermionField> tmp(2, Deofa_L.FermionGrid());
    gaussian(RNG5, Phi);
    Phi = Phi*scale;

    // evaluate -log(rw)
    // LH term
    for(int s=0; s<Ls; ++s){ axpby_ssp_pminus(spProj_Phi, 0.0, Phi, 1.0, Phi, s, s); }
    Deofa_L.Omega(spProj_Phi, tmp[0], -1, 0);
    G5R5(tmp[1], tmp[0]);
    tmp[0] = zero;
    SchurSolver(Deofa_L, tmp[1], tmp[0]);
    Deofa_L.Dtilde(tmp[0], tmp[1]);
    Deofa_L.Omega(tmp[1], tmp[0], -1, 1);
    rw_eofa[hit] = 2.0*Z - k*innerProduct(spProj_Phi,tmp[0]).real();

    // RH term
    for(int s=0; s<Ls; ++s){ axpby_ssp_pplus(spProj_Phi, 0.0, Phi, 1.0, Phi, s, s); }
    Deofa_R.Omega(spProj_Phi, tmp[0], 1, 0);
    G5R5(tmp[1], tmp[0]);
    tmp[0] = zero;
    SchurSolver(Deofa_R, tmp[1], tmp[0]);
    Deofa_R.Dtilde(tmp[0], tmp[1]);
    Deofa_R.Omega(tmp[1], tmp[0], 1, 1);
    rw_eofa[hit] += k*innerProduct(spProj_Phi,tmp[0]).real();
    std::cout << std::endl << "==================================================" << std::endl;
    std::cout << " --- EOFA: Hit " << hit << ": rw = " << rw_eofa[hit];
    std::cout << std::endl << "==================================================" << std::endl << std::endl;

  }

  std::vector<RealD> rhmc_result = jack_stats(rw_rhmc);
  std::vector<RealD> eofa_result = jack_stats(rw_eofa);
  std::cout << std::endl << "RHMC: rw = " << rhmc_result[0] << " +/- " << rhmc_result[1] << std::endl;
  std::cout << std::endl << "EOFA: rw = " << eofa_result[0] << " +/- " << eofa_result[1] << std::endl;

  Grid_finalize();
}
