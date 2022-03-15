/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/core/Test_wilson_exp_clover.cc

    Copyright (C) 2022

    Author: Guido Cossu <guido.cossu@ed.ac.uk>
            Fabian Joswig <fabian.joswig@ed.ac.uk>

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

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  auto latt_size = GridDefaultLatt();
  auto simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  auto mpi_layout = GridDefaultMpi();
  GridCartesian Grid(latt_size, simd_layout, mpi_layout);
  GridRedBlackCartesian RBGrid(&Grid);

  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;
  std::cout << GridLogMessage << "Grid floating point word size is REALF" << sizeof(RealF) << std::endl;
  std::cout << GridLogMessage << "Grid floating point word size is REALD" << sizeof(RealD) << std::endl;
  std::cout << GridLogMessage << "Grid floating point word size is REAL" << sizeof(Real) << std::endl;

  std::vector<int> seeds({1, 2, 3, 4});
  GridParallelRNG pRNG(&Grid);
  pRNG.SeedFixedIntegers(seeds);
  //  pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9});

  typedef typename WilsonExpCloverFermionR::FermionField FermionField;
  typename WilsonExpCloverFermionR::ImplParams params;
  WilsonAnisotropyCoefficients anis;

  FermionField src(&Grid);
  random(pRNG, src);
  FermionField result(&Grid);
  result = Zero();
  FermionField result2(&Grid);
  result2 = Zero();
  FermionField ref(&Grid);
  ref = Zero();
  FermionField tmp(&Grid);
  tmp = Zero();
  FermionField err(&Grid);
  err = Zero();
  FermionField phi(&Grid);
  random(pRNG, phi);
  FermionField chi(&Grid);
  random(pRNG, chi);
  LatticeGaugeField Umu(&Grid);
  SU<Nc>::HotConfiguration(pRNG, Umu);
  std::vector<LatticeColourMatrix> U(4, &Grid);

  double tolerance = 1e-4;

  double volume = 1;
  for (int mu = 0; mu < Nd; mu++)
  {
    volume = volume * latt_size[mu];
  }

  RealD mass = 0.1;
  RealD csw_r = 1.0;
  RealD csw_t = 1.0;

  WilsonExpCloverFermionR Dwc(Umu, Grid, RBGrid, mass, csw_r, csw_t, anis, params);
  CompactWilsonExpCloverFermionR Dwc_compact(Umu, Grid, RBGrid, mass, csw_r, csw_t, 1.0, anis, params);

  std::cout << GridLogMessage << "==========================================================" << std::endl;
  std::cout << GridLogMessage << "= Testing that Deo + Doe = Dunprec " << std::endl;
  std::cout << GridLogMessage << "==========================================================" << std::endl;

  FermionField src_e(&RBGrid);
  FermionField src_o(&RBGrid);
  FermionField r_e(&RBGrid);
  FermionField r_o(&RBGrid);
  FermionField r_eo(&Grid);
  pickCheckerboard(Even, src_e, src);
  pickCheckerboard(Odd, src_o, src);

  Dwc.Meooe(src_e, r_o);
  std::cout << GridLogMessage << "Applied Meo" << std::endl;
  Dwc.Meooe(src_o, r_e);
  std::cout << GridLogMessage << "Applied Moe" << std::endl;
  Dwc.Dhop(src, ref, DaggerNo);

  setCheckerboard(r_eo, r_o);
  setCheckerboard(r_eo, r_e);

  err = ref - r_eo;
  std::cout << GridLogMessage << "EO norm diff\t" << norm2(err) << " (" << norm2(ref) << " - " << norm2(r_eo) << ")" << std::endl;
  assert(fabs(norm2(err)) < tolerance);



  Dwc_compact.Meooe(src_e, r_o);
  std::cout << GridLogMessage << "Applied Meo" << std::endl;
  Dwc_compact.Meooe(src_o, r_e);
  std::cout << GridLogMessage << "Applied Moe" << std::endl;
  Dwc_compact.Dhop(src, ref, DaggerNo);

  setCheckerboard(r_eo, r_o);
  setCheckerboard(r_eo, r_e);

  err = ref - r_eo;
  std::cout << GridLogMessage << "EO norm diff compact\t" << norm2(err) << " (" << norm2(ref) << " - " << norm2(r_eo) << ")" << std::endl;
  assert(fabs(norm2(err)) < tolerance);


  std::cout << GridLogMessage << "==============================================================" << std::endl;
  std::cout << GridLogMessage << "= Test Ddagger is the dagger of D by requiring                " << std::endl;
  std::cout << GridLogMessage << "=  < phi | Deo | chi > * = < chi | Deo^dag| phi>  " << std::endl;
  std::cout << GridLogMessage << "==============================================================" << std::endl;

  FermionField chi_e(&RBGrid);
  FermionField chi_o(&RBGrid);

  FermionField dchi_e(&RBGrid);
  FermionField dchi_o(&RBGrid);

  FermionField phi_e(&RBGrid);
  FermionField phi_o(&RBGrid);

  FermionField dphi_e(&RBGrid);
  FermionField dphi_o(&RBGrid);

  pickCheckerboard(Even, chi_e, chi);
  pickCheckerboard(Odd, chi_o, chi);
  pickCheckerboard(Even, phi_e, phi);
  pickCheckerboard(Odd, phi_o, phi);

  Dwc.Meooe(chi_e, dchi_o);
  Dwc.Meooe(chi_o, dchi_e);
  Dwc.MeooeDag(phi_e, dphi_o);
  Dwc.MeooeDag(phi_o, dphi_e);

  ComplexD pDce = innerProduct(phi_e, dchi_e);
  ComplexD pDco = innerProduct(phi_o, dchi_o);
  ComplexD cDpe = innerProduct(chi_e, dphi_e);
  ComplexD cDpo = innerProduct(chi_o, dphi_o);

  std::cout << GridLogMessage << "e " << pDce << " " << cDpe << std::endl;
  std::cout << GridLogMessage << "o " << pDco << " " << cDpo << std::endl;

  std::cout << GridLogMessage << "pDce - conj(cDpo) " << pDce - conj(cDpo) << std::endl;
  std::cout << GridLogMessage << "pDco - conj(cDpe) " << pDco - conj(cDpe) << std::endl;

  Dwc_compact.Meooe(chi_e, dchi_o);
  Dwc_compact.Meooe(chi_o, dchi_e);
  Dwc_compact.MeooeDag(phi_e, dphi_o);
  Dwc_compact.MeooeDag(phi_o, dphi_e);

  pDce = innerProduct(phi_e, dchi_e);
  pDco = innerProduct(phi_o, dchi_o);
  cDpe = innerProduct(chi_e, dphi_e);
  cDpo = innerProduct(chi_o, dphi_o);

  std::cout << GridLogMessage << "e compact " << pDce << " " << cDpe << std::endl;
  std::cout << GridLogMessage << "o compact " << pDco << " " << cDpo << std::endl;

  std::cout << GridLogMessage << "pDce - conj(cDpo) compact " << pDce - conj(cDpo) << std::endl;
  std::cout << GridLogMessage << "pDco - conj(cDpe) compact " << pDco - conj(cDpe) << std::endl;

  std::cout << GridLogMessage << "==============================================================" << std::endl;
  std::cout << GridLogMessage << "= Test MeeInv Mee = 1   (if csw!=0)                           " << std::endl;
  std::cout << GridLogMessage << "==============================================================" << std::endl;

  pickCheckerboard(Even, chi_e, chi);
  pickCheckerboard(Odd, chi_o, chi);

  Dwc.Mooee(chi_e, src_e);
  Dwc.MooeeInv(src_e, phi_e);

  Dwc.Mooee(chi_o, src_o);
  Dwc.MooeeInv(src_o, phi_o);

  setCheckerboard(phi, phi_e);
  setCheckerboard(phi, phi_o);

  err = phi - chi;
  std::cout << GridLogMessage << "norm diff " << norm2(err) << std::endl;
  assert(fabs(norm2(err)) < tolerance);

  Dwc_compact.Mooee(chi_e, src_e);
  Dwc_compact.MooeeInv(src_e, phi_e);

  Dwc_compact.Mooee(chi_o, src_o);
  Dwc_compact.MooeeInv(src_o, phi_o);

  setCheckerboard(phi, phi_e);
  setCheckerboard(phi, phi_o);

  err = phi - chi;
  std::cout << GridLogMessage << "norm diff compact " << norm2(err) << std::endl;
  assert(fabs(norm2(err)) < tolerance);

  std::cout << GridLogMessage << "==============================================================" << std::endl;
  std::cout << GridLogMessage << "= Test MeeDag MeeInvDag = 1    (if csw!=0)                    " << std::endl;
  std::cout << GridLogMessage << "==============================================================" << std::endl;

  pickCheckerboard(Even, chi_e, chi);
  pickCheckerboard(Odd, chi_o, chi);

  Dwc.MooeeDag(chi_e, src_e);
  Dwc.MooeeInvDag(src_e, phi_e);

  Dwc.MooeeDag(chi_o, src_o);
  Dwc.MooeeInvDag(src_o, phi_o);

  setCheckerboard(phi, phi_e);
  setCheckerboard(phi, phi_o);

  err = phi - chi;
  std::cout << GridLogMessage << "norm diff " << norm2(err) << std::endl;
  assert(fabs(norm2(err)) < tolerance);

  Dwc_compact.MooeeDag(chi_e, src_e);
  Dwc_compact.MooeeInvDag(src_e, phi_e);

  Dwc_compact.MooeeDag(chi_o, src_o);
  Dwc_compact.MooeeInvDag(src_o, phi_o);

  setCheckerboard(phi, phi_e);
  setCheckerboard(phi, phi_o);

  err = phi - chi;
  std::cout << GridLogMessage << "norm diff compact " << norm2(err) << std::endl;
  assert(fabs(norm2(err)) < tolerance);

  std::cout << GridLogMessage << "==============================================================" << std::endl;
  std::cout << GridLogMessage << "= Test MeeInv MeeDag = 1      (if csw!=0)                     " << std::endl;
  std::cout << GridLogMessage << "==============================================================" << std::endl;

  pickCheckerboard(Even, chi_e, chi);
  pickCheckerboard(Odd, chi_o, chi);

  Dwc.MooeeDag(chi_e, src_e);
  Dwc.MooeeInv(src_e, phi_e);

  Dwc.MooeeDag(chi_o, src_o);
  Dwc.MooeeInv(src_o, phi_o);

  setCheckerboard(phi, phi_e);
  setCheckerboard(phi, phi_o);

  err = phi - chi;
  std::cout << GridLogMessage << "norm diff " << norm2(err) << std::endl;
  assert(fabs(norm2(err)) < tolerance);

  Dwc_compact.MooeeDag(chi_e, src_e);
  Dwc_compact.MooeeInv(src_e, phi_e);

  Dwc_compact.MooeeDag(chi_o, src_o);
  Dwc_compact.MooeeInv(src_o, phi_o);

  setCheckerboard(phi, phi_e);
  setCheckerboard(phi, phi_o);

  err = phi - chi;
  std::cout << GridLogMessage << "norm diff compact " << norm2(err) << std::endl;
  assert(fabs(norm2(err)) < tolerance);

  std::cout << GridLogMessage << "================================================================" << std::endl;
  std::cout << GridLogMessage << "= Testing gauge covariance Clover term with EO preconditioning  " << std::endl;
  std::cout << GridLogMessage << "================================================================" << std::endl;

  chi = Zero();
  phi = Zero();
  tmp = Zero();
  pickCheckerboard(Even, chi_e, chi);
  pickCheckerboard(Odd, chi_o, chi);
  pickCheckerboard(Even, phi_e, phi);
  pickCheckerboard(Odd, phi_o, phi);

  Dwc.Mooee(src_e, chi_e);
  Dwc.Mooee(src_o, chi_o);
  setCheckerboard(chi, chi_e);
  setCheckerboard(chi, chi_o);
  setCheckerboard(src, src_e);
  setCheckerboard(src, src_o);

  ////////////////////// Gauge Transformation
  std::vector<int> seeds2({5, 6, 7, 8});
  GridParallelRNG pRNG2(&Grid);
  pRNG2.SeedFixedIntegers(seeds2);
  LatticeColourMatrix Omega(&Grid);
  LatticeColourMatrix ShiftedOmega(&Grid);
  LatticeGaugeField U_prime(&Grid);
  U_prime = Zero();
  LatticeColourMatrix U_prime_mu(&Grid);
  U_prime_mu = Zero();
  SU<Nc>::LieRandomize(pRNG2, Omega, 1.0);
  for (int mu = 0; mu < Nd; mu++)
  {
    U[mu] = peekLorentz(Umu, mu);
    ShiftedOmega = Cshift(Omega, mu, 1);
    U_prime_mu = Omega * U[mu] * adj(ShiftedOmega);
    pokeLorentz(U_prime, U_prime_mu, mu);
  }
  /////////////////

  WilsonExpCloverFermionR Dwc_prime(U_prime, Grid, RBGrid, mass, csw_r, csw_t, anis, params);
  CompactWilsonExpCloverFermionR Dwc_compact_prime(U_prime, Grid, RBGrid, mass, csw_r, csw_t, 1.0, anis, params);

  tmp = Omega * src;
  pickCheckerboard(Even, src_e, tmp);
  pickCheckerboard(Odd, src_o, tmp);

  Dwc_prime.Mooee(src_e, phi_e);
  Dwc_prime.Mooee(src_o, phi_o);

  setCheckerboard(phi, phi_e);
  setCheckerboard(phi, phi_o);

  err = chi - adj(Omega) * phi;
  std::cout << GridLogMessage << "norm diff " << norm2(err) << std::endl;
  assert(fabs(norm2(err)) < tolerance);

  chi = Zero();
  phi = Zero();
  tmp = Zero();
  pickCheckerboard(Even, chi_e, chi);
  pickCheckerboard(Odd, chi_o, chi);
  pickCheckerboard(Even, phi_e, phi);
  pickCheckerboard(Odd, phi_o, phi);

  Dwc_compact.Mooee(src_e, chi_e);
  Dwc_compact.Mooee(src_o, chi_o);
  setCheckerboard(chi, chi_e);
  setCheckerboard(chi, chi_o);
  setCheckerboard(src, src_e);
  setCheckerboard(src, src_o);

  tmp = Omega * src;
  pickCheckerboard(Even, src_e, tmp);
  pickCheckerboard(Odd, src_o, tmp);

  Dwc_compact_prime.Mooee(src_e, phi_e);
  Dwc_compact_prime.Mooee(src_o, phi_o);

  setCheckerboard(phi, phi_e);
  setCheckerboard(phi, phi_o);

  err = chi - adj(Omega) * phi;
  std::cout << GridLogMessage << "norm diff compact " << norm2(err) << std::endl;
  assert(fabs(norm2(err)) < tolerance);

  std::cout << GridLogMessage << "=================================================================" << std::endl;
  std::cout << GridLogMessage << "= Testing gauge covariance Clover term w/o EO preconditioning  " << std::endl;
  std::cout << GridLogMessage << "================================================================" << std::endl;

  chi = Zero();
  phi = Zero();

  WilsonFermionR Dw(Umu, Grid, RBGrid, mass, params);

  Dw.M(src, result);
  Dwc.M(src, chi);

  Dwc_prime.M(Omega * src, phi);

  WilsonFermionR Dw_prime(U_prime, Grid, RBGrid, mass, params);
  Dw_prime.M(Omega * src, result2);

  err = result - adj(Omega) * result2;
  std::cout << GridLogMessage << "norm diff Wilson                 " << norm2(err) << std::endl;
  assert(fabs(norm2(err)) < tolerance);
  err = chi - adj(Omega) * phi;
  std::cout << GridLogMessage << "norm diff WilsonExpClover        " << norm2(err) << std::endl;
  assert(fabs(norm2(err)) < tolerance);

  chi = Zero();
  phi = Zero();

  Dwc_compact.M(src, chi);
  Dwc_compact_prime.M(Omega * src, phi);

  err = chi - adj(Omega) * phi;
  std::cout << GridLogMessage << "norm diff CompactWilsonExpClover " << norm2(err) << std::endl;
  assert(fabs(norm2(err)) < tolerance);

  std::cout << GridLogMessage << "==========================================================" << std::endl;
  std::cout << GridLogMessage << "= Testing Mooee(csw=0) Clover to reproduce Mooee Wilson   " << std::endl;
  std::cout << GridLogMessage << "==========================================================" << std::endl;

  chi = Zero();
  phi = Zero();
  err = Zero();
  WilsonExpCloverFermionR Dwc_csw0(Umu, Grid, RBGrid, mass, 0.0, 0.0, anis, params); //  <-- Notice: csw=0

  pickCheckerboard(Even, phi_e, phi);
  pickCheckerboard(Odd, phi_o, phi);
  pickCheckerboard(Even, chi_e, chi);
  pickCheckerboard(Odd, chi_o, chi);

  Dw.Mooee(src_e, chi_e);
  Dw.Mooee(src_o, chi_o);
  Dwc_csw0.Mooee(src_e, phi_e);
  Dwc_csw0.Mooee(src_o, phi_o);

  setCheckerboard(chi, chi_e);
  setCheckerboard(chi, chi_o);
  setCheckerboard(phi, phi_e);
  setCheckerboard(phi, phi_o);
  setCheckerboard(src, src_e);
  setCheckerboard(src, src_o);

  err = chi - phi;
  std::cout << GridLogMessage << "norm diff " << norm2(err) << std::endl;
  assert(fabs(norm2(err)) < tolerance);

  chi = Zero();
  phi = Zero();
  err = Zero();
  CompactWilsonExpCloverFermionR Dwc_compact_csw0(Umu, Grid, RBGrid, mass, 0.0, 0.0, 1.0, anis, params); //  <-- Notice: csw=0

  pickCheckerboard(Even, phi_e, phi);
  pickCheckerboard(Odd, phi_o, phi);
  pickCheckerboard(Even, chi_e, chi);
  pickCheckerboard(Odd, chi_o, chi);

  Dw.Mooee(src_e, chi_e);
  Dw.Mooee(src_o, chi_o);
  Dwc_compact_csw0.Mooee(src_e, phi_e);
  Dwc_compact_csw0.Mooee(src_o, phi_o);

  setCheckerboard(chi, chi_e);
  setCheckerboard(chi, chi_o);
  setCheckerboard(phi, phi_e);
  setCheckerboard(phi, phi_o);
  setCheckerboard(src, src_e);
  setCheckerboard(src, src_o);

  err = chi - phi;
  std::cout << GridLogMessage << "norm diff compact " << norm2(err) << std::endl;
  assert(fabs(norm2(err)) < tolerance);

  std::cout << GridLogMessage << "==========================================================" << std::endl;
  std::cout << GridLogMessage << "= Testing EO operator is equal to the unprec              " << std::endl;
  std::cout << GridLogMessage << "==========================================================" << std::endl;

  chi = Zero();
  phi = Zero();
  err = Zero();

  pickCheckerboard(Even, phi_e, phi);
  pickCheckerboard(Odd, phi_o, phi);
  pickCheckerboard(Even, chi_e, chi);
  pickCheckerboard(Odd, chi_o, chi);

  // M phi = (Mooee src_e + Meooe src_o , Meooe src_e + Mooee src_o)

  Dwc.M(src, ref); // Reference result from the unpreconditioned operator

  // EO matrix
  Dwc.Mooee(src_e, chi_e); 
  Dwc.Mooee(src_o, chi_o);
  Dwc.Meooe(src_o, phi_e);
  Dwc.Meooe(src_e, phi_o);

  phi_o += chi_o;
  phi_e += chi_e;

  setCheckerboard(phi, phi_e);
  setCheckerboard(phi, phi_o);

  err = ref - phi;
  std::cout << GridLogMessage << "ref (unpreconditioned operator) diff         : " << norm2(ref) << std::endl;
  std::cout << GridLogMessage << "phi (EO decomposition)          diff         : " << norm2(phi) << std::endl;
  std::cout << GridLogMessage << "norm diff                                    : " << norm2(err) << std::endl;
  assert(fabs(norm2(err)) < tolerance);

  chi = Zero();
  phi = Zero();
  err = Zero();

  pickCheckerboard(Even, phi_e, phi);
  pickCheckerboard(Odd, phi_o, phi);
  pickCheckerboard(Even, chi_e, chi);
  pickCheckerboard(Odd, chi_o, chi);

  // M phi = (Mooee src_e + Meooe src_o , Meooe src_e + Mooee src_o)

  Dwc_compact.M(src, ref); // Reference result from the unpreconditioned operator

  // EO matrix
  Dwc_compact.Mooee(src_e, chi_e); 
  Dwc_compact.Mooee(src_o, chi_o);
  Dwc_compact.Meooe(src_o, phi_e);
  Dwc_compact.Meooe(src_e, phi_o);

  phi_o += chi_o;
  phi_e += chi_e;

  setCheckerboard(phi, phi_e);
  setCheckerboard(phi, phi_o);

  err = ref - phi;
  std::cout << GridLogMessage << "ref (unpreconditioned operator) diff compact : " << norm2(ref) << std::endl;
  std::cout << GridLogMessage << "phi (EO decomposition)          diff compact : " << norm2(phi) << std::endl;
  std::cout << GridLogMessage << "norm diff compact                            : " << norm2(err) << std::endl;
  assert(fabs(norm2(err)) < tolerance);

  Grid_finalize();
}
