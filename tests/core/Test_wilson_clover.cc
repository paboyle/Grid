/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./benchmarks/Benchmark_wilson.cc

    Copyright (C) 2015

    Author: Guido Cossu <guido.cossu@ed.ac.uk>

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

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  std::vector<int> latt_size = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  std::vector<int> mpi_layout = GridDefaultMpi();
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

  typedef typename WilsonCloverFermionR::FermionField FermionField;
  typename WilsonCloverFermionR::ImplParams params;
  WilsonAnisotropyCoefficients anis;

  FermionField src(&Grid);
  random(pRNG, src);
  FermionField result(&Grid);
  result = zero;
  FermionField result2(&Grid);
  result2 = zero;
  FermionField ref(&Grid);
  ref = zero;
  FermionField tmp(&Grid);
  tmp = zero;
  FermionField err(&Grid);
  err = zero;
  FermionField err2(&Grid);
  err2 = zero;
  FermionField phi(&Grid);
  random(pRNG, phi);
  FermionField chi(&Grid);
  random(pRNG, chi);
  LatticeGaugeField Umu(&Grid);
  SU3::HotConfiguration(pRNG, Umu);
  std::vector<LatticeColourMatrix> U(4, &Grid);

  double volume = 1;
  for (int mu = 0; mu < Nd; mu++)
  {
    volume = volume * latt_size[mu];
  }

  RealD mass = 0.1;
  RealD csw_r = 1.0;
  RealD csw_t = 1.0;

  WilsonCloverFermionR Dwc(Umu, Grid, RBGrid, mass, csw_r, csw_t, anis, params);
  //Dwc.ImportGauge(Umu); // not necessary, included in the constructor

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
  std::cout << GridLogMessage << "EO norm diff   " << norm2(err) << " " << norm2(ref) << " " << norm2(r_eo) << std::endl;

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
  std::cout << GridLogMessage << "norm diff   " << norm2(err) << std::endl;

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
  std::cout << GridLogMessage << "norm diff   " << norm2(err) << std::endl;

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
  std::cout << GridLogMessage << "norm diff   " << norm2(err) << std::endl;

  std::cout << GridLogMessage << "================================================================" << std::endl;
  std::cout << GridLogMessage << "= Testing gauge covariance Clover term with EO preconditioning  " << std::endl;
  std::cout << GridLogMessage << "================================================================" << std::endl;

  chi = zero;
  phi = zero;
  tmp = zero;
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
  U_prime = zero;
  LatticeColourMatrix U_prime_mu(&Grid);
  U_prime_mu = zero;
  SU<Nc>::LieRandomize(pRNG2, Omega, 1.0);
  for (int mu = 0; mu < Nd; mu++)
  {
    U[mu] = peekLorentz(Umu, mu);
    ShiftedOmega = Cshift(Omega, mu, 1);
    U_prime_mu = Omega * U[mu] * adj(ShiftedOmega);
    pokeLorentz(U_prime, U_prime_mu, mu);
  }
  /////////////////

  WilsonCloverFermionR Dwc_prime(U_prime, Grid, RBGrid, mass, csw_r, csw_t, anis, params);
  Dwc_prime.ImportGauge(U_prime);

  tmp = Omega * src;
  pickCheckerboard(Even, src_e, tmp);
  pickCheckerboard(Odd, src_o, tmp);

  Dwc_prime.Mooee(src_e, phi_e);
  Dwc_prime.Mooee(src_o, phi_o);

  setCheckerboard(phi, phi_e);
  setCheckerboard(phi, phi_o);

  err = chi - adj(Omega) * phi;
  std::cout << GridLogMessage << "norm diff   " << norm2(err) << std::endl;

  std::cout << GridLogMessage << "=================================================================" << std::endl;
  std::cout << GridLogMessage << "= Testing gauge covariance Clover term w/o EO preconditioning  " << std::endl;
  std::cout << GridLogMessage << "================================================================" << std::endl;

  chi = zero;
  phi = zero;

  WilsonFermionR Dw(Umu, Grid, RBGrid, mass, params);
  Dw.ImportGauge(Umu);

  Dw.M(src, result);
  Dwc.M(src, chi);

  Dwc_prime.M(Omega * src, phi);

  WilsonFermionR Dw_prime(U_prime, Grid, RBGrid, mass, params);
  Dw_prime.ImportGauge(U_prime);
  Dw_prime.M(Omega * src, result2);

  err = chi - adj(Omega) * phi;
  err2 = result - adj(Omega) * result2;
  std::cout << GridLogMessage << "norm diff Wilson   " << norm2(err) << std::endl;
  std::cout << GridLogMessage << "norm diff WilsonClover  " << norm2(err2) << std::endl;

  std::cout << GridLogMessage << "==========================================================" << std::endl;
  std::cout << GridLogMessage << "= Testing Mooee(csw=0) Clover to reproduce Mooee Wilson   " << std::endl;
  std::cout << GridLogMessage << "==========================================================" << std::endl;

  chi = zero;
  phi = zero;
  err = zero;
  WilsonCloverFermionR Dwc_csw0(Umu, Grid, RBGrid, mass, 0.0, 0.0, anis, params); //  <-- Notice: csw=0
  Dwc_csw0.ImportGauge(Umu);

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
  std::cout << GridLogMessage << "norm diff  " << norm2(err) << std::endl;

  std::cout << GridLogMessage << "==========================================================" << std::endl;
  std::cout << GridLogMessage << "= Testing EO operator is equal to the unprec              " << std::endl;
  std::cout << GridLogMessage << "==========================================================" << std::endl;

  chi = zero;
  phi = zero;
  err = zero;

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
  std::cout << GridLogMessage << "ref (unpreconditioned operator) diff  :" << norm2(ref) << std::endl;
  std::cout << GridLogMessage << "phi (EO decomposition)          diff  :" << norm2(phi) << std::endl;
  std::cout << GridLogMessage << "norm diff                             :" << norm2(err) << std::endl;

  Grid_finalize();
}
