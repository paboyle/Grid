/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/solver/Test_coarse_even_odd.cc

    Copyright (C) 2015-2020

    Author: Daniel Richtmann <daniel.richtmann@ur.de>

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

using namespace Grid;

#ifndef NBASIS
#define NBASIS 40
#endif

// NOTE: The tests in this file are written in analogy to
// - tests/core/Test_wilson_even_odd.cc
// - tests/core/Test_wilson_clover.cc

std::vector<int> readFromCommandlineIvec(int*                    argc,
                                         char***                 argv,
                                         std::string&&           option,
                                         const std::vector<int>& defaultValue) {
  std::string      arg;
  std::vector<int> ret(defaultValue);
  if(GridCmdOptionExists(*argv, *argv + *argc, option)) {
    arg = GridCmdOptionPayload(*argv, *argv + *argc, option);
    GridCmdOptionIntVector(arg, ret);
  }
  return ret;
}

int main(int argc, char** argv) {
  Grid_init(&argc, &argv);

  /////////////////////////////////////////////////////////////////////////////
  //                          Read from command line                         //
  /////////////////////////////////////////////////////////////////////////////

  const int  nbasis    = NBASIS; static_assert((nbasis & 0x1) == 0, "");
  const int  nb        = nbasis/2;
  Coordinate blockSize = readFromCommandlineIvec(&argc, &argv, "--blocksize", {2, 2, 2, 2});

  std::cout << GridLogMessage << "Compiled with nbasis = " << nbasis << " -> nb = " << nb << std::endl;

  /////////////////////////////////////////////////////////////////////////////
  //                              General setup                              //
  /////////////////////////////////////////////////////////////////////////////

  Coordinate clatt = GridDefaultLatt();
  for(int d=0; d<clatt.size(); d++) clatt[d] = clatt[d] / blockSize[d];

  GridCartesian*         Grid_f   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  GridCartesian*         Grid_c   = SpaceTimeGrid::makeFourDimGrid(clatt, GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* RBGrid_f = SpaceTimeGrid::makeFourDimRedBlackGrid(Grid_f);
  GridRedBlackCartesian* RBGrid_c = SpaceTimeGrid::makeFourDimRedBlackGrid(Grid_c);

  std::cout << GridLogMessage << "Grid_f:" << std::endl; Grid_f->show_decomposition();
  std::cout << GridLogMessage << "Grid_c:" << std::endl; Grid_c->show_decomposition();
  std::cout << GridLogMessage << "RBGrid_f:" << std::endl; RBGrid_f->show_decomposition();
  std::cout << GridLogMessage << "RBGrid_c:" << std::endl; RBGrid_c->show_decomposition();

  GridParallelRNG pRNG_f(Grid_f);
  GridParallelRNG pRNG_c(Grid_c);

  std::vector<int> seeds({1, 2, 3, 4});

  pRNG_f.SeedFixedIntegers(seeds);
  pRNG_c.SeedFixedIntegers(seeds);

  /////////////////////////////////////////////////////////////////////////////
  //                    Setup of Dirac Matrix and Operator                   //
  /////////////////////////////////////////////////////////////////////////////

  LatticeGaugeField Umu(Grid_f);
#if (Nc==2)
  SU2::HotConfiguration(pRNG_f, Umu);
#elif (defined Nc==3)
  SU3::HotConfiguration(pRNG_f, Umu);
#elif (defined Nc==4)
  SU4::HotConfiguration(pRNG_f, Umu);
#elif (defined Nc==5)
  SU5::HotConfiguration(pRNG_f, Umu);
#endif
  RealD checkTolerance = (getPrecision<LatticeFermion>::value == 1) ? 1e-7 : 1e-15;

  RealD mass = -0.30;
  RealD csw  = 1.9192;

  WilsonCloverFermionR Dwc(Umu, *Grid_f, *RBGrid_f, mass, csw, csw);
  MdagMLinearOperator<WilsonCloverFermionR, LatticeFermion> MdagMOp_Dwc(Dwc);

  /////////////////////////////////////////////////////////////////////////////
  //                             Type definitions                            //
  /////////////////////////////////////////////////////////////////////////////

  typedef Aggregation<vSpinColourVector, vTComplex, nbasis>     Aggregates;
  typedef CoarsenedMatrix<vSpinColourVector, vTComplex, nbasis> CoarseDiracMatrix;
  typedef CoarseDiracMatrix::CoarseVector                       CoarseVector;

  /////////////////////////////////////////////////////////////////////////////
  //                           Setup of Aggregation                          //
  /////////////////////////////////////////////////////////////////////////////

  Aggregates Aggs(Grid_c, Grid_f, 0);
  {
    LatticeFermion tmp(Aggs.subspace[0].Grid());
    for(int n = 0; n < nb; n++) {
      gaussian(pRNG_f, Aggs.subspace[n]);
      G5C(tmp, Aggs.subspace[n]);
      axpby(Aggs.subspace[n + nb], 0.5, -0.5, Aggs.subspace[n], tmp);
      axpby(Aggs.subspace[n], 0.5, 0.5, Aggs.subspace[n], tmp);
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  //                  Setup of CoarsenedMatrix and Operator                  //
  /////////////////////////////////////////////////////////////////////////////

  const int hermitian = 0;
  CoarseDiracMatrix Dc(*Grid_c, *RBGrid_c, hermitian);
  Dc.CoarsenOperator(Grid_f, MdagMOp_Dwc, Aggs);
  MdagMLinearOperator<CoarseDiracMatrix, CoarseVector> MdagMOp_Dc(Dc);

  /////////////////////////////////////////////////////////////////////////////
  //                     Setup vectors used in all tests                     //
  /////////////////////////////////////////////////////////////////////////////

  CoarseVector src(Grid_c);  random(pRNG_c, src);
  CoarseVector diff(Grid_c); diff = Zero();

  /////////////////////////////////////////////////////////////////////////////
  //                              Start of tests                             //
  /////////////////////////////////////////////////////////////////////////////

  {
    std::cout << GridLogMessage << "===========================================================================" << std::endl;
    std::cout << GridLogMessage << "= Test Dhop + Mdiag = Munprec" << std::endl;
    std::cout << GridLogMessage << "===========================================================================" << std::endl;

    CoarseVector phi(Grid_c); phi = Zero();
    CoarseVector chi(Grid_c); chi = Zero();
    CoarseVector res(Grid_c); res = Zero();
    CoarseVector ref(Grid_c); ref = Zero();

    Dc.Mdiag(src, phi);          std::cout << GridLogMessage << "Applied Mdiag" << std::endl;
    Dc.Dhop(src, chi, DaggerNo); std::cout << GridLogMessage << "Applied Dhop"  << std::endl;
    Dc.M(src, ref);              std::cout << GridLogMessage << "Applied M"     << std::endl;

    res = phi + chi;

    diff = ref - res;
    auto absDev = norm2(diff);
    auto relDev = absDev / norm2(ref);
    std::cout << GridLogMessage << "norm2(Munprec), norm2(Dhop + Mdiag), abs. deviation, rel. deviation: "
              << norm2(ref) << " " << norm2(res) << " " << absDev << " " << relDev << " -> check "
              << ((relDev < checkTolerance) ? "passed" : "failed") << std::endl;
    assert(relDev <= checkTolerance);
  }

  {
    std::cout << GridLogMessage << "===========================================================================" << std::endl;
    std::cout << GridLogMessage << "= Test Meo + Moe = Dhop" << std::endl;
    std::cout << GridLogMessage << "===========================================================================" << std::endl;

    CoarseVector src_e(RBGrid_c); src_e = Zero();
    CoarseVector src_o(RBGrid_c); src_o = Zero();
    CoarseVector res_e(RBGrid_c); res_e = Zero();
    CoarseVector res_o(RBGrid_c); res_o = Zero();
    CoarseVector res(Grid_c);     res = Zero();
    CoarseVector ref(Grid_c);     ref = Zero();

    pickCheckerboard(Even, src_e, src);
    pickCheckerboard(Odd,  src_o, src);

    Dc.Meooe(src_e, res_o);        std::cout << GridLogMessage << "Applied Meo"  << std::endl;
    Dc.Meooe(src_o, res_e);        std::cout << GridLogMessage << "Applied Moe"  << std::endl;
    Dc.Dhop(src, ref, DaggerNo); std::cout << GridLogMessage << "Applied Dhop" << std::endl;

    setCheckerboard(res, res_o);
    setCheckerboard(res, res_e);

    diff = ref - res;
    auto absDev = norm2(diff);
    auto relDev = absDev / norm2(ref);
    std::cout << GridLogMessage << "norm2(Dhop), norm2(Meo + Moe), abs. deviation, rel. deviation: "
              << norm2(ref) << " " << norm2(res) << " " << absDev << " " << relDev
              << " -> check " << ((relDev < checkTolerance) ? "passed" : "failed") << std::endl;
    assert(relDev <= checkTolerance);
  }

  {
    std::cout << GridLogMessage << "===========================================================================" << std::endl;
    std::cout << GridLogMessage << "= Test |(Im(v^dag M^dag M v)| = 0" << std::endl;
    std::cout << GridLogMessage << "===========================================================================" << std::endl;

    CoarseVector tmp(Grid_c); tmp = Zero();
    CoarseVector phi(Grid_c); phi = Zero();

    Dc.M(src, tmp);    std::cout << GridLogMessage << "Applied M"    << std::endl;
    Dc.Mdag(tmp, phi); std::cout << GridLogMessage << "Applied Mdag" << std::endl;

    std::cout << GridLogMessage << "src = " << norm2(src) << " tmp = " << norm2(tmp) << " phi = " << norm2(phi) << std::endl;

    ComplexD dot = innerProduct(src, phi);

    auto relDev = abs(imag(dot)) / abs(real(dot));
    std::cout << GridLogMessage << "Re(v^dag M^dag M v), Im(v^dag M^dag M v), rel.deviation: "
              << real(dot) << " " << imag(dot) << " " << relDev
              << " -> check " << ((relDev < checkTolerance) ? "passed" : "failed") << std::endl;
    assert(relDev <= checkTolerance);
  }

  {
    std::cout << GridLogMessage << "===========================================================================" << std::endl;
    std::cout << GridLogMessage << "= Test |(Im(v^dag Mooee^dag Mooee v)| = 0 (full grid)" << std::endl;
    std::cout << GridLogMessage << "===========================================================================" << std::endl;

    CoarseVector tmp(Grid_c); tmp = Zero();
    CoarseVector phi(Grid_c); phi = Zero();

    Dc.Mooee(src, tmp);    std::cout << GridLogMessage << "Applied Mooee"    << std::endl;
    Dc.MooeeDag(tmp, phi); std::cout << GridLogMessage << "Applied MooeeDag" << std::endl;

    ComplexD dot = innerProduct(src, phi);

    auto relDev = abs(imag(dot)) / abs(real(dot));
    std::cout << GridLogMessage << "Re(v^dag Mooee^dag Mooee v), Im(v^dag Mooee^dag Mooee v), rel.deviation: "
              << real(dot) << " " << imag(dot) << " " << relDev
              << " -> check " << ((relDev < checkTolerance) ? "passed" : "failed") << std::endl;
    assert(relDev <= checkTolerance);
  }

  {
    std::cout << GridLogMessage << "===========================================================================" << std::endl;
    std::cout << GridLogMessage << "= Test MooeeInv Mooee = 1 (full grid)" << std::endl;
    std::cout << GridLogMessage << "===========================================================================" << std::endl;

    CoarseVector tmp(Grid_c); tmp = Zero();
    CoarseVector phi(Grid_c); phi = Zero();

    Dc.Mooee(src, tmp);    std::cout << GridLogMessage << "Applied Mooee"    << std::endl;
    Dc.MooeeInv(tmp, phi); std::cout << GridLogMessage << "Applied MooeeInv" << std::endl;

    diff        = src - phi;
    auto absDev = norm2(diff);
    auto relDev = absDev / norm2(src);
    std::cout << GridLogMessage << "norm2(src), norm2(MooeeInv Mooee src), abs. deviation, rel. deviation: "
              << norm2(src) << " " << norm2(phi) << " " << absDev << " " << relDev
              << " -> check " << ((relDev < checkTolerance) ? "passed" : "failed") << std::endl;
    assert(relDev <= checkTolerance);
  }

  {
    std::cout << GridLogMessage << "===========================================================================" << std::endl;
    std::cout << GridLogMessage << "= Test MeooeDagger is the dagger of Meooe by requiring" << std::endl;
    std::cout << GridLogMessage << "=  < phi | Meooe | chi > * = < chi | Meooe^dag| phi>" << std::endl;
    std::cout << GridLogMessage << "===========================================================================" << std::endl;

    // clang-format off
    CoarseVector phi(Grid_c); random(pRNG_c, phi);
    CoarseVector chi(Grid_c); random(pRNG_c, chi);
    CoarseVector chi_e(RBGrid_c);   chi_e = Zero();
    CoarseVector chi_o(RBGrid_c);   chi_o = Zero();
    CoarseVector dchi_e(RBGrid_c); dchi_e = Zero();
    CoarseVector dchi_o(RBGrid_c); dchi_o = Zero();
    CoarseVector phi_e(RBGrid_c);   phi_e = Zero();
    CoarseVector phi_o(RBGrid_c);   phi_o = Zero();
    CoarseVector dphi_e(RBGrid_c); dphi_e = Zero();
    CoarseVector dphi_o(RBGrid_c); dphi_o = Zero();
    // clang-format on

    pickCheckerboard(Even, chi_e, chi);
    pickCheckerboard(Odd,  chi_o, chi);
    pickCheckerboard(Even, phi_e, phi);
    pickCheckerboard(Odd,  phi_o, phi);

    Dc.Meooe(chi_e, dchi_o);    std::cout << GridLogMessage << "Applied Meo"    << std::endl;
    Dc.Meooe(chi_o, dchi_e);    std::cout << GridLogMessage << "Applied Moe"    << std::endl;
    Dc.MeooeDag(phi_e, dphi_o); std::cout << GridLogMessage << "Applied MeoDag" << std::endl;
    Dc.MeooeDag(phi_o, dphi_e); std::cout << GridLogMessage << "Applied MoeDag" << std::endl;

    ComplexD phiDchi_e = innerProduct(phi_e, dchi_e);
    ComplexD phiDchi_o = innerProduct(phi_o, dchi_o);
    ComplexD chiDphi_e = innerProduct(chi_e, dphi_e);
    ComplexD chiDphi_o = innerProduct(chi_o, dphi_o);

    std::cout << GridLogDebug << "norm dchi_e = " << norm2(dchi_e) << " norm dchi_o = " << norm2(dchi_o) << " norm dphi_e = " << norm2(dphi_e)
              << " norm dphi_o = " << norm2(dphi_e) << std::endl;

    std::cout << GridLogMessage << "e " << phiDchi_e << " " << chiDphi_e << std::endl;
    std::cout << GridLogMessage << "o " << phiDchi_o << " " << chiDphi_o << std::endl;

    std::cout << GridLogMessage << "phiDchi_e - conj(chiDphi_o) " << phiDchi_e - conj(chiDphi_o) << std::endl;
    std::cout << GridLogMessage << "phiDchi_o - conj(chiDphi_e) " << phiDchi_o - conj(chiDphi_e) << std::endl;
  }

  {
    std::cout << GridLogMessage << "===========================================================================" << std::endl;
    std::cout << GridLogMessage << "= Test MooeeInv Mooee = 1 (checkerboards separately)" << std::endl;
    std::cout << GridLogMessage << "===========================================================================" << std::endl;

    CoarseVector chi(Grid_c); random(pRNG_c, chi);
    CoarseVector tmp(Grid_c); tmp = Zero();
    CoarseVector phi(Grid_c); phi = Zero();
    CoarseVector chi_e(RBGrid_c); chi_e = Zero();
    CoarseVector chi_o(RBGrid_c); chi_o = Zero();
    CoarseVector phi_e(RBGrid_c); phi_e = Zero();
    CoarseVector phi_o(RBGrid_c); phi_o = Zero();
    CoarseVector tmp_e(RBGrid_c); tmp_e = Zero();
    CoarseVector tmp_o(RBGrid_c); tmp_o = Zero();

    pickCheckerboard(Even, chi_e, chi);
    pickCheckerboard(Odd,  chi_o, chi);
    pickCheckerboard(Even, tmp_e, tmp);
    pickCheckerboard(Odd,  tmp_o, tmp);

    Dc.Mooee(chi_e, tmp_e);    std::cout << GridLogMessage << "Applied Mee"    << std::endl;
    Dc.MooeeInv(tmp_e, phi_e); std::cout << GridLogMessage << "Applied MeeInv" << std::endl;
    Dc.Mooee(chi_o, tmp_o);    std::cout << GridLogMessage << "Applied Moo"    << std::endl;
    Dc.MooeeInv(tmp_o, phi_o); std::cout << GridLogMessage << "Applied MooInv" << std::endl;

    setCheckerboard(phi, phi_e);
    setCheckerboard(phi, phi_o);

    diff = chi - phi;
    auto absDev = norm2(diff);
    auto relDev = absDev / norm2(chi);
    std::cout << GridLogMessage << "norm2(chi), norm2(MeeInv Mee chi), abs. deviation, rel. deviation: "
              << norm2(chi) << " " << norm2(phi) << " " << absDev << " " << relDev
              << " -> check " << ((relDev < checkTolerance) ? "passed" : "failed") << std::endl;
    assert(relDev <= checkTolerance);
  }

  {
    std::cout << GridLogMessage << "===========================================================================" << std::endl;
    std::cout << GridLogMessage << "= Test MooeeDag MooeeInvDag = 1 (checkerboards separately)" << std::endl;
    std::cout << GridLogMessage << "===========================================================================" << std::endl;

    CoarseVector chi(Grid_c); random(pRNG_c, chi);
    CoarseVector tmp(Grid_c); tmp = Zero();
    CoarseVector phi(Grid_c); phi = Zero();
    CoarseVector chi_e(RBGrid_c); chi_e = Zero();
    CoarseVector chi_o(RBGrid_c); chi_o = Zero();
    CoarseVector phi_e(RBGrid_c); phi_e = Zero();
    CoarseVector phi_o(RBGrid_c); phi_o = Zero();
    CoarseVector tmp_e(RBGrid_c); tmp_e = Zero();
    CoarseVector tmp_o(RBGrid_c); tmp_o = Zero();

    pickCheckerboard(Even, chi_e, chi);
    pickCheckerboard(Odd,  chi_o, chi);
    pickCheckerboard(Even, tmp_e, tmp);
    pickCheckerboard(Odd,  tmp_o, tmp);

    Dc.MooeeDag(chi_e, tmp_e);    std::cout << GridLogMessage << "Applied MeeDag"    << std::endl;
    Dc.MooeeInvDag(tmp_e, phi_e); std::cout << GridLogMessage << "Applied MeeInvDag" << std::endl;
    Dc.MooeeDag(chi_o, tmp_o);    std::cout << GridLogMessage << "Applied MooDag"    << std::endl;
    Dc.MooeeInvDag(tmp_o, phi_o); std::cout << GridLogMessage << "Applied MooInvDag" << std::endl;

    setCheckerboard(phi, phi_e);
    setCheckerboard(phi, phi_o);

    diff = chi - phi;
    auto absDev = norm2(diff);
    auto relDev = absDev / norm2(chi);
    std::cout << GridLogMessage << "norm2(chi), norm2(MeeDag MeeInvDag chi), abs. deviation, rel. deviation: "
              << norm2(chi) << " " << norm2(phi) << " " << absDev << " " << relDev
              << " -> check " << ((relDev < checkTolerance) ? "passed" : "failed") << std::endl;
    assert(relDev <= checkTolerance);
  }

  {
    std::cout << GridLogMessage << "===========================================================================" << std::endl;
    std::cout << GridLogMessage << "= Test Meo + Moe + Moo + Mee = Munprec" << std::endl;
    std::cout << GridLogMessage << "===========================================================================" << std::endl;

    CoarseVector chi(Grid_c); chi = Zero();
    CoarseVector phi(Grid_c); phi = Zero();
    CoarseVector ref(Grid_c); ref = Zero();
    CoarseVector src_e(RBGrid_c); src_e = Zero();
    CoarseVector src_o(RBGrid_c); src_o = Zero();
    CoarseVector phi_e(RBGrid_c); phi_e = Zero();
    CoarseVector phi_o(RBGrid_c); phi_o = Zero();
    CoarseVector chi_e(RBGrid_c); chi_e = Zero();
    CoarseVector chi_o(RBGrid_c); chi_o = Zero();

    pickCheckerboard(Even, src_e, src);
    pickCheckerboard(Odd,  src_o, src);
    pickCheckerboard(Even, phi_e, phi);
    pickCheckerboard(Odd,  phi_o, phi);
    pickCheckerboard(Even, chi_e, chi);
    pickCheckerboard(Odd,  chi_o, chi);

    // M phi = (Mooee src_e + Meooe src_o , Mooee src_o + Meooe src_e)

    Dc.M(src, ref); // Reference result from the unpreconditioned operator

    // EO matrix
    Dc.Mooee(src_e, chi_e); std::cout << GridLogMessage << "Applied Mee" << std::endl;
    Dc.Mooee(src_o, chi_o); std::cout << GridLogMessage << "Applied Moo" << std::endl;
    Dc.Meooe(src_o, phi_e); std::cout << GridLogMessage << "Applied Moe" << std::endl;
    Dc.Meooe(src_e, phi_o); std::cout << GridLogMessage << "Applied Meo" << std::endl;

    phi_o += chi_o;
    phi_e += chi_e;

    setCheckerboard(phi, phi_e);
    setCheckerboard(phi, phi_o);

    std::cout << GridLogDebug << "norm phi_e = " << norm2(phi_e) << " norm phi_o = " << norm2(phi_o) << " norm phi = " << norm2(phi) << std::endl;

    diff = ref - phi;
    auto absDev = norm2(diff);
    auto relDev = absDev / norm2(ref);
    std::cout << GridLogMessage << "norm2(Dunprec), norm2(Deoprec), abs. deviation, rel. deviation: "
              << norm2(ref) << " " << norm2(phi) << " " << absDev << " " << relDev
              << " -> check " << ((relDev < checkTolerance) ? "passed" : "failed") << std::endl;
    assert(relDev <= checkTolerance);
  }

  {
    std::cout << GridLogMessage << "===========================================================================" << std::endl;
    std::cout << GridLogMessage << "= Test MpcDagMpc is hermitian" << std::endl;
    std::cout << GridLogMessage << "===========================================================================" << std::endl;

    CoarseVector phi(Grid_c); random(pRNG_c, phi);
    CoarseVector chi(Grid_c); random(pRNG_c, chi);
    CoarseVector chi_e(RBGrid_c);   chi_e = Zero();
    CoarseVector chi_o(RBGrid_c);   chi_o = Zero();
    CoarseVector dchi_e(RBGrid_c); dchi_e = Zero();
    CoarseVector dchi_o(RBGrid_c); dchi_o = Zero();
    CoarseVector phi_e(RBGrid_c);   phi_e = Zero();
    CoarseVector phi_o(RBGrid_c);   phi_o = Zero();
    CoarseVector dphi_e(RBGrid_c); dphi_e = Zero();
    CoarseVector dphi_o(RBGrid_c); dphi_o = Zero();

    pickCheckerboard(Even, chi_e, chi);
    pickCheckerboard(Odd,  chi_o, chi);
    pickCheckerboard(Even, phi_e, phi);
    pickCheckerboard(Odd,  phi_o, phi);

    SchurDiagMooeeOperator<CoarseDiracMatrix,CoarseVector> HermOpEO(Dc);

    HermOpEO.MpcDagMpc(chi_e, dchi_e); std::cout << GridLogMessage << "Applied MpcDagMpc to chi_e" << std::endl;
    HermOpEO.MpcDagMpc(chi_o, dchi_o); std::cout << GridLogMessage << "Applied MpcDagMpc to chi_o" << std::endl;
    HermOpEO.MpcDagMpc(phi_e, dphi_e); std::cout << GridLogMessage << "Applied MpcDagMpc to phi_e" << std::endl;
    HermOpEO.MpcDagMpc(phi_o, dphi_o); std::cout << GridLogMessage << "Applied MpcDagMpc to phi_o" << std::endl;

    ComplexD phiDchi_e = innerProduct(phi_e, dchi_e);
    ComplexD phiDchi_o = innerProduct(phi_o, dchi_o);
    ComplexD chiDphi_e = innerProduct(chi_e, dphi_e);
    ComplexD chiDphi_o = innerProduct(chi_o, dphi_o);

    std::cout << GridLogMessage << "e " << phiDchi_e << " " << chiDphi_e << std::endl;
    std::cout << GridLogMessage << "o " << phiDchi_o << " " << chiDphi_o << std::endl;

    std::cout << GridLogMessage << "phiDchi_e - conj(chiDphi_e) " << phiDchi_e - conj(chiDphi_e) << std::endl;
    std::cout << GridLogMessage << "phiDchi_o - conj(chiDphi_o) " << phiDchi_o - conj(chiDphi_o) << std::endl;
  }

  Grid_finalize();
}
