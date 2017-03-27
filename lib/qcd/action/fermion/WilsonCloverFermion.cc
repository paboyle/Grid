/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/action/fermion/WilsonCloverFermion.cc

    Copyright (C) 2017

    Author: paboyle <paboyle@ph.ed.ac.uk>
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
#include <Grid/qcd/spin/Dirac.h>

namespace Grid {
namespace QCD {

  template <class Impl>
  void WilsonCloverFermion<Impl>::AddCloverTerm(const FermionField& in,
                                                FermionField& out){
      FermionField tmp(out._grid);
      tmp = zero;
      // the product sigma_munu Fmunu is hermitian
      tmp += Bx*(Gamma(Gamma::Algebra::SigmaYZ)*in);
      tmp += By*(Gamma(Gamma::Algebra::MinusSigmaXZ)*in);
      tmp += Bz*(Gamma(Gamma::Algebra::SigmaXY)*in);
      tmp += Ex*(Gamma(Gamma::Algebra::MinusSigmaXT)*in);
      tmp += Ey*(Gamma(Gamma::Algebra::MinusSigmaYT)*in);
      tmp += Ez*(Gamma(Gamma::Algebra::MinusSigmaZT)*in);
      out += tmp*csw; // check signs

  }


  template <class Impl>
  RealD WilsonCloverFermion<Impl>::M(const FermionField& in, FermionField& out) {
    // Wilson term
    out.checkerboard = in.checkerboard;
    this->Dhop(in, out, DaggerNo);
    // Clover term
    // apply the sigma and Fmunu
    AddCloverTerm(in, out);
    // overall factor
    return axpy_norm(out, 4 + this->mass, in, out);
  }

  template <class Impl>
  RealD WilsonCloverFermion<Impl>::Mdag(const FermionField& in, FermionField& out) {
    // Wilson term
    out.checkerboard = in.checkerboard;
    this->Dhop(in, out, DaggerYes);
    // Clover term
    // apply the sigma and Fmunu
    AddCloverTerm(in, out);
    return axpy_norm(out, 4 + this->mass, in, out);
  }

  template <class Impl>
  void WilsonCloverFermion<Impl>::ImportGauge(const GaugeField& _Umu) {
    this->ImportGauge(_Umu);
    // Compute the field strength terms
    WilsonLoops<Impl>::FieldStrength(Bx, _Umu, Ydir, Zdir);
    WilsonLoops<Impl>::FieldStrength(By, _Umu, Zdir, Xdir);
    WilsonLoops<Impl>::FieldStrength(Bz, _Umu, Xdir, Ydir);
    WilsonLoops<Impl>::FieldStrength(Ex, _Umu, Tdir, Xdir);
    WilsonLoops<Impl>::FieldStrength(Ey, _Umu, Tdir, Ydir);
    WilsonLoops<Impl>::FieldStrength(Ez, _Umu, Tdir, Zdir);

    // Invert the Moo, Mee terms (?)
  }


  template<class Impl>
  void WilsonCloverFermion<Impl>::Mooee(const FermionField &in, FermionField &out) {
    out.checkerboard = in.checkerboard;
    assert(0); // to be completed
  }

  template<class Impl>
  void WilsonCloverFermion<Impl>::MooeeDag(const FermionField &in, FermionField &out) {
    assert(0); // not implemented yet
  }
  template<class Impl>
  void WilsonCloverFermion<Impl>::MooeeInv(const FermionField &in, FermionField &out) {
    assert(0); // not implemented yet
  }
  template<class Impl>
  void WilsonCloverFermion<Impl>::MooeeInvDag(const FermionField &in, FermionField &out) {
    assert(0); // not implemented yet
  }

  // Derivative parts
  template<class Impl>
  void WilsonCloverFermion<Impl>::MDeriv(GaugeField&mat, const FermionField&U, const FermionField&V, int dag){
    GaugeField tmp(mat._grid);
    this->DhopDeriv(mat, U, V, dag);
    MooDeriv(tmp, U, V, dag);
    mat += tmp;
  }

 // Derivative parts
  template<class Impl>
  void WilsonCloverFermion<Impl>::MooDeriv(GaugeField&mat, const FermionField&U, const FermionField&V, int dag){
    assert(0); // not implemented yet
  }

   // Derivative parts
  template<class Impl>
  void WilsonCloverFermion<Impl>::MeeDeriv(GaugeField&mat, const FermionField&U, const FermionField&V, int dag){
    assert(0); // not implemented yet
  }

  FermOpTemplateInstantiate(WilsonCloverFermion);

}
}
