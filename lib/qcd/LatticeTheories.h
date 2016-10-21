/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/QCD.h

Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: neo <cossu@post.kek.jp>
Author: paboyle <paboyle@ph.ed.ac.uk>

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#ifndef GRID_LT_H
#define GRID_LT_H
namespace Grid{

// First steps in the complete generalization of the Physics part
// Design not final
namespace LatticeTheories {

template <int Dimensions>
struct LatticeTheory {
  static const int Nd = Dimensions;
  static const int Nds = Dimensions * 2;  // double stored field
  template <typename vtype>
  using iSinglet = iScalar<iScalar<iScalar<vtype> > >;
};

template <int Dimensions, int Colours>
struct LatticeGaugeTheory : public LatticeTheory<Dimensions> {
  static const int Nds = Dimensions * 2;
  static const int Nd = Dimensions;
  static const int Nc = Colours;

  template <typename vtype> 
  using iColourMatrix = iScalar<iScalar<iMatrix<vtype, Nc> > >;
  template <typename vtype>
  using iLorentzColourMatrix = iVector<iScalar<iMatrix<vtype, Nc> >, Nd>;
  template <typename vtype>
  using iDoubleStoredColourMatrix = iVector<iScalar<iMatrix<vtype, Nc> >, Nds>;
  template <typename vtype>
  using iColourVector = iScalar<iScalar<iVector<vtype, Nc> > >;
};

template <int Dimensions, int Colours, int Spin>
struct FermionicLatticeGaugeTheory
    : public LatticeGaugeTheory<Dimensions, Colours> {
  static const int Nd = Dimensions;
  static const int Nds = Dimensions * 2;
  static const int Nc = Colours;
  static const int Ns = Spin;

  template <typename vtype>
  using iSpinMatrix = iScalar<iMatrix<iScalar<vtype>, Ns> >;
  template <typename vtype>
  using iSpinColourMatrix = iScalar<iMatrix<iMatrix<vtype, Nc>, Ns> >;
  template <typename vtype>
  using iSpinVector = iScalar<iVector<iScalar<vtype>, Ns> >;
  template <typename vtype>
  using iSpinColourVector = iScalar<iVector<iVector<vtype, Nc>, Ns> >;
  // These 2 only if Spin is a multiple of 2
  static const int Nhs = Spin / 2;
  template <typename vtype>
  using iHalfSpinVector = iScalar<iVector<iScalar<vtype>, Nhs> >;
  template <typename vtype>
  using iHalfSpinColourVector = iScalar<iVector<iVector<vtype, Nc>, Nhs> >;

  //tests
  typedef iColourMatrix<Complex> ColourMatrix;
  typedef iColourMatrix<ComplexF> ColourMatrixF;
  typedef iColourMatrix<ComplexD> ColourMatrixD;


};

// Examples, not complete now.
struct QCD : public FermionicLatticeGaugeTheory<4, 3, 4> {
    static const int Xp = 0;
    static const int Yp = 1;
    static const int Zp = 2;
    static const int Tp = 3;
    static const int Xm = 4;
    static const int Ym = 5;
    static const int Zm = 6;
    static const int Tm = 7;

    typedef FermionicLatticeGaugeTheory FLGT;

    typedef FLGT::iSpinMatrix<Complex  >          SpinMatrix;
    typedef FLGT::iSpinMatrix<ComplexF >          SpinMatrixF;
    typedef FLGT::iSpinMatrix<ComplexD >          SpinMatrixD;

};
struct QED : public FermionicLatticeGaugeTheory<4, 1, 4> {//fill
};

template <int Dimensions>
struct Scalar : public LatticeTheory<Dimensions> {};

};  // LatticeTheories

} // Grid

#endif
