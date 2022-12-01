/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/utils/SUn.h

Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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
#ifndef QCD_UTIL_SUN_H
#define QCD_UTIL_SUN_H

#define ONLY_IF_SU                                                       \
  typename dummy_name = group_name,                                      \
           typename = std::enable_if_t <                                 \
                          std::is_same<dummy_name, group_name>::value && \
                      is_su<dummy_name>::value >

#define ONLY_IF_Sp                                                       \
  typename dummy_name = group_name,                                      \
           typename = std::enable_if_t <                                 \
                          std::is_same<dummy_name, group_name>::value && \
                      is_sp<dummy_name>::value >

NAMESPACE_BEGIN(Grid);
namespace GroupName {
class SU {};
class Sp {};
}  // namespace GroupName

template <typename group_name>
struct is_su {
  static const bool value = false;
};

template <>
struct is_su<GroupName::SU> {
  static const bool value = true;
};

template <typename group_name>
struct is_sp {
  static const bool value = false;
};

template <>
struct is_sp<GroupName::Sp> {
  static const bool value = true;
};

template <typename group_name>
constexpr int compute_adjoint_dimension(int ncolour);

template <>
constexpr int compute_adjoint_dimension<GroupName::SU>(int ncolour) {
  return ncolour * ncolour - 1;
}

template <>
constexpr int compute_adjoint_dimension<GroupName::Sp>(int ncolour) {
  return ncolour / 2 * (ncolour + 1);
}

template <int ncolour, class group_name>
class GaugeGroup {
 public:
  static const int Dimension = ncolour;
  static const int AdjointDimension =
      compute_adjoint_dimension<group_name>(ncolour);
  static const int AlgebraDimension =
      compute_adjoint_dimension<group_name>(ncolour);

  template <typename vtype>
  using iSU2Matrix = iScalar<iScalar<iMatrix<vtype, 2> > >;
  template <typename vtype>
  using iGroupMatrix = iScalar<iScalar<iMatrix<vtype, ncolour> > >;
  template <typename vtype>
  using iAlgebraVector = iScalar<iScalar<iVector<vtype, AdjointDimension> > >;
  static int su2subgroups(void) { return su2subgroups(group_name()); }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Types can be accessed as SU<2>::Matrix , SU<2>::vSUnMatrix,
  // SU<2>::LatticeMatrix etc...
  //////////////////////////////////////////////////////////////////////////////////////////////////
  typedef iGroupMatrix<Complex> Matrix;
  typedef iGroupMatrix<ComplexF> MatrixF;
  typedef iGroupMatrix<ComplexD> MatrixD;

  typedef iGroupMatrix<vComplex> vMatrix;
  typedef iGroupMatrix<vComplexF> vMatrixF;
  typedef iGroupMatrix<vComplexD> vMatrixD;

  // For the projectors to the algebra
  // these should be real...
  // keeping complex for consistency with the SIMD vector types
  typedef iAlgebraVector<Complex> AlgebraVector;
  typedef iAlgebraVector<ComplexF> AlgebraVectorF;
  typedef iAlgebraVector<ComplexD> AlgebraVectorD;

  typedef iAlgebraVector<vComplex> vAlgebraVector;
  typedef iAlgebraVector<vComplexF> vAlgebraVectorF;
  typedef iAlgebraVector<vComplexD> vAlgebraVectorD;

  typedef Lattice<vMatrix> LatticeMatrix;
  typedef Lattice<vMatrixF> LatticeMatrixF;
  typedef Lattice<vMatrixD> LatticeMatrixD;

  typedef Lattice<vAlgebraVector> LatticeAlgebraVector;
  typedef Lattice<vAlgebraVectorF> LatticeAlgebraVectorF;
  typedef Lattice<vAlgebraVectorD> LatticeAlgebraVectorD;

  typedef iSU2Matrix<Complex> SU2Matrix;
  typedef iSU2Matrix<ComplexF> SU2MatrixF;
  typedef iSU2Matrix<ComplexD> SU2MatrixD;

  typedef iSU2Matrix<vComplex> vSU2Matrix;
  typedef iSU2Matrix<vComplexF> vSU2MatrixF;
  typedef iSU2Matrix<vComplexD> vSU2MatrixD;

  typedef Lattice<vSU2Matrix> LatticeSU2Matrix;
  typedef Lattice<vSU2MatrixF> LatticeSU2MatrixF;
  typedef Lattice<vSU2MatrixD> LatticeSU2MatrixD;

#include "Grid/qcd/utils/SUn.h"
#include "Grid/qcd/utils/Sp2n.h"

 public:
  template <class cplx>
  static void generator(int lieIndex, iGroupMatrix<cplx> &ta) {
    return generator(lieIndex, ta, group_name());
  }

  static void su2SubGroupIndex(int &i1, int &i2, int su2_index) {
    return su2SubGroupIndex(i1, i2, su2_index, group_name());
  }

  static void testGenerators(void) { testGenerators(group_name()); }

  static void printGenerators(void) {
    for (int gen = 0; gen < AdjointDimension; gen++) {
      Matrix ta;
      generator(gen, ta);
      std::cout << GridLogMessage << "Nc = " << ncolour << " t_" << gen
                << std::endl;
      std::cout << GridLogMessage << ta << std::endl;
    }
  }

  // reunitarise??
  template <typename LatticeMatrixType>
  static void LieRandomize(GridParallelRNG &pRNG, LatticeMatrixType &out,
                           double scale = 1.0) {
    GridBase *grid = out.Grid();

    typedef typename LatticeMatrixType::vector_type vector_type;
    typedef typename LatticeMatrixType::scalar_type scalar_type;

    typedef iSinglet<vector_type> vTComplexType;

    typedef Lattice<vTComplexType> LatticeComplexType;
    typedef typename GridTypeMapper<
        typename LatticeMatrixType::vector_object>::scalar_object MatrixType;

    LatticeComplexType ca(grid);
    LatticeMatrixType lie(grid);
    LatticeMatrixType la(grid);
    ComplexD ci(0.0, scale);
    //    ComplexD cone(1.0, 0.0);
    MatrixType ta;

    lie = Zero();

    for (int a = 0; a < AdjointDimension; a++) {
      random(pRNG, ca);

      ca = (ca + conjugate(ca)) * 0.5;
      ca = ca - 0.5;

      generator(a, ta);

      la = ci * ca * ta;

      lie = lie + la;  // e^{i la ta}
    }
    taExp(lie, out);
  }

  static void GaussianFundamentalLieAlgebraMatrix(GridParallelRNG &pRNG,
                                                  LatticeMatrix &out,
                                                  Real scale = 1.0) {
    GridBase *grid = out.Grid();
    LatticeReal ca(grid);
    LatticeMatrix la(grid);
    Complex ci(0.0, scale);
    Matrix ta;

    out = Zero();
    for (int a = 0; a < AdjointDimension; a++) {
      gaussian(pRNG, ca);
      generator(a, ta);
      la = toComplex(ca) * ta;
      out += la;
    }
    out *= ci;
  }

  static void FundamentalLieAlgebraMatrix(const LatticeAlgebraVector &h,
                                          LatticeMatrix &out,
                                          Real scale = 1.0) {
    conformable(h, out);
    GridBase *grid = out.Grid();
    LatticeMatrix la(grid);
    Matrix ta;

    out = Zero();
    for (int a = 0; a < AdjointDimension; a++) {
      generator(a, ta);
      la = peekColour(h, a) * timesI(ta) * scale;
      out += la;
    }
  }

  // Projects the algebra components a lattice matrix (of dimension ncol*ncol -1
  // ) inverse operation: FundamentalLieAlgebraMatrix
  static void projectOnAlgebra(LatticeAlgebraVector &h_out,
                               const LatticeMatrix &in, Real scale = 1.0) {
    conformable(h_out, in);
    h_out = Zero();
    Matrix Ta;

    for (int a = 0; a < AdjointDimension; a++) {
      generator(a, Ta);
      pokeColour(h_out, -2.0 * (trace(timesI(Ta) * in)) * scale, a);
    }
  }

  template <typename GaugeField>
  static void HotConfiguration(GridParallelRNG &pRNG, GaugeField &out) {
    typedef typename GaugeField::vector_type vector_type;
    typedef iGroupMatrix<vector_type> vMatrixType;
    typedef Lattice<vMatrixType> LatticeMatrixType;

    LatticeMatrixType Umu(out.Grid());
    for (int mu = 0; mu < Nd; mu++) {
      LieRandomize(pRNG, Umu, 1.0);
      PokeIndex<LorentzIndex>(out, Umu, mu);
    }
  }
  template <typename GaugeField>
  static void TepidConfiguration(GridParallelRNG &pRNG, GaugeField &out) {
    typedef typename GaugeField::vector_type vector_type;
    typedef iGroupMatrix<vector_type> vMatrixType;
    typedef Lattice<vMatrixType> LatticeMatrixType;

    LatticeMatrixType Umu(out.Grid());
    for (int mu = 0; mu < Nd; mu++) {
      LieRandomize(pRNG, Umu, 0.01);
      PokeIndex<LorentzIndex>(out, Umu, mu);
    }
  }
  template <typename GaugeField>
  static void ColdConfiguration(GaugeField &out) {
    typedef typename GaugeField::vector_type vector_type;
    typedef iGroupMatrix<vector_type> vMatrixType;
    typedef Lattice<vMatrixType> LatticeMatrixType;

    LatticeMatrixType Umu(out.Grid());
    Umu = 1.0;
    for (int mu = 0; mu < Nd; mu++) {
      PokeIndex<LorentzIndex>(out, Umu, mu);
    }
  }
  template <typename GaugeField>
  static void ColdConfiguration(GridParallelRNG &pRNG, GaugeField &out) {
    ColdConfiguration(out);
  }

  template <typename LatticeMatrixType, ONLY_IF_SU>
  static void taProj(const LatticeMatrixType &in, LatticeMatrixType &out) {
    out = Ta(in);
  }
  template <typename LatticeMatrixType>
  static void taExp(const LatticeMatrixType &x, LatticeMatrixType &ex) {
    typedef typename LatticeMatrixType::scalar_type ComplexType;

    LatticeMatrixType xn(x.Grid());
    RealD nfac = 1.0;

    xn = x;
    ex = xn + ComplexType(1.0);  // 1+x

    // Do a 12th order exponentiation
    for (int i = 2; i <= 12; ++i) {
      nfac = nfac / RealD(i);  // 1/2, 1/2.3 ...
      xn = xn * x;             // x2, x3,x4....
      ex = ex + xn * nfac;     // x2/2!, x3/3!....
    }
  }
};

template <int N>
LatticeComplexD Determinant(
    const Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > &Umu) {
  GridBase *grid = Umu.Grid();
  auto lvol = grid->lSites();
  LatticeComplexD ret(grid);

  autoView(Umu_v, Umu, CpuRead);
  autoView(ret_v, ret, CpuWrite);
  thread_for(site, lvol, {
    Eigen::MatrixXcd EigenU = Eigen::MatrixXcd::Zero(N, N);
    Coordinate lcoor;
    grid->LocalIndexToLocalCoor(site, lcoor);
    iScalar<iScalar<iMatrix<ComplexD, N> > > Us;
    peekLocalSite(Us, Umu_v, lcoor);
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        EigenU(i, j) = Us()()(i, j);
      }
    }
    ComplexD det = EigenU.determinant();
    pokeLocalSite(det, ret_v, lcoor);
  });
  return ret;
}
template <int N>
static void ProjectSUn(
    Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > &Umu) {
  Umu = ProjectOnGroup(Umu);
  auto det = Determinant(Umu);

  det = conjugate(det);

  for (int i = 0; i < N; i++) {
    auto element = PeekIndex<ColourIndex>(Umu, N - 1, i);
    element = element * det;
    PokeIndex<ColourIndex>(Umu, element, Nc - 1, i);
  }
}
template <int N>
static void ProjectSUn(
    Lattice<iVector<iScalar<iMatrix<vComplexD, N> >, Nd> > &U) {
  GridBase *grid = U.Grid();
  // Reunitarise
  for (int mu = 0; mu < Nd; mu++) {
    auto Umu = PeekIndex<LorentzIndex>(U, mu);
    Umu = ProjectOnGroup(Umu);
    ProjectSUn(Umu);
    PokeIndex<LorentzIndex>(U, Umu, mu);
  }
}
// Explicit specialisation for SU(3).
// Explicit specialisation for SU(3).
static void ProjectSU3(
    Lattice<iScalar<iScalar<iMatrix<vComplexD, 3> > > > &Umu) {
  GridBase *grid = Umu.Grid();
  const int x = 0;
  const int y = 1;
  const int z = 2;
  // Reunitarise
  Umu = ProjectOnGroup(Umu);
  autoView(Umu_v, Umu, CpuWrite);
  thread_for(ss, grid->oSites(), {
    auto cm = Umu_v[ss];
    cm()()(2, x) = adj(cm()()(0, y) * cm()()(1, z) -
                       cm()()(0, z) * cm()()(1, y));  // x= yz-zy
    cm()()(2, y) = adj(cm()()(0, z) * cm()()(1, x) -
                       cm()()(0, x) * cm()()(1, z));  // y= zx-xz
    cm()()(2, z) = adj(cm()()(0, x) * cm()()(1, y) -
                       cm()()(0, y) * cm()()(1, x));  // z= xy-yx
    Umu_v[ss] = cm;
  });
}
static void ProjectSU3(
    Lattice<iVector<iScalar<iMatrix<vComplexD, 3> >, Nd> > &U) {
  GridBase *grid = U.Grid();
  // Reunitarise
  for (int mu = 0; mu < Nd; mu++) {
    auto Umu = PeekIndex<LorentzIndex>(U, mu);
    Umu = ProjectOnGroup(Umu);
    ProjectSU3(Umu);
    PokeIndex<LorentzIndex>(U, Umu, mu);
  }
}

template <int ncolour>
using SU = GaugeGroup<ncolour, GroupName::SU>;

template <int ncolour>
using Sp = GaugeGroup<ncolour, GroupName::Sp>;

typedef SU<2> SU2;
typedef SU<3> SU3;
typedef SU<4> SU4;
typedef SU<5> SU5;

typedef SU<Nc> FundamentalMatrices;

template <int N>
static void ProjectSp2n(
    Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > &Umu) {
  Umu = ProjectOnSpGroup(Umu);
  auto det = Determinant(Umu);  // ok ?

  det = conjugate(det);

  for (int i = 0; i < N; i++) {
    auto element = PeekIndex<ColourIndex>(Umu, N - 1, i);
    element = element * det;
    PokeIndex<ColourIndex>(Umu, element, Nc - 1, i);
  }
}
template <int N>
static void ProjectSp2n(
    Lattice<iVector<iScalar<iMatrix<vComplexD, N> >, Nd> > &U) {
  GridBase *grid = U.Grid();
  for (int mu = 0; mu < Nd; mu++) {
    auto Umu = PeekIndex<LorentzIndex>(U, mu);
    Umu = ProjectOnSpGroup(Umu);
    ProjectSp2n(Umu);
    PokeIndex<LorentzIndex>(U, Umu, mu);
  }
}

typedef Sp<2> Sp2;
typedef Sp<4> Sp4;
typedef Sp<6> Sp6;
typedef Sp<8> Sp8;


NAMESPACE_END(Grid);
#endif
