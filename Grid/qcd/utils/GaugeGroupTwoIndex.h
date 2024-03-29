////////////////////////////////////////////////////////////////////////
//
// * Two index representation generators
//
// * Normalisation for the fundamental generators:
//   trace ta tb = 1/2 delta_ab = T_F delta_ab
//   T_F = 1/2  for SU(N) groups
//
//
//   base for NxN two index (anti-symmetric) matrices
//   normalized to 1 (d_ij is the kroenecker delta)
//
//   (e^(ij)_{kl} = 1 / sqrt(2) (d_ik d_jl +/- d_jk d_il)
//
//   Then the generators are written as
//
//   (iT_a)^(ij)(lk) = i * ( tr[e^(ij)^dag e^(lk) T^trasp_a] +
//   tr[e^(lk)e^(ij)^dag T_a] )  //
//
//
////////////////////////////////////////////////////////////////////////

// Authors: David Preti, Guido Cossu

#ifndef QCD_UTIL_GAUGEGROUPTWOINDEX_H
#define QCD_UTIL_GAUGEGROUPTWOINDEX_H

NAMESPACE_BEGIN(Grid);

enum TwoIndexSymmetry { Symmetric = 1, AntiSymmetric = -1 };

constexpr inline Real delta(int a, int b) { return (a == b) ? 1.0 : 0.0; }

namespace detail {

template <class cplx, int nc, TwoIndexSymmetry S>
struct baseOffDiagonalSpHelper;

template <class cplx, int nc>
struct baseOffDiagonalSpHelper<cplx, nc, AntiSymmetric> {
  static const int ngroup = nc / 2;
  static void baseOffDiagonalSp(int i, int j, iScalar<iScalar<iMatrix<cplx, nc> > > &eij) {
    eij = Zero();
    RealD tmp;

    if ((i == ngroup + j) && (1 <= j) && (j < ngroup)) {
      for (int k = 0; k < j+1; k++) {
        if (k < j) {
          tmp = 1 / sqrt(j * (j + 1));
          eij()()(k, k + ngroup) = tmp;
          eij()()(k + ngroup, k) = -tmp;
        }
        if (k == j) {
          tmp = -j / sqrt(j * (j + 1));
          eij()()(k, k + ngroup) = tmp;
          eij()()(k + ngroup, k) = -tmp;
        }
      }

    }

    else if (i != ngroup + j) {
      for (int k = 0; k < nc; k++)
        for (int l = 0; l < nc; l++) {
          eij()()(l, k) =
              delta(i, k) * delta(j, l) - delta(j, k) * delta(i, l);
        }
    }
    RealD nrm = 1. / std::sqrt(2.0);
    eij = eij * nrm;
  }
};

template <class cplx, int nc>
struct baseOffDiagonalSpHelper<cplx, nc, Symmetric> {
  static void baseOffDiagonalSp(int i, int j, iScalar<iScalar<iMatrix<cplx, nc> > > &eij) {
    eij = Zero();
    for (int k = 0; k < nc; k++)
      for (int l = 0; l < nc; l++)
        eij()()(l, k) =
            delta(i, k) * delta(j, l) + delta(j, k) * delta(i, l);

    RealD nrm = 1. / std::sqrt(2.0);
    eij = eij * nrm;
  }
};

}   // closing detail namespace

template <int ncolour, TwoIndexSymmetry S, class group_name>
class GaugeGroupTwoIndex : public GaugeGroup<ncolour, group_name> {
 public:
  // The chosen convention is that we are taking ncolour to be N in SU<N> but 2N
  // in Sp(2N). ngroup is equal to N for SU but 2N/2 = N for Sp(2N).
  static_assert(std::is_same<group_name, GroupName::SU>::value or
                    std::is_same<group_name, GroupName::Sp>::value,
                "ngroup is only implemented for SU and Sp currently.");
  static const int ngroup =
      std::is_same<group_name, GroupName::SU>::value ? ncolour : ncolour / 2;
  static const int Dimension =
      (ncolour * (ncolour + S) / 2) + (std::is_same<group_name, GroupName::Sp>::value ? (S - 1) / 2 : 0);
  static const int DimensionAS =
      (ncolour * (ncolour - 1) / 2) + (std::is_same<group_name, GroupName::Sp>::value ? (- 1) : 0);
  static const int DimensionS =
      ncolour * (ncolour + 1) / 2;
  static const int NumGenerators =
      GaugeGroup<ncolour, group_name>::AlgebraDimension;

  template <typename vtype>
  using iGroupTwoIndexMatrix = iScalar<iScalar<iMatrix<vtype, Dimension> > >;

  typedef iGroupTwoIndexMatrix<Complex> TIMatrix;
  typedef iGroupTwoIndexMatrix<ComplexF> TIMatrixF;
  typedef iGroupTwoIndexMatrix<ComplexD> TIMatrixD;

  typedef iGroupTwoIndexMatrix<vComplex> vTIMatrix;
  typedef iGroupTwoIndexMatrix<vComplexF> vTIMatrixF;
  typedef iGroupTwoIndexMatrix<vComplexD> vTIMatrixD;

  typedef Lattice<vTIMatrix> LatticeTwoIndexMatrix;
  typedef Lattice<vTIMatrixF> LatticeTwoIndexMatrixF;
  typedef Lattice<vTIMatrixD> LatticeTwoIndexMatrixD;

  typedef Lattice<iVector<iScalar<iMatrix<vComplex, Dimension> >, Nd> >
      LatticeTwoIndexField;
  typedef Lattice<iVector<iScalar<iMatrix<vComplexF, Dimension> >, Nd> >
      LatticeTwoIndexFieldF;
  typedef Lattice<iVector<iScalar<iMatrix<vComplexD, Dimension> >, Nd> >
      LatticeTwoIndexFieldD;

  template <typename vtype>
  using iGroupMatrix = iScalar<iScalar<iMatrix<vtype, ncolour> > >;

  typedef iGroupMatrix<Complex> Matrix;
  typedef iGroupMatrix<ComplexF> MatrixF;
  typedef iGroupMatrix<ComplexD> MatrixD;
    
private:
  template <class cplx>
  static void baseDiagonal(int Index, iGroupMatrix<cplx> &eij) {
    eij = Zero();
    eij()()(Index - ncolour * (ncolour - 1) / 2,
            Index - ncolour * (ncolour - 1) / 2) = 1.0;
  }
    
  template <class cplx>
  static void baseOffDiagonal(int i, int j, iGroupMatrix<cplx> &eij, GroupName::SU) {
    eij = Zero();
    for (int k = 0; k < ncolour; k++)
      for (int l = 0; l < ncolour; l++)
        eij()()(l, k) =
            delta(i, k) * delta(j, l) + S * delta(j, k) * delta(i, l);

    RealD nrm = 1. / std::sqrt(2.0);
    eij = eij * nrm;
  }
    
  template <class cplx>
  static void baseOffDiagonal(int i, int j, iGroupMatrix<cplx> &eij, GroupName::Sp) {
    detail::baseOffDiagonalSpHelper<cplx, ncolour, S>::baseOffDiagonalSp(i, j, eij);
  }

public:
    
  template <class cplx>
  static void base(int Index, iGroupMatrix<cplx> &eij) {
  // returns (e)^(ij)_{kl} necessary for change of base U_F -> U_R
    assert(Index < Dimension);
    eij = Zero();
  // for the linearisation of the 2 indexes
    static int a[ncolour * (ncolour - 1) / 2][2];  // store the a <-> i,j
    static bool filled = false;
    if (!filled) {
      int counter = 0;
      for (int i = 1; i < ncolour; i++) {
      for (int j = 0; j < i; j++) {
        if (std::is_same<group_name, GroupName::Sp>::value)
          {
            if (j==0 && i==ngroup+j && S==-1) {
            //std::cout << "skipping" << std::endl; // for Sp2n this vanishes identically.
              j = j+1;
            }
          }
          a[counter][0] = i;
          a[counter][1] = j;
          counter++;
          }
      }
      filled = true;
    }
    if (Index < ncolour*ncolour - DimensionS)
    {
      baseOffDiagonal(a[Index][0], a[Index][1], eij, group_name());
    } else {
      baseDiagonal(Index, eij);
    }
  }
    
  static void printBase(void) {
    for (int gen = 0; gen < Dimension; gen++) {
      Matrix tmp;
      base(gen, tmp);
      std::cout << GridLogMessage << "Nc = " << ncolour << " t_" << gen
                << std::endl;
      std::cout << GridLogMessage << tmp << std::endl;
    }
  }

  template <class cplx>
  static void generator(int Index, iGroupTwoIndexMatrix<cplx> &i2indTa) {
    Vector<iGroupMatrix<cplx> > ta(NumGenerators);
    Vector<iGroupMatrix<cplx> > eij(Dimension);
    iGroupMatrix<cplx> tmp;

    for (int a = 0; a < NumGenerators; a++)
      GaugeGroup<ncolour, group_name>::generator(a, ta[a]);

    for (int a = 0; a < Dimension; a++) base(a, eij[a]);

    for (int a = 0; a < Dimension; a++) {
      tmp = transpose(eij[a]*ta[Index]) + transpose(eij[a]) * ta[Index];
      for (int b = 0; b < Dimension; b++) {
        Complex iTr = TensorRemove(timesI(trace(tmp * eij[b])));
        i2indTa()()(a, b) = iTr;
      }
    }
  }

  static void printGenerators(void) {
    for (int gen = 0; gen < NumGenerators; gen++) {
      TIMatrix i2indTa;
      generator(gen, i2indTa);
      std::cout << GridLogMessage << "Nc = " << ncolour << " t_" << gen
                << std::endl;
      std::cout << GridLogMessage << i2indTa << std::endl;
    }
  }

  static void testGenerators(void) {
    TIMatrix i2indTa, i2indTb;
    std::cout << GridLogMessage << "2IndexRep - Checking if traceless"
              << std::endl;
    for (int a = 0; a < NumGenerators; a++) {
      generator(a, i2indTa);
      std::cout << GridLogMessage << a << std::endl;
      assert(norm2(trace(i2indTa)) < 1.0e-6);
    }
    std::cout << GridLogMessage << std::endl;

    std::cout << GridLogMessage << "2IndexRep - Checking if antihermitean"
              << std::endl;
    for (int a = 0; a < NumGenerators; a++) {
      generator(a, i2indTa);
      std::cout << GridLogMessage << a << std::endl;
      assert(norm2(adj(i2indTa) + i2indTa) < 1.0e-6);
    }

    std::cout << GridLogMessage << std::endl;
    std::cout << GridLogMessage
              << "2IndexRep - Checking Tr[Ta*Tb]=delta(a,b)*(N +- 2)/2"
              << std::endl;
    for (int a = 0; a < NumGenerators; a++) {
      for (int b = 0; b < NumGenerators; b++) {
        generator(a, i2indTa);
        generator(b, i2indTb);

        // generator returns iTa, so we need a minus sign here
        Complex Tr = -TensorRemove(trace(i2indTa * i2indTb));
        std::cout << GridLogMessage << "a=" << a << "b=" << b << "Tr=" << Tr
                  << std::endl;
        if (a == b) {
          assert(real(Tr) - ((ncolour + S * 2) * 0.5) < 1e-8);
        } else {
          assert(real(Tr) < 1e-8);
        }
        assert(imag(Tr) < 1e-8);
      }
    }
    std::cout << GridLogMessage << std::endl;
  }

  static void TwoIndexLieAlgebraMatrix(
      const typename GaugeGroup<ncolour, group_name>::LatticeAlgebraVector &h,
      LatticeTwoIndexMatrix &out, Real scale = 1.0) {
    conformable(h, out);
    GridBase *grid = out.Grid();
    LatticeTwoIndexMatrix la(grid);
    TIMatrix i2indTa;

    out = Zero();
    for (int a = 0; a < NumGenerators; a++) {
      generator(a, i2indTa);
      la = peekColour(h, a) * i2indTa;
      out += la;
    }
    out *= scale;
  }

  // Projects the algebra components
  // of a lattice matrix ( of dimension ncol*ncol -1 )
  static void projectOnAlgebra(
      typename GaugeGroup<ncolour, group_name>::LatticeAlgebraVector &h_out,
      const LatticeTwoIndexMatrix &in, Real scale = 1.0) {
    conformable(h_out, in);
    h_out = Zero();
    TIMatrix i2indTa;
    Real coefficient = -2.0 / (ncolour + 2 * S) * scale;
    // 2/(Nc +/- 2) for the normalization of the trace in the two index rep
    for (int a = 0; a < NumGenerators; a++) {
      generator(a, i2indTa);
      pokeColour(h_out, real(trace(i2indTa * in)) * coefficient, a);
    }
  }

  // a projector that keeps the generators stored to avoid the overhead of
  // recomputing them
  static void projector(
      typename GaugeGroup<ncolour, group_name>::LatticeAlgebraVector &h_out,
      const LatticeTwoIndexMatrix &in, Real scale = 1.0) {
    conformable(h_out, in);
    // to store the generators
    static std::vector<TIMatrix> i2indTa(NumGenerators);
    h_out = Zero();
    static bool precalculated = false;
    if (!precalculated) {
      precalculated = true;
      for (int a = 0; a < NumGenerators; a++) generator(a, i2indTa[a]);
    }

    Real coefficient =
        -2.0 / (ncolour + 2 * S) * scale;  // 2/(Nc +/- 2) for the normalization
    // of the trace in the two index rep

    for (int a = 0; a < NumGenerators; a++) {
      auto tmp = real(trace(i2indTa[a] * in)) * coefficient;
      pokeColour(h_out, tmp, a);
    }
  }
};

template <int ncolour, TwoIndexSymmetry S>
using SU_TwoIndex = GaugeGroupTwoIndex<ncolour, S, GroupName::SU>;

// Some useful type names
typedef SU_TwoIndex<Nc, Symmetric> TwoIndexSymmMatrices;
typedef SU_TwoIndex<Nc, AntiSymmetric> TwoIndexAntiSymmMatrices;

typedef SU_TwoIndex<2, Symmetric> SU2TwoIndexSymm;
typedef SU_TwoIndex<3, Symmetric> SU3TwoIndexSymm;
typedef SU_TwoIndex<4, Symmetric> SU4TwoIndexSymm;
typedef SU_TwoIndex<5, Symmetric> SU5TwoIndexSymm;

typedef SU_TwoIndex<2, AntiSymmetric> SU2TwoIndexAntiSymm;
typedef SU_TwoIndex<3, AntiSymmetric> SU3TwoIndexAntiSymm;
typedef SU_TwoIndex<4, AntiSymmetric> SU4TwoIndexAntiSymm;
typedef SU_TwoIndex<5, AntiSymmetric> SU5TwoIndexAntiSymm;

template <int ncolour, TwoIndexSymmetry S>
using Sp_TwoIndex = GaugeGroupTwoIndex<ncolour, S, GroupName::Sp>;

typedef Sp_TwoIndex<Nc, Symmetric> SpTwoIndexSymmMatrices;
typedef Sp_TwoIndex<Nc, AntiSymmetric> SpTwoIndexAntiSymmMatrices;

typedef Sp_TwoIndex<2, Symmetric> Sp2TwoIndexSymm;
typedef Sp_TwoIndex<4, Symmetric> Sp4TwoIndexSymm;

typedef Sp_TwoIndex<4, AntiSymmetric> Sp4TwoIndexAntiSymm;

NAMESPACE_END(Grid);

#endif
