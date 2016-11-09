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

#ifndef QCD_UTIL_SUN2INDEX_H
#define QCD_UTIL_SUN2INDEX_H


namespace Grid {
namespace QCD {

enum TwoIndexSymmetry { Symmetric = 1, AntiSymmetric = -1 };

inline Real delta(int a, int b) { return (a == b) ? 1.0 : 0.0; }

template <int ncolour, TwoIndexSymmetry S>
class SU_TwoIndex : public SU<ncolour> {
 public:
  static const int Dimension = ncolour * (ncolour + S) / 2;
  static const int NumGenerators = SU<ncolour>::AdjointDimension;

  template <typename vtype>
  using iSUnTwoIndexMatrix = iScalar<iScalar<iMatrix<vtype, Dimension> > >;

  typedef iSUnTwoIndexMatrix<Complex> TIMatrix;
  typedef iSUnTwoIndexMatrix<ComplexF> TIMatrixF;
  typedef iSUnTwoIndexMatrix<ComplexD> TIMatrixD;

  typedef iSUnTwoIndexMatrix<vComplex> vTIMatrix;
  typedef iSUnTwoIndexMatrix<vComplexF> vTIMatrixF;
  typedef iSUnTwoIndexMatrix<vComplexD> vTIMatrixD;

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
  using iSUnMatrix = iScalar<iScalar<iMatrix<vtype, ncolour> > >;

  typedef iSUnMatrix<Complex> Matrix;
  typedef iSUnMatrix<ComplexF> MatrixF;
  typedef iSUnMatrix<ComplexD> MatrixD;

  template <class cplx>
  static void base(int Index, iSUnMatrix<cplx> &eij) {
    // returns (e)^(ij)_{kl} necessary for change of base U_F -> U_R
    assert(Index < NumGenerators);
    eij = zero;

    // for the linearisation of the 2 indexes 
    static int a[ncolour * (ncolour - 1) / 2][2]; // store the a <-> i,j
    static bool filled = false;
    if (!filled) {
      int counter = 0;
      for (int i = 1; i < ncolour; i++) {
        for (int j = 0; j < i; j++) {
          a[counter][0] = i;
          a[counter][1] = j;
          counter++;
        }
      }
      filled = true;
    }

    if (Index < ncolour * (ncolour - 1) / 2) {
      baseOffDiagonal(a[Index][0], a[Index][1], eij);
    } else {
      baseDiagonal(Index, eij);
    }
  }

  template <class cplx>
  static void baseDiagonal(int Index, iSUnMatrix<cplx> &eij) {
    eij = zero;
    eij()()(Index - ncolour * (ncolour - 1) / 2,
            Index - ncolour * (ncolour - 1) / 2) = 1.0;
  }

  template <class cplx>
  static void baseOffDiagonal(int i, int j, iSUnMatrix<cplx> &eij) {
    eij = zero;
    for (int k = 0; k < ncolour; k++)
      for (int l = 0; l < ncolour; l++)
        eij()()(l, k) = delta(i, k) * delta(j, l) +
                        S * delta(j, k) * delta(i, l);

    RealD nrm = 1. / std::sqrt(2.0);
    eij = eij * nrm;
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
  static void generator(int Index, iSUnTwoIndexMatrix<cplx> &i2indTa) {
    Vector<typename SU<ncolour>::template iSUnMatrix<cplx> > ta(
        ncolour * ncolour - 1);
    Vector<typename SU<ncolour>::template iSUnMatrix<cplx> > eij(Dimension);
    typename SU<ncolour>::template iSUnMatrix<cplx> tmp;
    i2indTa = zero;
    
    for (int a = 0; a < ncolour * ncolour - 1; a++)
      SU<ncolour>::generator(a, ta[a]);
    
    for (int a = 0; a < Dimension; a++) base(a, eij[a]);

    for (int a = 0; a < Dimension; a++) {
      tmp = transpose(ta[Index]) * adj(eij[a]) + adj(eij[a]) * ta[Index];
      for (int b = 0; b < Dimension; b++) {
        typename SU<ncolour>::template iSUnMatrix<cplx> tmp1 =
            tmp * eij[b]; 
        Complex iTr = TensorRemove(timesI(trace(tmp1)));
        i2indTa()()(a, b) = iTr;
      }
    }
  }

  static void printGenerators(void) {
    for (int gen = 0; gen < ncolour * ncolour - 1; gen++) {
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
    for (int a = 0; a < ncolour * ncolour - 1; a++) {
      generator(a, i2indTa);
      std::cout << GridLogMessage << a << std::endl;
      assert(norm2(trace(i2indTa)) < 1.0e-6);
    }
    std::cout << GridLogMessage << std::endl;

    std::cout << GridLogMessage << "2IndexRep - Checking if antihermitean"
              << std::endl;
    for (int a = 0; a < ncolour * ncolour - 1; a++) {
      generator(a, i2indTa);
      std::cout << GridLogMessage << a << std::endl;
      assert(norm2(adj(i2indTa) + i2indTa) < 1.0e-6);
    }

    std::cout << GridLogMessage << std::endl;
    std::cout << GridLogMessage
              << "2IndexRep - Checking Tr[Ta*Tb]=delta(a,b)*(N +- 2)/2"
              << std::endl;
    for (int a = 0; a < ncolour * ncolour - 1; a++) {
      for (int b = 0; b < ncolour * ncolour - 1; b++) {
        generator(a, i2indTa);
        generator(b, i2indTb);

        // generator returns iTa, so we need a minus sign here
        Complex Tr = -TensorRemove(trace(i2indTa * i2indTb));
        std::cout << GridLogMessage << "a=" << a << "b=" << b << "Tr=" << Tr
                  << std::endl;
      }
    }
    std::cout << GridLogMessage << std::endl;
  }

  static void TwoIndexLieAlgebraMatrix(
      const typename SU<ncolour>::LatticeAlgebraVector &h,
      LatticeTwoIndexMatrix &out, Real scale = 1.0) {
    conformable(h, out);
    GridBase *grid = out._grid;
    LatticeTwoIndexMatrix la(grid);
    TIMatrix i2indTa;

    out = zero;
    for (int a = 0; a < ncolour * ncolour - 1; a++) {
      generator(a, i2indTa);
      la = peekColour(h, a) * i2indTa;
      out += la;
    }
    out *= scale;
  }

  // Projects the algebra components 
  // of a lattice matrix ( of dimension ncol*ncol -1 )
  static void projectOnAlgebra(
      typename SU<ncolour>::LatticeAlgebraVector &h_out,
      const LatticeTwoIndexMatrix &in, Real scale = 1.0) {
    conformable(h_out, in);
    h_out = zero;
    TIMatrix i2indTa;
    Real coefficient = -2.0 / (ncolour + 2 * S) * scale;
    // 2/(Nc +/- 2) for the normalization of the trace in the two index rep
    for (int a = 0; a < ncolour * ncolour - 1; a++) {
      generator(a, i2indTa);
      auto tmp = real(trace(i2indTa * in)) * coefficient;
      pokeColour(h_out, tmp, a);
    }
  }

  // a projector that keeps the generators stored to avoid the overhead of
  // recomputing them
  static void projector(typename SU<ncolour>::LatticeAlgebraVector &h_out,
                        const LatticeTwoIndexMatrix &in, Real scale = 1.0) {
    conformable(h_out, in);
    // to store the generators
    static std::vector<TIMatrix> i2indTa(ncolour * ncolour -1); 
    h_out = zero;
    static bool precalculated = false;
    if (!precalculated) {
      precalculated = true;
      for (int a = 0; a < ncolour * ncolour - 1; a++) generator(a, i2indTa[a]);
    }

    Real coefficient =
        -2.0 / (ncolour + 2 * S) * scale;  // 2/(Nc +/- 2) for the normalization
                                           // of the trace in the two index rep

    for (int a = 0; a < ncolour * ncolour - 1; a++) {
      auto tmp = real(trace(i2indTa[a] * in)) * coefficient;
      pokeColour(h_out, tmp, a);
    }
  }
};

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


}
}

#endif
