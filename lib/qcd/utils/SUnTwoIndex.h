#ifndef QCD_UTIL_SUNADJOINT_H
#define QCD_UTIL_SUNADJOINT_H

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
//   (iT^(ij))_lk = i 
//
////////////////////////////////////////////////////////////////////////

namespace Grid {
  namespace QCD {

    enum TwoIndexSymmetry {Symmetric = 1, AntiSymmetric = -1};
    
    template <int ncolour, TwoIndexSymmetry S>
    class SU_TwoIndex : public SU<ncolour> {
    public:
      static const int Dimension = ncolour * (ncolour + S) / 2;
    
      template <typename vtype>
      using iSUnTwoIndexMatrix =
	iScalar<iScalar<iMatrix<vtype, Dimension > > >;
    
      typedef iSUnTwoIndexMatrix<Complex> TIMatrix;
      typedef iSUnTwoIndexMatrix<ComplexF> TIMatrixF;
      typedef iSUnTwoIndexMatrix<ComplexD> TIMatrixD;
    
      typedef iSUnTwoIndexMatrix<vComplex> vTIMatrix;
      typedef iSUnTwoIndexMatrix<vComplexF> vTIMatrixF;
      typedef iSUnTwoIndexMatrix<vComplexD> vTIMatrixD;
    
      typedef Lattice<vAMatrix>  LatticeTwoIndexMatrix;
      typedef Lattice<vAMatrixF> LatticeTwoIndexMatrixF;
      typedef Lattice<vAMatrixD> LatticeTwoIndexMatrixD;
    
      typedef Lattice<iVector<iScalar<iMatrix<vComplex, Dimension> >, Nd> >
      LatticeTwoIndexField;
      typedef Lattice<iVector<iScalar<iMatrix<vComplexF, Dimension> >, Nd> >
      LatticeTwoIndexFieldF;
      typedef Lattice<iVector<iScalar<iMatrix<vComplexD, Dimension> >, Nd> >
      LatticeTwoIndexFieldD;
    
    


      template <class cplx>
      static void generator(int Index, iSUnTwoIndexMatrix<cplx> &iTwoIdxTa) {
	// returns i(T)^(ij) necessary for the projectors
	// see definitions above
	iTwoIdxTa = zero;
	Vector<typename SU<ncolour>::template iSUnMatrix<cplx> > tij(Dimension);
	typename SU<ncolour>::template iSUnMatrix<cplx> tmp;


	for (int a = 0; a < Dimension; a++) {

	}
      }

      static void printGenerators(void) {
	for (int gen = 0; gen < Dimension; gen++) {
	  AMatrix ta;
	  generator(gen, ta);
	  std::cout << GridLogMessage << "Nc = " << ncolour << " t_" << gen
		    << std::endl;
	  std::cout << GridLogMessage << ta << std::endl;
	}
      }

      static void testGenerators(void) {
	TIMatrix TwoIndexTa;

      }

      static void TwoIndexLieAlgebraMatrix(const typename SU<ncolour>::LatticeAlgebraVector &h,
					   LatticeTwoIndexMatrix &out, Real scale = 1.0) {
	conformable(h, out);
	GridBase *grid = out._grid;
	LatticeAdjMatrix la(grid);
	TIMatrix iTa;

	out = zero;
	for (int a = 0; a < Dimension; a++) {
	  generator(a, iTa);
	  la = peekColour(h, a) * iTa;
	  out += la;
	}
	out *= scale;
      }

      // Projects the algebra components a lattice matrix (of dimension ncol*ncol -1 )
      static void projectOnAlgebra(typename SU<ncolour>::LatticeAlgebraVector &h_out, const LatticeTwoIndexMatrix &in, Real scale = 1.0) {
	conformable(h_out, in);
	h_out = zero;
	TIMatrix iTa;
	Real coefficient = - 2.0/(ncolour + 2*S) * scale;// 2/(Nc +/- 2) for the normalization of the trace in the two index rep

	for (int a = 0; a < Dimension; a++) {
	  generator(a, iTa);
	  auto tmp = real(trace(iTa * in)) * coefficient;
	  pokeColour(h_out, tmp, a);
	}
      }

      // a projector that keeps the generators stored to avoid the overhead of recomputing them 
      static void projector(typename SU<ncolour>::LatticeAlgebraVector &h_out, const LatticeTwoIndexMatrix &in, Real scale = 1.0) {
	conformable(h_out, in);
	static std::vector<TIMatrix> iTa(Dimension);  // to store the generators
	h_out = zero;
	static bool precalculated = false; 
	if (!precalculated){
	  precalculated = true;
	  for (int a = 0; a < Dimension; a++) generator(a, iTa[a]);
	}

	Real coefficient = - 2.0/(ncolour + 2*S) * scale; // 2/(Nc +/- 2) for the normalization of the trace in the two index rep

	for (int a = 0; a < Dimension; a++) {
	  auto tmp = real(trace(iTa[a] * in)) * coefficient; 
	  pokeColour(h_out, tmp, a);
	}
      }


    };




    // Some useful type names
    typedef SU_TwoIndex<2, Symmetric> SU2TwoIndexSymm;
    typedef SU_TwoIndex<3, Symmetric> SU3TwoIndexSymm;
    typedef SU_TwoIndex<4, Symmetric> SU4TwoIndexSymm;
    typedef SU_TwoIndex<5, Symmetric> SU5TwoIndexSymm;

    typedef SU_TwoIndex<Nc, Symmetric> TwoIndexSymmMatrices;

    typedef SU_TwoIndex<2, AntiSymmetric> SU2TwoIndexAntiSymm;
    typedef SU_TwoIndex<3, AntiSymmetric> SU3TwoIndexAntiSymm;
    typedef SU_TwoIndex<4, AntiSymmetric> SU4TwoIndexAntiSymm;
    typedef SU_TwoIndex<5, AntiSymmetric> SU5TwoIndexAntiSymm;

    typedef SU_TwoIndex<Nc, AntiSymmetric> TwoIndexAntiSymmMatrices;

    
  }
}

#endif
