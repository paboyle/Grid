#ifndef QCD_UTIL_SUNADJOINT_H
#define QCD_UTIL_SUNADJOINT_H

////////////////////////////////////////////////////////////////////////
//
// * Adjoint representation generators
//
// * Normalisation for the fundamental generators: 
//   trace ta tb = 1/2 delta_ab = T_F delta_ab
//   T_F = 1/2  for SU(N) groups
//
//
//   base for NxN hermitian traceless matrices
//   normalized to 1:
//
//   (e_Adj)^a = t^a / sqrt(T_F)
//
//   then the real, antisymmetric generators for the adjoint representations
//   are computed ( shortcut: e^a == (e_Adj)^a )
//
//   (iT_adj)^d_ba = i tr[e^a t^d e^b - t^d e^a e^b]
//
////////////////////////////////////////////////////////////////////////

NAMESPACE_BEGIN(Grid);

template <int ncolour>
class SU_Adjoint : public SU<ncolour> {
public:
  static const int Dimension = ncolour * ncolour - 1;

  template <typename vtype>
  using iSUnAdjointMatrix =
    iScalar<iScalar<iMatrix<vtype, Dimension > > >;

  // Actually the adjoint matrices are real...
  // Consider this overhead... FIXME
  typedef iSUnAdjointMatrix<Complex> AMatrix;
  typedef iSUnAdjointMatrix<ComplexF> AMatrixF;
  typedef iSUnAdjointMatrix<ComplexD> AMatrixD;

  typedef iSUnAdjointMatrix<vComplex>  vAMatrix;
  typedef iSUnAdjointMatrix<vComplexF> vAMatrixF;
  typedef iSUnAdjointMatrix<vComplexD> vAMatrixD;

  typedef Lattice<vAMatrix>  LatticeAdjMatrix;
  typedef Lattice<vAMatrixF> LatticeAdjMatrixF;
  typedef Lattice<vAMatrixD> LatticeAdjMatrixD;

  typedef Lattice<iVector<iScalar<iMatrix<vComplex, Dimension> >, Nd> >  LatticeAdjField;
  typedef Lattice<iVector<iScalar<iMatrix<vComplexF, Dimension> >, Nd> > LatticeAdjFieldF;
  typedef Lattice<iVector<iScalar<iMatrix<vComplexD, Dimension> >, Nd> > LatticeAdjFieldD;


  template <typename vtype>
  using iSUnMatrix = iScalar<iScalar<iMatrix<vtype, ncolour> > >;

  typedef Lattice<iScalar<iScalar<iVector<vComplex, Dimension> > > >  LatticeAdjVector;

  template <class cplx>
  static accelerator_inline void generator(int Index, iSUnAdjointMatrix<cplx> &iAdjTa) {
    // returns i(T_Adj)^index necessary for the projectors
    // see definitions above
    iAdjTa = Zero();

    iVector<iSUnMatrix<cplx>,Dimension> ta;
    iSUnMatrix<cplx> tmp;

    // FIXME not very efficient to get all the generators everytime
    for (int a = 0; a < Dimension; a++) SU<ncolour>::generator(a, ta(a));

    for (int a = 0; a < Dimension; a++) {
      tmp = ta(a) * ta(Index) - ta(Index) * ta(a);
      for (int b = 0; b < Dimension; b++) {
        iSUnMatrix<cplx> tmp1 = 2.0 * tmp * ta(b);  // 2.0 from the normalization
        Complex iTr = TensorRemove(timesI(trace(tmp1)));
	//cplx iTr = TensorRemove(timesI(trace(tmp1)));// does not work if cplx is Grid_simd type, as no constructor for cplx i(0.0, 1.0); in generatorSigmaX is available
        //iAdjTa()()(b, a) = iTr;
        iAdjTa()()(a, b) = iTr;
      }
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
    AMatrix adjTa;
    std::cout << GridLogMessage << "Adjoint - Checking if real" << std::endl;
    for (int a = 0; a < Dimension; a++) {
      generator(a, adjTa);
      std::cout << GridLogMessage << a << std::endl;
      assert(norm2(adjTa - conjugate(adjTa)) < 1.0e-6);
    }
    std::cout << GridLogMessage << std::endl;

    std::cout << GridLogMessage << "Adjoint - Checking if antisymmetric"
              << std::endl;
    for (int a = 0; a < Dimension; a++) {
      generator(a, adjTa);
      std::cout << GridLogMessage << a << std::endl;
      assert(norm2(adjTa + transpose(adjTa)) < 1.0e-6);
    }
    std::cout << GridLogMessage << std::endl;
  }

  static void AdjointLieAlgebraMatrix(
				      const typename SU<ncolour>::LatticeAlgebraVector &h,
				      LatticeAdjMatrix &out, Real scale = 1.0) {
    conformable(h, out);
    GridBase *grid = out.Grid();
    LatticeAdjMatrix la(grid);
    AMatrix iTa;

    out = Zero();
    for (int a = 0; a < Dimension; a++) {
      generator(a, iTa);
      la = peekColour(h, a) * iTa;
      out += la;
    }
    out *= scale;
  }

  // Projects the algebra components a lattice matrix (of dimension ncol*ncol -1 )
  static void projectOnAlgebra(typename SU<ncolour>::LatticeAlgebraVector &h_out, const LatticeAdjMatrix &in, Real scale = 1.0) 
  {
    
    conformable(h_out, in);
    h_out = Zero();
    AMatrix iTa;
    Real coefficient = - 1.0/(ncolour) * scale;// 1/Nc for the normalization of the trace in the adj rep

    for (int a = 0; a < Dimension; a++) {
      generator(a, iTa);
      pokeColour(h_out, real(trace(iTa * in)) * coefficient, a);
    }
  }

  // a projector that keeps the generators stored to avoid the overhead of recomputing them 
  static void projector(typename SU<ncolour>::LatticeAlgebraVector &h_out, const LatticeAdjMatrix &in, Real scale = 1.0) {
    conformable(h_out, in);
    static std::vector<AMatrix> iTa(Dimension);  // to store the generators
    h_out = Zero();
    static bool precalculated = false; 
    if (!precalculated){
      precalculated = true;
      for (int a = 0; a < Dimension; a++) generator(a, iTa[a]);
    }

    Real coefficient = -1.0 / (ncolour) * scale;  // 1/Nc for the normalization of
    // the trace in the adj rep

    for (int a = 0; a < Dimension; a++) {
      auto tmp = real(trace(iTa[a] * in)) * coefficient; 
      pokeColour(h_out, tmp, a);
    }
  }
#if 1
  // Turn complex 3x3 Lie algebra elem into its real 8x8 adj rep
  static void make_adjoint_rep(LatticeAdjMatrix &out, const typename SU<ncolour>::LatticeMatrix &in) {
    // Use real and totally anti-symmetric matrices as adjoint rep
    //  => basis matrices for su(N) are anti-hermitian
    // Take: T^a = i t^a where t^a: basis elem in Grid's convention; T^a satisfies Luscher's normalization convention -2tr[T^a,T^b] = \delta_{ab}
    // Recall: adj rep generated above, iAdjTa, is multiplied by i
    //   This makes adj rep here is consistent with the one written in terms of T^a
    // So for adj rep, we follow Luscher's convention
    
    // in: traceless-anti-hermitian

    GridBase *grid = out.Grid();
    AMatrix iTa;
      
    autoView(out_v,out,AcceleratorWrite);
    autoView(in_v,in,AcceleratorRead);
    int N = ncolour;
    int NNm1 = N * (N - 1);
    int hNNm1= NNm1/2;
    const int nsimd = vAMatrix::Nsimd();
    accelerator_for(ss,grid->oSites(),nsimd,{
	typedef decltype(coalescedRead(out_v[0])) adj_mat;
	typedef decltype(coalescedRead(vTComplex())) computeComplex;
	
	adj_mat Tadj=Zero();
	computeComplex c;
	
        for(int su2Index=0;su2Index<hNNm1;su2Index++){
          int i1, i2;
          SU<ncolour>::su2SubGroupIndex(i1, i2, su2Index);
          int ax = su2Index*2;
          int ay = su2Index*2+1;

	  // multiply t^a by -i
	  // cplx = -2.0*trace(ci*tb*Zx);
	  // ZxAd = ZxAd + cplx * TRb;
	  generator(ax,iTa);
	  c()()() = 2.0*real(in_v(ss)()()(i2,i1));
	  Tadj = Tadj + c*iTa;
	  generator(ay,iTa);
          c()()() = 2.0*imag(in_v(ss)()()(i1,i2));
	  Tadj = Tadj + c*iTa;
        }
        for(int diagIndex=0;diagIndex<N-1;diagIndex++){
          int k = diagIndex + 1; // diagIndex starts from 0
          int a = NNm1+diagIndex;
          RealD scale = 1.0/sqrt(2.0*k*(k+1));

	  generator(a,iTa);
          auto tmp = in_v(ss)()()(0,0);
          for(int i=1;i<k;i++){
            tmp=tmp+in_v(ss)()()(i,i);
          }
          tmp = tmp - in_v(ss)()()(k,k)*k;
	  c()()() = 2.0 * scale * imag(tmp);
          Tadj = Tadj + c*iTa;
        }
	coalescedWrite(out_v[ss],Tadj);
      });
  }
#endif
};




// Some useful type names

typedef SU_Adjoint<2> SU2Adjoint;
typedef SU_Adjoint<3> SU3Adjoint;
typedef SU_Adjoint<4> SU4Adjoint;
typedef SU_Adjoint<5> SU5Adjoint;

typedef SU_Adjoint<Nc> AdjointMatrices;

NAMESPACE_END(Grid);

#endif
