/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/utils/GaugeGroup.h

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
#ifndef QCD_UTIL_GAUGEGROUP_H
#define QCD_UTIL_GAUGEGROUP_H

// Important detail: nvcc requires all template parameters to have names.
// This is the only reason why the second template parameter has a name.
#define ONLY_IF_SU                                                       \
  typename dummy_name = group_name,                                      \
           typename named_dummy = std::enable_if_t <                                 \
                          std::is_same<dummy_name, group_name>::value && \
                      is_su<dummy_name>::value >

#define ONLY_IF_Sp                                                       \
  typename dummy_name = group_name,                                      \
           typename named_dummy = std::enable_if_t <                                 \
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
  template <typename vtype>
  using iSUnAlgebraMatrix =
    iScalar<iScalar<iMatrix<vtype, AdjointDimension> > >;
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
   
  typedef iSUnAlgebraMatrix<vComplex>  vAlgebraMatrix;
  typedef iSUnAlgebraMatrix<vComplexF> vAlgebraMatrixF;
  typedef iSUnAlgebraMatrix<vComplexD> vAlgebraMatrixD;

  typedef Lattice<vAlgebraMatrix>  LatticeAlgebraMatrix;
  typedef Lattice<vAlgebraMatrixF> LatticeAlgebraMatrixF;
  typedef Lattice<vAlgebraMatrixD> LatticeAlgebraMatrixD;
  

  typedef iSU2Matrix<Complex> SU2Matrix;
  typedef iSU2Matrix<ComplexF> SU2MatrixF;
  typedef iSU2Matrix<ComplexD> SU2MatrixD;

  typedef iSU2Matrix<vComplex> vSU2Matrix;
  typedef iSU2Matrix<vComplexF> vSU2MatrixF;
  typedef iSU2Matrix<vComplexD> vSU2MatrixD;

  typedef Lattice<vSU2Matrix> LatticeSU2Matrix;
  typedef Lattice<vSU2MatrixF> LatticeSU2MatrixF;
  typedef Lattice<vSU2MatrixD> LatticeSU2MatrixD;

  // Private implementation details are specified in the following files:
  // Grid/qcd/utils/SUn.impl
  // Grid/qcd/utils/SUn.impl
  // The public part of the interface follows below and refers to these
  // private member functions.

#include <Grid/qcd/utils/SUn.impl.h>
#include <Grid/qcd/utils/Sp2n.impl.h>

 public:
  template <class cplx>
  static void generator(int lieIndex, iGroupMatrix<cplx> &ta) {
    return generator(lieIndex, ta, group_name());
  }

  static accelerator_inline void su2SubGroupIndex(int &i1, int &i2, int su2_index) {
    return su2SubGroupIndex(i1, i2, su2_index, group_name());
  }

  static void testGenerators(void) { testGenerators(group_name()); }

  static void printGenerators(void) {
    for (int gen = 0; gen < AlgebraDimension; gen++) {
      Matrix ta;
      generator(gen, ta);
      std::cout << GridLogMessage << "Nc = " << ncolour << " t_" << gen
                << std::endl;
      std::cout << GridLogMessage << ta << std::endl;
    }
  }

  template <typename LatticeMatrixType>
  static void LieRandomize(GridParallelRNG &pRNG, LatticeMatrixType &out,
                           double scale = 1.0) {
    GridBase *grid = out.Grid();

    typedef typename LatticeMatrixType::vector_type vector_type;

    typedef iSinglet<vector_type> vTComplexType;

    typedef Lattice<vTComplexType> LatticeComplexType;
    typedef typename GridTypeMapper<
        typename LatticeMatrixType::vector_object>::scalar_object MatrixType;

    LatticeComplexType ca(grid);
    LatticeMatrixType lie(grid);
    LatticeMatrixType la(grid);
    ComplexD ci(0.0, scale);
    MatrixType ta;

    lie = Zero();

    for (int a = 0; a < AlgebraDimension; a++) {
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
    for (int a = 0; a < AlgebraDimension; a++) {
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
    for (int a = 0; a < AlgebraDimension; a++) {
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

    for (int a = 0; a < AlgebraDimension; a++) {
      generator(a, Ta);
      pokeColour(h_out, -2.0 * (trace(timesI(Ta) * in)) * scale, a);
    }
  }

   
  template <class vtype>
  accelerator_inline static iScalar<vtype> ProjectOnGeneralGroup(const iScalar<vtype> &r) {
    return ProjectOnGeneralGroup(r, group_name());
  }

  template <class vtype, int N>
  accelerator_inline static iVector<vtype,N> ProjectOnGeneralGroup(const iVector<vtype,N> &r) {
    return ProjectOnGeneralGroup(r, group_name());
  }

  template <class vtype,int N, typename std::enable_if< GridTypeMapper<vtype>::TensorLevel == 0 >::type * =nullptr>
  accelerator_inline static iMatrix<vtype,N> ProjectOnGeneralGroup(const iMatrix<vtype,N> &arg) {
    return ProjectOnGeneralGroup(arg, group_name());
  }

  template <int N,class vComplex_t>                  // Projects on the general groups U(N), Sp(2N)xZ2 i.e. determinant is allowed a complex phase.
  static void ProjectOnGeneralGroup(Lattice<iVector<iScalar<iMatrix<vComplex_t, N> >, Nd> > &U) {
    for (int mu = 0; mu < Nd; mu++) {
      auto Umu = PeekIndex<LorentzIndex>(U, mu);
      Umu = ProjectOnGeneralGroup(Umu);
    }
  }
       

  
  template <int N,class vComplex_t>
  static Lattice<iScalar<iScalar<iMatrix<vComplex_t, N> > > > ProjectOnGeneralGroup(const Lattice<iScalar<iScalar<iMatrix<vComplex_t, N> > > > &Umu) {
    return ProjectOnGeneralGroup(Umu, group_name());
  }

  template <int N,class vComplex_t>       // Projects on SU(N), Sp(2N), with unit determinant, by first projecting on general group and then enforcing unit determinant
  static void ProjectOnSpecialGroup(Lattice<iScalar<iScalar<iMatrix<vComplex_t, N> > > > &Umu) {
       Umu = ProjectOnGeneralGroup(Umu);
       auto det = Determinant(Umu);

       det = conjugate(det);

       for (int i = 0; i < N; i++) {
           auto element = PeekIndex<ColourIndex>(Umu, N - 1, i);
           element = element * det;
           PokeIndex<ColourIndex>(Umu, element, Nc - 1, i);
       }
   }

  template <int N,class vComplex_t>    // reunitarise, resimplectify... previously ProjectSUn
    static void ProjectOnSpecialGroup(Lattice<iVector<iScalar<iMatrix<vComplex_t, N> >, Nd> > &U) {
      // Reunitarise
      for (int mu = 0; mu < Nd; mu++) {
        auto Umu = PeekIndex<LorentzIndex>(U, mu);
        ProjectOnSpecialGroup(Umu);
        PokeIndex<LorentzIndex>(U, Umu, mu);
      }
    }
    
  template <typename GaugeField>
  static void HotConfiguration(GridParallelRNG &pRNG, GaugeField &out) {
    typedef typename GaugeField::vector_type vector_type;
    typedef iGroupMatrix<vector_type> vMatrixType;
    typedef Lattice<vMatrixType> LatticeMatrixType;

    LatticeMatrixType Umu(out.Grid());
    LatticeMatrixType tmp(out.Grid());
    for (int mu = 0; mu < Nd; mu++) {
      //      LieRandomize(pRNG, Umu, 1.0);
      //      PokeIndex<LorentzIndex>(out, Umu, mu);
      gaussian(pRNG,Umu);
      tmp = Ta(Umu);
      taExp(tmp,Umu);
      ProjectOnSpecialGroup(Umu);
      //      ProjectSUn(Umu);
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

  template <typename LatticeMatrixType>
  static void taProj(const LatticeMatrixType &in, LatticeMatrixType &out) {
    taProj(in, out, group_name());
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

// Ta are hermitian (?)
// Anti herm is i Ta basis
static void LieAlgebraProject(LatticeAlgebraMatrix &out,const LatticeMatrix &in, int b)
{
  conformable(in, out);
  GridBase *grid = out.Grid();
  LatticeComplex tmp(grid);
  Matrix ta;
  // Using Luchang's projection convention
  //  2 Tr{Ta Tb} A_b= 2/2 delta ab A_b = A_a
  autoView(out_v,out,AcceleratorWrite);
  autoView(in_v,in,AcceleratorRead);
  int N = ncolour;
  int NNm1 = N * (N - 1);
  int hNNm1= NNm1/2;
  RealD sqrt_2 = sqrt(2.0);
  Complex ci(0.0,1.0);

  const int nsimd=  Matrix::Nsimd();
  accelerator_for(ss,grid->oSites(),nsimd,{
      for(int su2Index=0;su2Index<hNNm1;su2Index++){
	int i1, i2;
	su2SubGroupIndex(i1, i2, su2Index);
	int ax = su2Index*2;
	int ay = su2Index*2+1;
	// in is traceless ANTI-hermitian whereas Grid generators are Hermitian.
	// trace( Ta x Ci in)
	// Bet I need to move to real part with mult by -i
	coalescedWrite(out_v[ss]()()(ax,b),0.5*(real(in_v(ss)()()(i2,i1)) - real(in_v(ss)()()(i1,i2))));
	coalescedWrite(out_v[ss]()()(ay,b),0.5*(imag(in_v(ss)()()(i1,i2)) + imag(in_v(ss)()()(i2,i1))));
      }
      for(int diagIndex=0;diagIndex<N-1;diagIndex++){
	int k = diagIndex + 1; // diagIndex starts from 0
	int a = NNm1+diagIndex;
	RealD scale = 1.0/sqrt(2.0*k*(k+1));
	auto tmp = in_v(ss)()()(0,0);
	for(int i=1;i<k;i++){
	  tmp=tmp+in_v(ss)()()(i,i);
	}
	tmp = tmp - in_v(ss)()()(k,k)*k;
	coalescedWrite(out_v[ss]()()(a,b),imag(tmp) * scale);
      }
    });
}

  
};
    
template <int ncolour>
using SU = GaugeGroup<ncolour, GroupName::SU>;

template <int ncolour>
using Sp = GaugeGroup<ncolour, GroupName::Sp>;

typedef SU<2> SU2;
typedef SU<3> SU3;
typedef SU<4> SU4;
typedef SU<5> SU5;

typedef SU<Nc> FundamentalMatrices;
    
typedef Sp<2> Sp2;
typedef Sp<4> Sp4;
typedef Sp<6> Sp6;
typedef Sp<8> Sp8;

template <int N,class vComplex_t>
static void ProjectSUn(Lattice<iScalar<iScalar<iMatrix<vComplex_t, N> > > > &Umu)
{
    GaugeGroup<N,GroupName::SU>::ProjectOnSpecialGroup(Umu);
}
  
template <int N,class vComplex_t>
static void ProjectSUn(Lattice<iVector<iScalar<iMatrix<vComplex_t, N> >,Nd> > &U)
{
    GaugeGroup<N,GroupName::SU>::ProjectOnSpecialGroup(U);
}
    
template <int N,class vComplex_t>
static void ProjectSpn(Lattice<iScalar<iScalar<iMatrix<vComplex_t, N> > > > &Umu)
{
    GaugeGroup<N,GroupName::Sp>::ProjectOnSpecialGroup(Umu);
}
    
template <int N,class vComplex_t>
static void ProjectSpn(Lattice<iVector<iScalar<iMatrix<vComplex_t, N> >,Nd> > &U)
{
    GaugeGroup<N,GroupName::Sp>::ProjectOnSpecialGroup(U);
}

// Explicit specialisation for SU(3).
static void ProjectSU3(Lattice<iScalar<iScalar<iMatrix<vComplexD, 3> > > > &Umu)
{
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
static void ProjectSU3(Lattice<iVector<iScalar<iMatrix<vComplexD, 3> >, Nd> > &U)
{
  GridBase *grid = U.Grid();
  // Reunitarise
  for (int mu = 0; mu < Nd; mu++) {
    auto Umu = PeekIndex<LorentzIndex>(U, mu);
    Umu = ProjectOnGroup(Umu);
    ProjectSU3(Umu);
    PokeIndex<LorentzIndex>(U, Umu, mu);
  }
}

NAMESPACE_END(Grid);
#endif
