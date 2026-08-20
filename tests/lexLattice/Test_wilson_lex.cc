/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/lexLattice/Test_wilson_lex.cc

    Copyright (C) 2026

Author: Peter Boyle <pboyle@bnl.gov>

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

//
// Wilson fermion operator on a lexicographic (Nsimd()==1) lattice.
//
// Exercises the spin-projected half spinor path -- WilsonCompressor and
// WilsonStencil -- which the staggered and stencil tests do not reach.
//
//   - lexicographic Dhop against a covariant-shift reference in the same chart
//   - both charts from one gauge field and source, compared elementwise
//   - gamma5 hermiticity and the even-odd hopping term in the lex chart
//
#include <Grid/Grid.h>

using namespace Grid;

const RealD tol = 1.0e-10;

Gamma::Algebra Gmu [] = {
  Gamma::Algebra::GammaX,
  Gamma::Algebra::GammaY,
  Gamma::Algebra::GammaZ,
  Gamma::Algebra::GammaT
};

typedef WilsonFermion<WilsonImplD>    vWilsonOp;
typedef WilsonFermion<lexWilsonImplD> lexWilsonOp;

////////////////////////////////////////////////////////////////////////
// Layout interchange; both lattices share the same scalar_object.
////////////////////////////////////////////////////////////////////////
template<class vobjOut,class vobjIn>
void transfer(Lattice<vobjOut> &out,const Lattice<vobjIn> &in)
{
  typedef typename vobjIn::scalar_object sobj;
  static_assert(std::is_same<sobj,typename vobjOut::scalar_object>::value,
		"transfer: lattices must share scalar_object");

  GRID_ASSERT(out.Grid()->gSites() == in.Grid()->gSites());

  std::vector<sobj> buf;
  unvectorizeToLexOrdArray(buf,in);
  vectorizeFromLexOrdArray(buf,out);
}

////////////////////////////////////////////////////////////////////////
// Wilson hopping term from covariant shifts, in whichever chart Impl names
////////////////////////////////////////////////////////////////////////
template<class Impl>
void ReferenceDhop(typename Impl::FermionField &ref,
		   const typename Impl::GaugeField &Umu,
		   const typename Impl::FermionField &src)
{
  typedef typename Impl::GaugeLinkField LinkField;
  typedef typename Impl::FermionField   FermionField;

  GridBase *grid = src.Grid();

  std::vector<LinkField> U(Nd,grid);
  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }

  FermionField tmp(grid);
  ref = Zero();

  for(int mu=0;mu<Nd;mu++){
    tmp = U[mu]*Cshift(src,mu,1);
    ref = ref + tmp - Gamma(Gmu[mu])*tmp;

    tmp = adj(U[mu])*src;
    tmp = Cshift(tmp,mu,-1);
    ref = ref + tmp + Gamma(Gmu[mu])*tmp;
  }
  ref = -0.5*ref;
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Coordinate latt  = GridDefaultLatt();
  Coordinate vsimd = GridDefaultSimd(Nd,vComplexD::Nsimd());
  Coordinate lsimd({1,1,1,1});
  Coordinate mpi   = GridDefaultMpi();

  GridCartesian         vGrid(latt,vsimd,mpi);
  GridRedBlackCartesian vRBGrid(&vGrid);
  GridCartesian         lGrid(latt,lsimd,mpi);
  GridRedBlackCartesian lRBGrid(&lGrid);

  std::cout << GridLogMessage << "vectorised    Nsimd = " << vGrid.Nsimd() << std::endl;
  std::cout << GridLogMessage << "lexicographic Nsimd = " << lGrid.Nsimd() << std::endl;
  GRID_ASSERT( lGrid.Nsimd() == 1 );

  RealD mass = 0.1;

  GridParallelRNG pRNG(&vGrid);
  pRNG.SeedFixedIntegers(std::vector<int>({1,2,3,4}));

  //////////////////////////////////////////////////
  // One gauge field and one source, both charts
  //////////////////////////////////////////////////
  typename vWilsonOp::GaugeField   Umu(&vGrid);  SU<Nc>::HotConfiguration(pRNG,Umu);
  typename vWilsonOp::FermionField src(&vGrid);  random(pRNG,src);
  typename vWilsonOp::FermionField phi(&vGrid);  random(pRNG,phi);

  typename lexWilsonOp::GaugeField   Ulex(&lGrid);   transfer(Ulex,Umu);
  typename lexWilsonOp::FermionField srclex(&lGrid); transfer(srclex,src);
  typename lexWilsonOp::FermionField philex(&lGrid); transfer(philex,phi);

  vWilsonOp   Dw   (Umu ,vGrid,vRBGrid,mass);
  lexWilsonOp Dwlex(Ulex,lGrid,lRBGrid,mass);

  //////////////////////////////////////////////////
  // Lexicographic Dhop against covariant shifts
  //////////////////////////////////////////////////
  {
    typename lexWilsonOp::FermionField ref(&lGrid);
    typename lexWilsonOp::FermionField res(&lGrid);
    typename lexWilsonOp::FermionField err(&lGrid);

    ReferenceDhop<lexWilsonImplD>(ref,Ulex,srclex);
    Dwlex.Dhop(srclex,res,DaggerNo);

    err = res - ref;
    RealD n = norm2(err)/norm2(ref);
    std::cout << GridLogMessage << "lex Dhop vs covariant shifts: relative " << n
	      << " (|Dhop|^2 " << norm2(res) << ")" << std::endl;
    GRID_ASSERT( n < tol );
  }

  //////////////////////////////////////////////////
  // Chart equivalence, elementwise
  //////////////////////////////////////////////////
  {
    typename vWilsonOp::FermionField   vres(&vGrid);
    typename lexWilsonOp::FermionField lres(&lGrid);
    typename lexWilsonOp::FermionField vres_lex(&lGrid);
    typename lexWilsonOp::FermionField err(&lGrid);

    Dw.M(src,vres);   Dwlex.M(srclex,lres);
    transfer(vres_lex,vres);
    err = vres_lex - lres;
    RealD n = norm2(err)/norm2(lres);
    std::cout << GridLogMessage << "M simd vs lex: relative " << n
	      << "  |simd|^2 " << norm2(vres) << "  |lex|^2 " << norm2(lres) << std::endl;
    GRID_ASSERT( n < tol );

    Dw.Mdag(src,vres);   Dwlex.Mdag(srclex,lres);
    transfer(vres_lex,vres);
    err = vres_lex - lres;
    n = norm2(err)/norm2(lres);
    std::cout << GridLogMessage << "Mdag simd vs lex: relative " << n << std::endl;
    GRID_ASSERT( n < tol );

    Dw.Dhop(src,vres,DaggerNo);   Dwlex.Dhop(srclex,lres,DaggerNo);
    transfer(vres_lex,vres);
    err = vres_lex - lres;
    n = norm2(err)/norm2(lres);
    std::cout << GridLogMessage << "Dhop simd vs lex: relative " << n << std::endl;
    GRID_ASSERT( n < tol );
  }

  //////////////////////////////////////////////////
  // gamma5 hermiticity in the lexicographic chart
  //   g5 M g5 = Mdag
  //////////////////////////////////////////////////
  {
    typename lexWilsonOp::FermionField tmp(&lGrid);
    typename lexWilsonOp::FermionField g5Mg5(&lGrid);
    typename lexWilsonOp::FermionField Mdag(&lGrid);
    typename lexWilsonOp::FermionField err(&lGrid);

    tmp = Gamma(Gamma::Algebra::Gamma5)*srclex;
    Dwlex.M(tmp,g5Mg5);
    g5Mg5 = Gamma(Gamma::Algebra::Gamma5)*g5Mg5;

    Dwlex.Mdag(srclex,Mdag);

    err = g5Mg5 - Mdag;
    RealD n = norm2(err)/norm2(Mdag);
    std::cout << GridLogMessage << "lex g5 M g5 - Mdag: relative " << n << std::endl;
    GRID_ASSERT( n < tol );
  }

  //////////////////////////////////////////////////
  // Checkerboarded hopping term, both charts
  //////////////////////////////////////////////////
  {
    typename vWilsonOp::FermionField   vsrc_o(&vRBGrid);   vsrc_o.Checkerboard()=Odd;
    typename vWilsonOp::FermionField   vres_e(&vRBGrid);   vres_e.Checkerboard()=Even;
    typename lexWilsonOp::FermionField lsrc_o(&lRBGrid);   lsrc_o.Checkerboard()=Odd;
    typename lexWilsonOp::FermionField lres_e(&lRBGrid);   lres_e.Checkerboard()=Even;

    pickCheckerboard(Odd,vsrc_o,src);
    pickCheckerboard(Odd,lsrc_o,srclex);

    Dw.Meooe   (vsrc_o,vres_e);
    Dwlex.Meooe(lsrc_o,lres_e);

    RealD vn = norm2(vres_e);
    RealD ln = norm2(lres_e);
    RealD n  = fabs(vn-ln)/vn;
    std::cout << GridLogMessage << "Meooe |simd|^2 " << vn << "  |lex|^2 " << ln
	      << "  relative " << n << std::endl;
    GRID_ASSERT( n < tol );
  }

  std::cout << GridLogMessage << "Test_wilson_lex: ALL PASS" << std::endl;

  Grid_finalize();
}
