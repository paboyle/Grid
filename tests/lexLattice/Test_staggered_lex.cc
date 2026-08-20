/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/lexLattice/Test_staggered_lex.cc

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
// Naive staggered fermion operator on a lexicographic (Nsimd()==1) lattice.
//
// Two independent checks:
//   - lexicographic Dhop against a covariant-shift reference built in the
//     same chart
//   - both charts driven from the same gauge field and source, compared
//     elementwise through the layout interchange
//
#include <Grid/Grid.h>

using namespace Grid;

const RealD tol = 1.0e-10;

typedef NaiveStaggeredFermion<StaggeredImplD>    vStagOp;
typedef NaiveStaggeredFermion<lexStaggeredImplD> lexStagOp;

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
// Dhop from covariant shifts, in whichever chart Gimpl names
////////////////////////////////////////////////////////////////////////
template<class Impl>
void ReferenceDhop(typename Impl::FermionField &ref,
		   const typename Impl::GaugeField &Umu,
		   const typename Impl::FermionField &src,
		   RealD c1,RealD u0)
{
  typedef typename Impl::GaugeLinkField LinkField;
  typedef typename Impl::ComplexField   ComplexField;
  typedef typename Impl::FermionField   FermionField;

  GridBase *grid = src.Grid();
  RealD c1tad = 0.5*c1/u0;

  std::vector<LinkField> U(Nd,grid);
  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }

  typedef Lattice<iScalar<typename GridTypeMapper<typename Impl::Simd>::Integerified> > IntField;

  IntField x(grid); LatticeCoordinate(x,0);
  IntField y(grid); LatticeCoordinate(y,1);
  IntField z(grid); LatticeCoordinate(z,2);
  IntField lin_z(grid); lin_z = x+y;
  IntField lin_t(grid); lin_t = x+y+z;

  FermionField tmp(grid);
  ref = Zero();

  for(int mu=0;mu<Nd;mu++){

    ComplexField phases(grid); phases=1.0;
    if ( mu == 1 ) phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
    if ( mu == 2 ) phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
    if ( mu == 3 ) phases = where( mod(lin_t,2)==(Integer)0, phases,-phases);

    tmp = PeriodicBC::CovShiftForward(U[mu],mu,src);
    ref = ref + c1tad*tmp*phases;

    tmp = PeriodicBC::CovShiftBackward(U[mu],mu,src);
    ref = ref - c1tad*tmp*phases;
  }
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
  RealD c1   = 9.0/8.0;
  RealD u0   = 1.0;

  GridParallelRNG pRNG(&vGrid);
  pRNG.SeedFixedIntegers(std::vector<int>({1,2,3,4}));

  //////////////////////////////////////////////////
  // One gauge field and one source, both charts
  //////////////////////////////////////////////////
  typename vStagOp::GaugeField   Umu(&vGrid);   SU<Nc>::HotConfiguration(pRNG,Umu);
  typename vStagOp::FermionField src(&vGrid);   random(pRNG,src);
  typename vStagOp::FermionField phi(&vGrid);   random(pRNG,phi);

  typename lexStagOp::GaugeField   Ulex(&lGrid);  transfer(Ulex,Umu);
  typename lexStagOp::FermionField srclex(&lGrid); transfer(srclex,src);
  typename lexStagOp::FermionField philex(&lGrid); transfer(philex,phi);

  typename vStagOp::ImplParams   vparams;
  typename lexStagOp::ImplParams lparams;

  vStagOp   Ds   (Umu ,vGrid,vRBGrid,mass,c1,u0,vparams);
  lexStagOp Dslex(Ulex,lGrid,lRBGrid,mass,c1,u0,lparams);

  //////////////////////////////////////////////////
  // Lexicographic Dhop against covariant shifts
  //////////////////////////////////////////////////
  {
    typename lexStagOp::FermionField ref(&lGrid);
    typename lexStagOp::FermionField res(&lGrid);
    typename lexStagOp::FermionField err(&lGrid);

    ReferenceDhop<lexStaggeredImplD>(ref,Ulex,srclex,c1,u0);
    Dslex.Dhop(srclex,res,DaggerNo);

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
    typename vStagOp::FermionField   vres(&vGrid);
    typename lexStagOp::FermionField lres(&lGrid);
    typename lexStagOp::FermionField vres_lex(&lGrid);
    typename lexStagOp::FermionField err(&lGrid);

    struct { const char *name; int dag; } ops[2] = { {"M",DaggerNo}, {"Mdag",DaggerYes} };

    for(int o=0;o<2;o++){
      if ( ops[o].dag == DaggerNo ) { Ds.M   (src,vres); Dslex.M   (srclex,lres); }
      else                          { Ds.Mdag(src,vres); Dslex.Mdag(srclex,lres); }

      transfer(vres_lex,vres);
      err = vres_lex - lres;
      RealD n = norm2(err)/norm2(lres);
      std::cout << GridLogMessage << ops[o].name << " simd vs lex: relative " << n
		<< "  |simd|^2 " << norm2(vres) << "  |lex|^2 " << norm2(lres) << std::endl;
      GRID_ASSERT( n < tol );
    }

    Ds.Dhop(src,vres,DaggerNo);  Dslex.Dhop(srclex,lres,DaggerNo);
    transfer(vres_lex,vres);
    err = vres_lex - lres;
    RealD n = norm2(err)/norm2(lres);
    std::cout << GridLogMessage << "Dhop simd vs lex: relative " << n << std::endl;
    GRID_ASSERT( n < tol );
  }

  //////////////////////////////////////////////////
  // Adjoint identity within the lexicographic chart
  //////////////////////////////////////////////////
  {
    typename lexStagOp::FermionField Mphi(&lGrid);
    typename lexStagOp::FermionField Mdagchi(&lGrid);

    Dslex.M   (srclex,Mphi);
    Dslex.Mdag(philex,Mdagchi);

    ComplexD lhs = innerProduct(philex,Mphi);
    ComplexD rhs = innerProduct(Mdagchi,srclex);
    RealD n = abs(lhs-rhs)/abs(lhs);
    std::cout << GridLogMessage << "lex <phi|M src> = " << lhs
	      << "  <Mdag phi|src> = " << rhs << "  relative " << n << std::endl;
    GRID_ASSERT( n < tol );
  }

  //////////////////////////////////////////////////
  // Checkerboarded hopping term, both charts
  //////////////////////////////////////////////////
  {
    typename vStagOp::FermionField   vsrc_o(&vRBGrid);   vsrc_o.Checkerboard()=Odd;
    typename vStagOp::FermionField   vres_e(&vRBGrid);   vres_e.Checkerboard()=Even;
    typename lexStagOp::FermionField lsrc_o(&lRBGrid);   lsrc_o.Checkerboard()=Odd;
    typename lexStagOp::FermionField lres_e(&lRBGrid);   lres_e.Checkerboard()=Even;

    pickCheckerboard(Odd,vsrc_o,src);
    pickCheckerboard(Odd,lsrc_o,srclex);

    Ds.Meooe   (vsrc_o,vres_e);
    Dslex.Meooe(lsrc_o,lres_e);

    RealD vn = norm2(vres_e);
    RealD ln = norm2(lres_e);
    RealD n  = fabs(vn-ln)/vn;
    std::cout << GridLogMessage << "Meooe |simd|^2 " << vn << "  |lex|^2 " << ln
	      << "  relative " << n << std::endl;
    GRID_ASSERT( n < tol );
  }

  std::cout << GridLogMessage << "Test_staggered_lex: ALL PASS" << std::endl;

  Grid_finalize();
}
