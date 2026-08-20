/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/lexLattice/Test_dwf_lex.cc

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
// Five dimensional fermions on a lexicographic (Nsimd()==1) lattice.
//
// Domain wall and Mobius, driven from one gauge field and one source in both
// charts and compared elementwise through the layout interchange.  Mobius uses
// b=1.5 c=0.5, the production choice.
//
#include <Grid/Grid.h>

using namespace Grid;

const RealD tol = 1.0e-10;

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
// Compare one operator application between the two charts
////////////////////////////////////////////////////////////////////////
template<class vOp,class lexOp>
void Compare(std::string name,int dag,
	     vOp &vD,   typename vOp::FermionField   &vsrc,
	     lexOp &lD, typename lexOp::FermionField &lsrc)
{
  typename vOp::FermionField     vres(vsrc.Grid());
  typename lexOp::FermionField   lres(lsrc.Grid());
  typename lexOp::FermionField   vres_lex(lsrc.Grid());
  typename lexOp::FermionField   err(lsrc.Grid());

  if ( dag == DaggerNo ) { vD.M   (vsrc,vres); lD.M   (lsrc,lres); }
  else                   { vD.Mdag(vsrc,vres); lD.Mdag(lsrc,lres); }

  transfer(vres_lex,vres);
  err = vres_lex - lres;
  RealD n = norm2(err)/norm2(lres);
  std::cout << GridLogMessage << name << " simd vs lex: relative " << n
	    << "  |simd|^2 " << norm2(vres) << "  |lex|^2 " << norm2(lres) << std::endl;
  GRID_ASSERT( n < tol );
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls = 8;

  Coordinate latt  = GridDefaultLatt();
  Coordinate vsimd = GridDefaultSimd(Nd,vComplexD::Nsimd());
  Coordinate lsimd({1,1,1,1});
  Coordinate mpi   = GridDefaultMpi();

  GridCartesian         *vUGrid   = new GridCartesian(latt,vsimd,mpi);
  GridRedBlackCartesian *vUrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(vUGrid);
  GridCartesian         *vFGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,vUGrid);
  GridRedBlackCartesian *vFrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,vUGrid);

  GridCartesian         *lUGrid   = new GridCartesian(latt,lsimd,mpi);
  GridRedBlackCartesian *lUrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(lUGrid);
  GridCartesian         *lFGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,lUGrid);
  GridRedBlackCartesian *lFrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,lUGrid);

  std::cout << GridLogMessage << "Ls = " << Ls
	    << "  vectorised Nsimd(4d) = " << vUGrid->Nsimd()
	    << " Nsimd(5d) = " << vFGrid->Nsimd() << std::endl;
  std::cout << GridLogMessage << "     lexicographic Nsimd(4d) = " << lUGrid->Nsimd()
	    << " Nsimd(5d) = " << lFGrid->Nsimd() << std::endl;
  GRID_ASSERT( lUGrid->Nsimd() == 1 );
  GRID_ASSERT( lFGrid->Nsimd() == 1 );

  RealD mass = 0.05;
  RealD M5   = 1.8;

  GridParallelRNG pRNG4(vUGrid); pRNG4.SeedFixedIntegers(std::vector<int>({1,2,3,4}));
  GridParallelRNG pRNG5(vFGrid); pRNG5.SeedFixedIntegers(std::vector<int>({5,6,7,8}));

  //////////////////////////////////////////////////
  // One gauge field and one source, both charts
  //////////////////////////////////////////////////
  LatticeGaugeFieldD Umu(vUGrid); SU<Nc>::HotConfiguration(pRNG4,Umu);
  lexLatticeGaugeFieldD Ulex(lUGrid); transfer(Ulex,Umu);

  typedef DomainWallFermion<WilsonImplD>    vDwf;
  typedef DomainWallFermion<lexWilsonImplD> lexDwf;
  typedef MobiusFermion<WilsonImplD>        vMob;
  typedef MobiusFermion<lexWilsonImplD>     lexMob;

  typename vDwf::FermionField   src(vFGrid);   random(pRNG5,src);
  typename vDwf::FermionField   phi(vFGrid);   random(pRNG5,phi);
  typename lexDwf::FermionField srclex(lFGrid); transfer(srclex,src);
  typename lexDwf::FermionField philex(lFGrid); transfer(philex,phi);

  //////////////////////////////////////////////////
  // Domain wall
  //////////////////////////////////////////////////
  {
    vDwf   Dv (Umu ,*vFGrid,*vFrbGrid,*vUGrid,*vUrbGrid,mass,M5);
    lexDwf Dl (Ulex,*lFGrid,*lFrbGrid,*lUGrid,*lUrbGrid,mass,M5);

    Compare("DWF M   ",DaggerNo ,Dv,src,Dl,srclex);
    Compare("DWF Mdag",DaggerYes,Dv,src,Dl,srclex);

    // Adjoint identity within the lexicographic chart
    typename lexDwf::FermionField Mphi(lFGrid),Mdagchi(lFGrid);
    Dl.M   (srclex,Mphi);
    Dl.Mdag(philex,Mdagchi);
    ComplexD lhs = innerProduct(philex,Mphi);
    ComplexD rhs = innerProduct(Mdagchi,srclex);
    RealD n = abs(lhs-rhs)/abs(lhs);
    std::cout << GridLogMessage << "DWF lex <phi|M src> = " << lhs
	      << "  <Mdag phi|src> = " << rhs << "  relative " << n << std::endl;
    GRID_ASSERT( n < tol );

    // Checkerboarded hopping term
    typename vDwf::FermionField   vo(vFrbGrid),ve(vFrbGrid);
    typename lexDwf::FermionField lo(lFrbGrid),le(lFrbGrid);
    vo.Checkerboard()=Odd; ve.Checkerboard()=Even;
    lo.Checkerboard()=Odd; le.Checkerboard()=Even;
    pickCheckerboard(Odd,vo,src);
    pickCheckerboard(Odd,lo,srclex);
    Dv.Meooe(vo,ve);
    Dl.Meooe(lo,le);
    RealD vn = norm2(ve), ln = norm2(le);
    std::cout << GridLogMessage << "DWF Meooe |simd|^2 " << vn << " |lex|^2 " << ln
	      << " relative " << fabs(vn-ln)/vn << std::endl;
    GRID_ASSERT( fabs(vn-ln)/vn < tol );
  }

  //////////////////////////////////////////////////
  // Mobius, production coefficients
  //////////////////////////////////////////////////
  {
    RealD b=1.5, c=0.5;
    vMob   Dv (Umu ,*vFGrid,*vFrbGrid,*vUGrid,*vUrbGrid,mass,M5,b,c);
    lexMob Dl (Ulex,*lFGrid,*lFrbGrid,*lUGrid,*lUrbGrid,mass,M5,b,c);

    Compare("Mobius M   ",DaggerNo ,Dv,src,Dl,srclex);
    Compare("Mobius Mdag",DaggerYes,Dv,src,Dl,srclex);

    // MooeeInv is the Ls-direction solve; check it inverts Mooee in the lex chart
    typename lexMob::FermionField lo(lFrbGrid),t1(lFrbGrid),t2(lFrbGrid),e(lFrbGrid);
    lo.Checkerboard()=Odd; t1.Checkerboard()=Odd; t2.Checkerboard()=Odd; e.Checkerboard()=Odd;
    pickCheckerboard(Odd,lo,srclex);
    Dl.Mooee(lo,t1);
    Dl.MooeeInv(t1,t2);
    e = t2 - lo;
    RealD n = norm2(e)/norm2(lo);
    std::cout << GridLogMessage << "Mobius lex MooeeInv(Mooee) - 1: relative " << n << std::endl;
    GRID_ASSERT( n < tol );
  }

  std::cout << GridLogMessage << "Test_dwf_lex: ALL PASS" << std::endl;

  Grid_finalize();
}
