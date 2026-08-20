/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/lexLattice/Test_RectPlaq_lex.cc

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
// Plaquette and 2x1 rectangle built from covariant shifts, checked against
// WilsonLoops, plus link trace and a blockSum coarsening.
//
// The measurement is written once against a gauge implementation and run in
// both the vectorised and the lexicographic chart; every number must agree.
// The configuration crosses layouts through NERSC IO, whose reader verifies
// the header written by the other chart.
//
#include <Grid/Grid.h>

using namespace Grid;

const RealD tol = 1.0e-10;

template<class Gimpl>
struct Measurements
{
  RealD plaq_shift;
  RealD plaq_loops;
  RealD plaq_staple;
  RealD rect_shift;
  RealD rect_loops;
  RealD link;
  RealD coarse_plaq;
};

template<class Gimpl>
void Measure(typename Gimpl::Field &Umu,GridBase *coarse,Measurements<Gimpl> &m)
{
  typedef typename Gimpl::LinkField    LinkField;
  typedef typename Gimpl::ComplexField ComplexField;
  typedef WilsonLoops<Gimpl>           WL;

  GridBase *grid = Umu.Grid();
  RealD vol = grid->gSites();

  std::vector<LinkField> U(Nd,grid);
  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }

  ///////////////////////////////////////////////////
  // Link trace
  ///////////////////////////////////////////////////
  ComplexField LinkTrace(grid);
  LinkTrace = Zero();
  for(int mu=0;mu<Nd;mu++){
    LinkTrace = LinkTrace + trace(U[mu]);
  }
  m.link = real(TensorRemove(sum(LinkTrace)))/vol/Nd/Nc;

  ///////////////////////////////////////////////////
  // Plaquette from covariant shifts
  ///////////////////////////////////////////////////
  ComplexField Plaq(grid);
  Plaq = Zero();
  for(int mu=1;mu<Nd;mu++){
  for(int nu=0;nu<mu;nu++){
    Plaq = Plaq + trace(PeriodicBC::CovShiftForward(U[mu],mu,U[nu])
		    *adj(PeriodicBC::CovShiftForward(U[nu],nu,U[mu])));
  }}
  m.plaq_shift = real(TensorRemove(sum(Plaq)))/vol/6.0/Nc;
  m.plaq_loops = WL::avgPlaquette(Umu);

  ///////////////////////////////////////////////////
  // 2x1 and 1x2 rectangles from covariant shifts
  ///////////////////////////////////////////////////
  ComplexField Rect(grid);
  Rect = Zero();
  for(int mu=1;mu<Nd;mu++){
  for(int nu=0;nu<mu;nu++){
    Rect = Rect + trace(
	   PeriodicBC::CovShiftForward(U[mu],mu,PeriodicBC::CovShiftForward(U[mu],mu,U[nu]))
       *adj(PeriodicBC::CovShiftForward(U[nu],nu,PeriodicBC::CovShiftForward(U[mu],mu,U[mu]))) );
    Rect = Rect + trace(
	   PeriodicBC::CovShiftForward(U[mu],mu,PeriodicBC::CovShiftForward(U[nu],nu,U[nu]))
       *adj(PeriodicBC::CovShiftForward(U[nu],nu,PeriodicBC::CovShiftForward(U[nu],nu,U[mu]))) );
  }}
  m.rect_shift = real(TensorRemove(sum(Rect)))/vol/12.0/Nc;
  m.rect_loops = WL::avgRectangle(Umu);

  ///////////////////////////////////////////////////
  // Plaquette through the staples
  ///////////////////////////////////////////////////
  {
    RealD stap = 0.0;
    LinkField    staple(grid);
    ComplexField stap_tr(grid);
    for(int mu=0;mu<Nd;mu++){
      WL::Staple(staple,Umu,mu);
      stap_tr = trace(U[mu]*staple);
      stap += real(TensorRemove(sum(stap_tr)));
    }
    m.plaq_staple = stap/vol/6.0/Nc/4.0;
  }

  ///////////////////////////////////////////////////
  // Coarsened plaquette; blockSum must conserve the sum
  ///////////////////////////////////////////////////
  {
    ComplexField cPlaq(coarse);
    blockSum(cPlaq,Plaq);
    m.coarse_plaq = real(TensorRemove(sum(cPlaq)))/vol/6.0/Nc;
  }
}

template<class Gimpl>
void Report(std::string name,Measurements<Gimpl> &m)
{
  std::cout << GridLogMessage << name << ": plaquette shifts " << m.plaq_shift
            << " loops " << m.plaq_loops << " staples " << m.plaq_staple << std::endl;
  std::cout << GridLogMessage << name << ": rectangle shifts " << m.rect_shift
            << " loops " << m.rect_loops << std::endl;
  std::cout << GridLogMessage << name << ": link trace " << m.link
            << " coarsened plaquette " << m.coarse_plaq << std::endl;
}

template<class Gimpl>
void SelfConsistent(Measurements<Gimpl> &m)
{
  GRID_ASSERT( fabs(m.plaq_shift  - m.plaq_loops) < 1.0e-8 );
  GRID_ASSERT( fabs(m.plaq_staple - m.plaq_loops) < 1.0e-8 );
  GRID_ASSERT( fabs(m.rect_shift  - m.rect_loops) < 1.0e-8 );
  GRID_ASSERT( fabs(m.coarse_plaq - m.plaq_shift) < 1.0e-8 );
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Coordinate latt  = GridDefaultLatt();
  Coordinate vsimd = GridDefaultSimd(Nd,vComplexD::Nsimd());
  Coordinate lsimd({1,1,1,1});
  Coordinate mpi   = GridDefaultMpi();

  Coordinate clatt(Nd);
  for(int d=0;d<Nd;d++){
    GRID_ASSERT( (latt[d]%2) == 0 );
    clatt[d] = latt[d]/2;
  }

  GridCartesian vGrid(latt,vsimd,mpi);
  GridCartesian lGrid(latt,lsimd,mpi);
  GridCartesian vCoarse(clatt,vsimd,mpi);
  GridCartesian lCoarse(clatt,lsimd,mpi);

  std::cout << GridLogMessage << "vectorised    Nsimd = " << vGrid.Nsimd() << std::endl;
  std::cout << GridLogMessage << "lexicographic Nsimd = " << lGrid.Nsimd() << std::endl;
  GRID_ASSERT( lGrid.Nsimd() == 1 );

  GridParallelRNG pRNG(&vGrid);
  pRNG.SeedFixedIntegers(std::vector<int>({1,2,3,4}));

  LatticeGaugeFieldD Umu(&vGrid);
  SU<Nc>::HotConfiguration(pRNG,Umu);

  std::string file("./ckpoint_rectplaq_lex.4000");
  NerscIO::writeConfiguration(Umu,file,0,0);

  lexLatticeGaugeFieldD Ulex(&lGrid);
  FieldMetaData header;
  NerscIO::readConfiguration(Ulex,header,file);

  Measurements<PeriodicGimplD>    vm;
  Measurements<lexPeriodicGimplD> lm;

  Measure<PeriodicGimplD>   (Umu ,&vCoarse,vm);
  Measure<lexPeriodicGimplD>(Ulex,&lCoarse,lm);

  Report("simd",vm);
  Report("lex ",lm);

  SelfConsistent(vm);
  SelfConsistent(lm);

  GRID_ASSERT( fabs(vm.plaq_shift  - lm.plaq_shift ) < tol );
  GRID_ASSERT( fabs(vm.plaq_loops  - lm.plaq_loops ) < tol );
  GRID_ASSERT( fabs(vm.plaq_staple - lm.plaq_staple) < tol );
  GRID_ASSERT( fabs(vm.rect_shift  - lm.rect_shift ) < tol );
  GRID_ASSERT( fabs(vm.rect_loops  - lm.rect_loops ) < tol );
  GRID_ASSERT( fabs(vm.link        - lm.link       ) < tol );
  GRID_ASSERT( fabs(vm.coarse_plaq - lm.coarse_plaq) < tol );

  std::cout << GridLogMessage << "Test_RectPlaq_lex: ALL PASS" << std::endl;

  Grid_finalize();
}
