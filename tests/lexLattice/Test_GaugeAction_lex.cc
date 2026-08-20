/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/lexLattice/Test_GaugeAction_lex.cc

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
// Wilson loops on a lexicographic gauge field.
//
// The configuration is generated on the vectorised grid, written as NERSC,
// and read back into the lexicographic layout: the reader recomputes the
// plaquette and link trace and checks them against the header written by
// the vectorised chart, so the crossing is validated inside Grid's own IO.
// Plaquette, link trace, rectangle and the staple identity are then compared
// between the two charts.
//
#include <Grid/Grid.h>

using namespace Grid;

const RealD tol = 1.0e-10;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  typedef WilsonLoops<PeriodicGimplD>    vWL;
  typedef WilsonLoops<lexPeriodicGimplD> lexWL;

  Coordinate latt  = GridDefaultLatt();
  Coordinate vsimd = GridDefaultSimd(Nd,vComplexD::Nsimd());
  Coordinate lsimd({1,1,1,1});
  Coordinate mpi   = GridDefaultMpi();

  GridCartesian vGrid(latt,vsimd,mpi);
  GridCartesian lGrid(latt,lsimd,mpi);

  std::cout << GridLogMessage << "vectorised    Nsimd = " << vGrid.Nsimd() << std::endl;
  std::cout << GridLogMessage << "lexicographic Nsimd = " << lGrid.Nsimd() << std::endl;
  GRID_ASSERT( lGrid.Nsimd() == 1 );

  GridParallelRNG pRNG(&vGrid);
  pRNG.SeedFixedIntegers(std::vector<int>({1,2,3,4}));

  LatticeGaugeFieldD Umu(&vGrid);
  SU<Nc>::HotConfiguration(pRNG,Umu);

  std::string file("./ckpoint_lex.4000");
  NerscIO::writeConfiguration(Umu,file,0,0);

  lexLatticeGaugeFieldD Ulex(&lGrid);
  FieldMetaData header;
  NerscIO::readConfiguration(Ulex,header,file);

  //////////////////////////////////////////////////
  // Plaquette, link trace, rectangle in both charts
  //////////////////////////////////////////////////
  RealD vplaq = vWL::avgPlaquette(Umu);
  RealD lplaq = lexWL::avgPlaquette(Ulex);

  RealD vlink = vWL::linkTrace(Umu);
  RealD llink = lexWL::linkTrace(Ulex);

  RealD vrect = vWL::avgRectangle(Umu);
  RealD lrect = lexWL::avgRectangle(Ulex);

  std::cout << GridLogMessage << "plaquette   simd " << vplaq << " lex " << lplaq
            << " header " << header.plaquette << std::endl;
  std::cout << GridLogMessage << "link trace  simd " << vlink << " lex " << llink
            << " header " << header.link_trace << std::endl;
  std::cout << GridLogMessage << "rectangle   simd " << vrect << " lex " << lrect << std::endl;

  GRID_ASSERT( fabs(vplaq-lplaq) < tol );
  GRID_ASSERT( fabs(vlink-llink) < tol );
  GRID_ASSERT( fabs(vrect-lrect) < tol );

  //////////////////////////////////////////////////
  // Plaquette via staples, in the lexicographic chart
  //////////////////////////////////////////////////
  {
    RealD vol = lGrid.gSites();
    RealD stap_plaq = 0.0;

    lexLatticeColourMatrixD stap(&lGrid);
    lexLatticeColourMatrixD Ul(&lGrid);
    lexLatticeComplexD      stap_tr(&lGrid);

    for(int mu=0;mu<Nd;mu++){
      Ul = PeekIndex<LorentzIndex>(Ulex,mu);
      lexWL::Staple(stap,Ulex,mu);
      stap_tr = trace(Ul*stap);
      auto Ts = sum(stap_tr);
      stap_plaq += real(TensorRemove(Ts));
    }
    RealD StapScale = 1.0/vol/6.0/Nc/4.0;
    RealD plaq_from_staples = stap_plaq*StapScale;

    std::cout << GridLogMessage << "plaquette via staples (lex) " << plaq_from_staples
              << " direct " << lplaq << std::endl;
    GRID_ASSERT( fabs(plaq_from_staples - lplaq) < 1.0e-8 );
  }

  std::cout << GridLogMessage << "Test_GaugeAction_lex: ALL PASS" << std::endl;

  Grid_finalize();
}
