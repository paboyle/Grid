    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_remez.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::cout<<GridLogMessage << "Testing Remez"<<std::endl;

  double     lo=1.0e-3;
  double     hi=5.0;
  int precision=64;
  int    degree=16;
  AlgRemez remez(lo,hi,precision);

  ////////////////////////////////////////
  // sqrt and inverse sqrt
  ////////////////////////////////////////

  std::cout<<GridLogMessage << "Generating degree "<<degree<<" for x^(1/2)"<<std::endl;
  remez.generateApprox(degree,1,2);
  MultiShiftFunction Root2(remez,1.0,false);
  MultiShiftFunction InvRoot2(remez,1.0,true);


  std::cout<<GridLogMessage << "Generating degree "<<degree<<" for x^(1/4)"<<std::endl;
  remez.generateApprox(degree,1,4);
  MultiShiftFunction Root4(remez,1.0,false);
  MultiShiftFunction InvRoot4(remez,1.0,true);

  std::cout<<GridLogMessage << "Generating degree "<<degree<<" for x^(1/8)"<<std::endl;
  remez.generateApprox(degree,1,8);
  MultiShiftFunction Root8(remez,1.0,false);
  MultiShiftFunction InvRoot8(remez,1.0,true);

  std::cout<<GridLogMessage << "Generating degree "<<degree<<" for x^(1/16)"<<std::endl;
  remez.generateApprox(degree,1,16);
  MultiShiftFunction Root16(remez,1.0,false);
  MultiShiftFunction InvRoot16(remez,1.0,true);

  std::cout<<GridLogMessage << "Generating degree "<<degree<<" for x^(1/32)"<<std::endl;
  remez.generateApprox(degree,1,32);
  MultiShiftFunction Root32(remez,1.0,false);
  MultiShiftFunction InvRoot32(remez,1.0,true);
  
  ofstream gnuplot(std::string("Root2.gnu"),std::ios::out|std::ios::trunc);
  Root2.gnuplot(gnuplot);

  ofstream gnuplot_i2(std::string("InvRoot2.gnu"),std::ios::out|std::ios::trunc);
  InvRoot2.gnuplot(gnuplot_i2);

  ofstream gnuplot_i4(std::string("InvRoot4.gnu"),std::ios::out|std::ios::trunc);
  InvRoot4.gnuplot(gnuplot_i4);

  ofstream gnuplot_i8(std::string("InvRoot8.gnu"),std::ios::out|std::ios::trunc);
  InvRoot8.gnuplot(gnuplot_i8);

  ofstream gnuplot_i16(std::string("InvRoot16.gnu"),std::ios::out|std::ios::trunc);
  InvRoot16.gnuplot(gnuplot_i16);

  ofstream gnuplot_i32(std::string("InvRoot32.gnu"),std::ios::out|std::ios::trunc);
  InvRoot32.gnuplot(gnuplot_i32);




  double x=0.6789;
  double sx=std::sqrt(x);
  double ssx=std::sqrt(sx);
  double isx=1.0/sx;
  double issx=1.0/ssx;

  double asx  =Root2.approx(x);
  double assx =Root4.approx(x);
  double aisx =InvRoot2.approx(x);
  double aissx=InvRoot4.approx(x);

  std::cout<<GridLogMessage << "x^(1/2) : "<<sx<<" "<<asx<<std::endl;
  std::cout<<GridLogMessage << "x^(1/4) : "<<ssx<<" "<<assx<<std::endl;
  std::cout<<GridLogMessage << "x^(-1/2): "<<isx<<" "<<aisx<<std::endl;
  std::cout<<GridLogMessage << "x^(-1/4): "<<issx<<" "<<aissx<<std::endl;
  assert(fabs(sx-asx)<1.0e-6);
  assert(fabs(ssx-assx)<1.0e-6);
  assert(fabs(isx-aisx)<1.0e-6);
  assert(fabs(issx-aissx)<1.0e-6);

  Grid_finalize();
}
