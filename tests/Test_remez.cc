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
#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::cout<<GridLogMessage << "Testing Remez"<<std::endl;

  double     lo=0.01;
  double     hi=1.0;
  int precision=64;
  int    degree=10;
  AlgRemez remez(0.001,1.0,precision);

  ////////////////////////////////////////
  // sqrt and inverse sqrt
  ////////////////////////////////////////

  std::cout<<GridLogMessage << "Generating degree "<<degree<<" for x^(1/2)"<<std::endl;
  remez.generateApprox(degree,1,2);
  MultiShiftFunction Sqrt(remez,1.0,false);
  MultiShiftFunction InvSqrt(remez,1.0,true);


  std::cout<<GridLogMessage << "Generating degree "<<degree<<" for x^(1/4)"<<std::endl;
  remez.generateApprox(degree,1,4);
  MultiShiftFunction SqrtSqrt(remez,1.0,false);
  MultiShiftFunction InvSqrtSqrt(remez,1.0,true);

  
  ofstream gnuplot(std::string("Sqrt.gnu"),std::ios::out|std::ios::trunc);
  Sqrt.gnuplot(gnuplot);

  ofstream gnuplot_inv(std::string("InvSqrt.gnu"),std::ios::out|std::ios::trunc);
  InvSqrt.gnuplot(gnuplot);

  double x=0.6789;
  double sx=std::sqrt(x);
  double ssx=std::sqrt(sx);
  double isx=1.0/sx;
  double issx=1.0/ssx;

  double asx  =Sqrt.approx(x);
  double assx =SqrtSqrt.approx(x);
  double aisx =InvSqrt.approx(x);
  double aissx=InvSqrtSqrt.approx(x);

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
