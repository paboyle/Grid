    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/approx/MultiShiftFunction.h

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
#ifndef MULTI_SHIFT_FUNCTION
#define MULTI_SHIFT_FUNCTION

namespace Grid {

class MultiShiftFunction {
public:
  int order;
  std::vector<RealD> poles;
  std::vector<RealD> residues;
  std::vector<RealD> tolerances;
  RealD norm;
  RealD lo,hi;

  MultiShiftFunction(int n,RealD _lo,RealD _hi): poles(n), residues(n), lo(_lo), hi(_hi) {;};
  RealD approx(RealD x);
  void csv(std::ostream &out);
  void gnuplot(std::ostream &out);

  void Init(AlgRemez & remez,double tol,bool inverse) 
  {
    order=remez.getDegree();
    tolerances.resize(remez.getDegree(),tol);
    poles.resize(remez.getDegree());
    residues.resize(remez.getDegree());
    remez.getBounds(lo,hi);
    if ( inverse ) remez.getIPFE (&residues[0],&poles[0],&norm);
    else           remez.getPFE (&residues[0],&poles[0],&norm);
  }
  // Allow deferred initialisation
  MultiShiftFunction(void){};
  MultiShiftFunction(AlgRemez & remez,double tol,bool inverse)
  {
    Init(remez,tol,inverse);
  }

};
}
#endif
