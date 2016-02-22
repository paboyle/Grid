    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/WilsonTMFermion.cc

    Copyright (C) 2015

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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#include <Grid.h>

namespace Grid {
namespace QCD {

    /*
     * BF sequence
     *
      void bfmbase<Float>::MooeeInv(Fermion_t psi, 
			       Fermion_t chi, 
			      int dag, int cb)

    double m    = this->mass;
    double tm   = this->twistedmass;
    double mtil = 4.0+this->mass;

    double sq = mtil*mtil + tm*tm;

    double a = mtil/sq;
    double b = -tm /sq;
    if(dag) b=-b;
    axpibg5x(chi,psi,a,b);

      void bfmbase<Float>::Mooee(Fermion_t psi, 
			   Fermion_t chi, 
			   int dag,int cb)
    double a = 4.0+this->mass;
    double b = this->twistedmass;
    if(dag) b=-b;
    axpibg5x(chi,psi,a,b);
    */

  template<class Impl>
  void WilsonTMFermion<Impl>::Mooee(const FermionField &in, FermionField &out) {
    RealD a = 4.0+this->mass;
    RealD b = this->mu;
    out.checkerboard = in.checkerboard;
    axpibg5x(out,in,a,b);
  }
  template<class Impl>
  void WilsonTMFermion<Impl>::MooeeDag(const FermionField &in, FermionField &out) {
    RealD a = 4.0+this->mass;
    RealD b = -this->mu;
    out.checkerboard = in.checkerboard;
    axpibg5x(out,in,a,b);
  }
  template<class Impl>
  void WilsonTMFermion<Impl>::MooeeInv(const FermionField &in, FermionField &out) {
    RealD m    = this->mass;
    RealD tm   = this->mu;
    RealD mtil = 4.0+this->mass;
    RealD sq   = mtil*mtil+tm*tm;
    RealD a    = mtil/sq;
    RealD b    = -tm /sq;
    axpibg5x(out,in,a,b);
  }
  template<class Impl>
  void WilsonTMFermion<Impl>::MooeeInvDag(const FermionField &in, FermionField &out) {
    RealD m    = this->mass;
    RealD tm   = this->mu;
    RealD mtil = 4.0+this->mass;
    RealD sq   = mtil*mtil+tm*tm;
    RealD a    = mtil/sq;
    RealD b    = tm /sq;
    axpibg5x(out,in,a,b);
  }

  FermOpTemplateInstantiate(WilsonTMFermion);

}
}
