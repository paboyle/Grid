/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/modules/topological_charge.h

Copyright (C) 2017

Author: Guido Cossu <guido.cossu@ed.ac.uk>

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

#ifndef HMC_TOP_CHARGE_H
#define HMC_TOP_CHARGE_H

namespace Grid {
namespace QCD {

// this is only defined for a gauge theory
template <class Impl>
class TopologicalCharge : public HmcObservable<typename Impl::Field> {
 public:
    // here forces the Impl to be of gauge fields
    // if not the compiler will complain
    INHERIT_GIMPL_TYPES(Impl);

    // necessary for HmcObservable compatibility
    typedef typename Impl::Field Field;

    void TrajectoryComplete(int traj,
                            Field &U,
                            GridSerialRNG &sRNG,
                            GridParallelRNG &pRNG) {

    // 4d topological charge
    // Bx = -iF(y,z), By = -iF(z,y), Bz = -iF(x,y)
    GaugeLinkField Bx(U._grid), By(U._grid), Bz(U._grid);
    WilsonLoops<Impl>::FieldStrength(Bx, U, Ydir, Zdir);
    WilsonLoops<Impl>::FieldStrength(By, U, Zdir, Xdir);
    WilsonLoops<Impl>::FieldStrength(Bz, U, Xdir, Ydir);

    // Ex = -iF(t,x), Ey = -iF(t,y), Ez = -iF(t,z)
    GaugeLinkField Ex(U._grid), Ey(U._grid), Ez(U._grid);
    WilsonLoops<Impl>::FieldStrength(Ex, U, Tdir, Xdir);
    WilsonLoops<Impl>::FieldStrength(Ey, U, Tdir, Ydir);
    WilsonLoops<Impl>::FieldStrength(Ez, U, Tdir, Zdir);

    double coeff = 8.0/(32.0*M_PI*M_PI);

    LatticeComplex qfield = coeff*trace(Bx*Ex + By*Ey + Bz*Ez);
    TComplex Tq = sum(qfield);
    Real q = TensorRemove(Tq).real();

    int def_prec = std::cout.precision();

    std::cout << GridLogMessage
        << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
        << "Topological Charge: [ " << traj << " ] "<< q << std::endl;

    std::cout.precision(def_prec);
    }

};
}
}

#endif  //  HMC_TOP_CHARGE_H
