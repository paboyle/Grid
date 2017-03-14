/*************************************************************************************

  Grid physics library, www.github.com/paboyle/Grid

  Source file: ./lib/qcd/action/gauge/WilsonGaugeAction.h

  Copyright (C) 2015

  Author: Guido Cossu <guido,cossu@ed.ac.uk>

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

#ifndef SCALAR_INT_ACTION_H
#define SCALAR_INT_ACTION_H

namespace Grid {
  // FIXME drop the QCD namespace everywhere here

template <class Impl>
class ScalarInteractionAction : public QCD::Action<typename Impl::Field> {
    RealD mass_square;
    RealD lambda;

 public:
    INHERIT_FIELD_TYPES(Impl);
    ScalarInteractionAction(RealD ms, RealD l) : mass_square(ms), lambda(l) {}

    virtual std::string LogParameters() {
      std::stringstream sstream;
      sstream << GridLogMessage << "[ScalarAction] lambda      : " << lambda      << std::endl;
      sstream << GridLogMessage << "[ScalarAction] mass_square : " << mass_square << std::endl;
      return sstream.str();
    }

    virtual std::string action_name() {return "ScalarAction";}

    virtual void refresh(const Field &U,
                         GridParallelRNG &pRNG) {}  // noop as no pseudoferms

    virtual RealD S(const Field &p) {
        Field action(p._grid);
        Field pshift(p._grid);
        Field phisquared(p._grid);
        phisquared = p*p;
        action = (2.0*QCD::Nd + mass_square)*phisquared + lambda*phisquared*phisquared;
        for (int mu = 0; mu < QCD::Nd; mu++) {
            pshift = Cshift(p, mu, +1);  // not efficient implement with stencils
            action -= pshift*p + p*pshift;
        }
        return -(TensorRemove(sum(trace(action)))).real();
    };

    virtual void deriv(const Field &p,
                       Field &force) {
        force = (2.0*QCD::Nd + mass_square)*p + 2.0*lambda*p*p*p;
        // following is inefficient
        for (int mu = 0; mu < QCD::Nd; mu++) force -= Cshift(p, mu, -1) + Cshift(p, mu, 1);
    }
};

}  // namespace Grid

#endif  // SCALAR_INT_ACTION_H
