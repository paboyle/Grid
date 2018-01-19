/*************************************************************************************

  Grid physics library, www.github.com/paboyle/Grid

  Source file: ./lib/qcd/action/gauge/shGordonAction.h

  Copyright (C) 2018

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

#ifndef SHGORDON_ACTION_H
#define SHGORDON_ACTION_H

namespace Grid {

template <class Impl>
class shGordonAction : public QCD::Action<typename Impl::Field> {
 public:
    INHERIT_FIELD_TYPES(Impl);

 private:
    RealD mass_square;
    RealD g;

 public:
    shGordonAction(RealD ms, RealD g) : mass_square(ms), g(g) {}

    virtual std::string LogParameters() {
      std::stringstream sstream;
      sstream << GridLogMessage << "[shGordonAction] g           : " << g           << std::endl;
      sstream << GridLogMessage << "[shGordonAction] mass_square : " << mass_square << std::endl;
      return sstream.str();
    }
    virtual std::string action_name() {return "shGordonAction";}

    virtual void refresh(const Field &U, GridParallelRNG &pRNG) {}  // noop as no pseudoferms

    virtual RealD S(const Field &phi) {
      return QCD::Nd * ScalarObs<Impl>::sumphisquared(phi) + ScalarObs<Impl>::sumphider(phi) + 0.5*mass_square/(g*g)*sum(trace(exp(g*phi) + exp(-g*phi)))   ;
    };

    virtual void deriv(const Field &phi,
                       Field &force) {
        Field tmp(phi._grid);
        tmp = 2.0*QCD::Nd*phi;
        for (int mu = 0; mu < QCD::Nd; mu++) tmp -= Cshift(phi, mu, 1) + Cshift(phi, mu, -1);

        std::cout << GridLogDebug << "Phi norm : " << norm2(phi) << std::endl;
        force+= tmp + 0.5*mass_square/g*(exp(g*phi) - exp(-g*phi));
    }
};



}  // namespace Grid

#endif // SHGORDON_ACTION_H
