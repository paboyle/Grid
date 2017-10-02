/*************************************************************************************

  Grid physics library, www.github.com/paboyle/Grid

  Source file: ./lib/qcd/action/gauge/WilsonGaugeAction.h

  Copyright (C) 2015

  Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
  Author: Peter Boyle <paboyle@ph.ed.ac.uk>
  Author: neo <cossu@post.kek.jp>
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

  See the full license in the file "LICENSE" in the top level distribution
directory
  *************************************************************************************/
/*  END LEGAL */

#ifndef SCALAR_ACTION_H
#define SCALAR_ACTION_H

namespace Grid {
  // FIXME drop the QCD namespace everywhere here

template <class Impl>
class ScalarAction : public QCD::Action<typename Impl::Field> {
 public:
    INHERIT_FIELD_TYPES(Impl);

 private:
    RealD mass_square;
    RealD lambda;

 public:
    ScalarAction(RealD ms, RealD l) : mass_square(ms), lambda(l) {}

    virtual std::string LogParameters() {
      std::stringstream sstream;
      sstream << GridLogMessage << "[ScalarAction] lambda      : " << lambda      << std::endl;
      sstream << GridLogMessage << "[ScalarAction] mass_square : " << mass_square << std::endl;
      return sstream.str();
    }
    virtual std::string action_name() {return "ScalarAction";}

    virtual void refresh(const Field &U, GridParallelRNG &pRNG) {}  // noop as no pseudoferms

    virtual RealD S(const Field &p) {
      return (mass_square * 0.5 + QCD::Nd) * ScalarObs<Impl>::sumphisquared(p) +
    (lambda / 24.) * ScalarObs<Impl>::sumphifourth(p) +
    ScalarObs<Impl>::sumphider(p);
    };

    virtual void deriv(const Field &p,
                       Field &force) {
      Field tmp(p._grid);
      Field p2(p._grid);
      ScalarObs<Impl>::phisquared(p2, p);
      tmp = -(Cshift(p, 0, -1) + Cshift(p, 0, 1));
      for (int mu = 1; mu < QCD::Nd; mu++) tmp -= Cshift(p, mu, -1) + Cshift(p, mu, 1);

      force =+(mass_square + 2. * QCD::Nd) * p + (lambda / 6.) * p2 * p + tmp;
    }
};



}  // namespace Grid

#endif // SCALAR_ACTION_H
