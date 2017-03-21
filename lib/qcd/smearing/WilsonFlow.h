/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/modules/plaquette.h

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

#ifndef WILSONFLOW_H
#define WILSONFLOW_H

namespace Grid {
namespace QCD {

template <class Gimpl>
class WilsonFlow: public Smear<Gimpl>{
    unsigned int Nstep;
    RealD epsilon;

    mutable WilsonGaugeAction<Gimpl> SG;

    void evolve_step(typename Gimpl::GaugeField&) const;
    RealD tau(unsigned int t)const {return epsilon*(t+1.0); }


 public:
    INHERIT_GIMPL_TYPES(Gimpl)

    explicit WilsonFlow(unsigned int Nstep, RealD epsilon):
        Nstep(Nstep),
        epsilon(epsilon),
        SG(WilsonGaugeAction<Gimpl>(3.0)) {
            // WilsonGaugeAction with beta 3.0
            assert(epsilon > 0.0);
            LogMessage();
    }

    void LogMessage() {
        std::cout << GridLogMessage
            << "[WilsonFlow] Nstep   : " << Nstep << std::endl;
        std::cout << GridLogMessage
            << "[WilsonFlow] epsilon : " << epsilon << std::endl;
        std::cout << GridLogMessage
            << "[WilsonFlow] full trajectory : " << Nstep * epsilon << std::endl;
    }

    virtual void smear(GaugeField&, const GaugeField&) const;

    virtual void derivative(GaugeField&, const GaugeField&, const GaugeField&) const {
        assert(0);
        // undefined for WilsonFlow
    }

    RealD energyDensityPlaquette(unsigned int step, const GaugeField& U) const;
};


////////////////////////////////////////////////////////////////////////////////
// Implementations
////////////////////////////////////////////////////////////////////////////////
template <class Gimpl>
void WilsonFlow<Gimpl>::evolve_step(typename Gimpl::GaugeField &U) const{
    GaugeField Z(U._grid);
    GaugeField tmp(U._grid);
    SG.deriv(U, Z);
    Z *= 0.25;                                  // Z0 = 1/4 * F(U)
    Gimpl::update_field(Z, U, -2.0*epsilon);    // U = W1 = exp(ep*Z0)*W0

    Z *= -17.0/8.0;
    SG.deriv(U, tmp); Z += tmp;                 // -17/32*Z0 +Z1
    Z *= 8.0/9.0;                               // Z = -17/36*Z0 +8/9*Z1
    Gimpl::update_field(Z, U, -2.0*epsilon);    // U_= W2 = exp(ep*Z)*W1

    Z *= -4.0/3.0;
    SG.deriv(U, tmp); Z += tmp;                 // 4/3*(17/36*Z0 -8/9*Z1) +Z2
    Z *= 3.0/4.0;                               // Z = 17/36*Z0 -8/9*Z1 +3/4*Z2
    Gimpl::update_field(Z, U, -2.0*epsilon);    // V(t+e) = exp(ep*Z)*W2
}

template <class Gimpl>
RealD WilsonFlow<Gimpl>::energyDensityPlaquette(unsigned int step, const GaugeField& U) const {
    RealD td = tau(step);
    return 2.0 * td * td * SG.S(U)/U._grid->gSites();
}

template <class Gimpl>
void WilsonFlow<Gimpl>::smear(GaugeField& out, const GaugeField& in) const {
    out = in;
    for (unsigned int step = 0; step < Nstep; step++) {
        evolve_step(out);
        std::cout << "[WilsonFlow] Energy density (plaq) : "
            << step << " "
            << energyDensityPlaquette(step,out) << std::endl;
    }

}

}  // namespace QCD
}  // namespace Grid

#endif   // WILSONFLOW_H
