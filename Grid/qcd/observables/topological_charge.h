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

struct TopologySmearingParameters : Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(TopologySmearingParameters,
    int, steps,
    float, step_size,
    int, meas_interval,
    float, maxTau);

    TopologySmearingParameters(int s = 0, float ss = 0.0f, int mi = 0, float mT = 0.0f):
        steps(s), step_size(ss), meas_interval(mi), maxTau(mT){}

    template < class ReaderClass >
    TopologySmearingParameters(Reader<ReaderClass>& Reader){
        read(Reader, "Smearing", *this);  
    }  
};



struct TopologyObsParameters : Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(TopologyObsParameters,
      int, interval,
      bool, do_smearing,
      TopologySmearingParameters, Smearing);  

    TopologyObsParameters(int interval = 1, bool smearing = false):
        interval(interval), Smearing(smearing){}

    template <class ReaderClass >
      TopologyObsParameters(Reader<ReaderClass>& Reader){
        read(Reader, "TopologyMeasurement", *this);
  }
};


// this is only defined for a gauge theory
template <class Impl>
class TopologicalCharge : public HmcObservable<typename Impl::Field> {
    TopologyObsParameters Pars;

 public:
    // here forces the Impl to be of gauge fields
    // if not the compiler will complain
    INHERIT_GIMPL_TYPES(Impl);

    // necessary for HmcObservable compatibility
    typedef typename Impl::Field Field;

    TopologicalCharge(int interval = 1, bool do_smearing = false):
        Pars(interval, do_smearing){}
    
    TopologicalCharge(TopologyObsParameters P):Pars(P){
        std::cout << GridLogDebug << "Creating TopologicalCharge " << std::endl;
    }

    void TrajectoryComplete(int traj,
                            Field &U,
                            GridSerialRNG &sRNG,
                            GridParallelRNG &pRNG) {

    if (traj%Pars.interval == 0){
        // Smearing
        Field Usmear = U;
        int def_prec = std::cout.precision();
        
        if (Pars.do_smearing){
            // using wilson flow by default here
            WilsonFlow<PeriodicGimplR> WF(Pars.Smearing.steps, Pars.Smearing.step_size, Pars.Smearing.meas_interval);
            WF.smear_adaptive(Usmear, U, Pars.Smearing.maxTau);
            Real T0   = WF.energyDensityPlaquette(Usmear);
            std::cout << GridLogMessage << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
                      << "T0                : [ " << traj << " ] "<< T0 << std::endl;
        }

        Real q    = WilsonLoops<Impl>::TopologicalCharge(Usmear);
        std::cout << GridLogMessage
            << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
            << "Topological Charge: [ " << traj << " ] "<< q << std::endl;

        std::cout.precision(def_prec);
        }
    }

};
}
}

#endif  //  HMC_TOP_CHARGE_H
