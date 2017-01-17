/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/hmc/HMC.h

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
//--------------------------------------------------------------------
/*! @file HMC.h
 * @brief Classes for Hybrid Monte Carlo update
 *
 * @author Guido Cossu
 */
//--------------------------------------------------------------------
#ifndef HMC_INCLUDED
#define HMC_INCLUDED

#include <string>

namespace Grid {
namespace QCD {

struct HMCparameters: Serializable {
	GRID_SERIALIZABLE_CLASS_MEMBERS(HMCparameters,
  Integer, StartTrajectory,
  Integer, Trajectories, /* @brief Number of sweeps in this run */
  bool, MetropolisTest,
  Integer, NoMetropolisUntil,
  )

  // nest here the MDparameters and make all serializable

  HMCparameters() {
    ////////////////////////////// Default values
    MetropolisTest    = true;
    NoMetropolisUntil = 10;
    StartTrajectory   = 0;
    Trajectories      = 10;
    /////////////////////////////////
  }

  void print_parameters() const {
    std::cout << GridLogMessage << "[HMC parameters] Trajectories            : " << Trajectories << "\n";
    std::cout << GridLogMessage << "[HMC parameters] Start trajectory        : " << StartTrajectory << "\n";
    std::cout << GridLogMessage << "[HMC parameters] Metropolis test (on/off): " << std::boolalpha << MetropolisTest << "\n";
    std::cout << GridLogMessage << "[HMC parameters] Thermalization trajs    : " << NoMetropolisUntil << "\n";
  }
  
};

template <class Field>
class HmcObservable {
 public:
  virtual void TrajectoryComplete(int traj, Field &U, GridSerialRNG &sRNG,
                                  GridParallelRNG &pRNG) = 0;
};

  // this is only defined for a gauge theory
template <class Gimpl>
class PlaquetteLogger : public HmcObservable<typename Gimpl::Field> {
 private:
  std::string Stem;

 public:
  INHERIT_GIMPL_TYPES(Gimpl);
  PlaquetteLogger(std::string cf) { Stem = cf; };

  void TrajectoryComplete(int traj, GaugeField &U, GridSerialRNG &sRNG,
                          GridParallelRNG &pRNG) {
    std::string file;
    {
      std::ostringstream os;
      os << Stem << "." << traj;
      file = os.str();
    }
    std::ofstream of(file);

    RealD peri_plaq = WilsonLoops<PeriodicGimplR>::avgPlaquette(U);
    RealD peri_rect = WilsonLoops<PeriodicGimplR>::avgRectangle(U);

    RealD impl_plaq = WilsonLoops<Gimpl>::avgPlaquette(U);
    RealD impl_rect = WilsonLoops<Gimpl>::avgRectangle(U);

    of << traj << " " << impl_plaq << " " << impl_rect << "  " << peri_plaq
       << " " << peri_rect << std::endl;
    std::cout << GridLogMessage << "traj"
              << " "
              << "plaq "
              << " "
              << " rect  "
              << "  "
              << "peri_plaq"
              << " "
              << "peri_rect" << std::endl;
    std::cout << GridLogMessage << traj << " " << impl_plaq << " " << impl_rect
              << "  " << peri_plaq << " " << peri_rect << std::endl;
  }
};

template <class IntegratorType>
class HybridMonteCarlo {
 private:
  const HMCparameters Params;

  typedef typename IntegratorType::Field Field;
  
  GridSerialRNG &sRNG;    // Fixme: need a RNG management strategy.
  GridParallelRNG &pRNG;  // Fixme: need a RNG management strategy.
  Field &Ucur;

  IntegratorType &TheIntegrator;
  std::vector<HmcObservable<Field> *> Observables;

  /////////////////////////////////////////////////////////
  // Metropolis step
  /////////////////////////////////////////////////////////
  bool metropolis_test(const RealD DeltaH) {
    RealD rn_test;

    RealD prob = std::exp(-DeltaH);

    random(sRNG, rn_test);

    std::cout << GridLogMessage
              << "--------------------------------------------------\n";
    std::cout << GridLogMessage << "exp(-dH) = " << prob
              << "  Random = " << rn_test << "\n";
    std::cout << GridLogMessage
              << "Acc. Probability = " << ((prob < 1.0) ? prob : 1.0) << "\n";

    if ((prob > 1.0) || (rn_test <= prob)) {  // accepted
      std::cout << GridLogMessage << "Metropolis_test -- ACCEPTED\n";
      std::cout << GridLogMessage
                << "--------------------------------------------------\n";
      return true;
    } else {  // rejected
      std::cout << GridLogMessage << "Metropolis_test -- REJECTED\n";
      std::cout << GridLogMessage
                << "--------------------------------------------------\n";
      return false;
    }
  }

  /////////////////////////////////////////////////////////
  // Evolution
  /////////////////////////////////////////////////////////
  RealD evolve_step(Field &U) {
    TheIntegrator.refresh(U, pRNG);  // set U and initialize P and phi's

    RealD H0 = TheIntegrator.S(U);  // initial state action

    std::streamsize current_precision = std::cout.precision();
    std::cout.precision(17);
    std::cout << GridLogMessage << "Total H before trajectory = " << H0 << "\n";
    std::cout.precision(current_precision);

    TheIntegrator.integrate(U);

    RealD H1 = TheIntegrator.S(U);  // updated state action

    std::cout.precision(17);
    std::cout << GridLogMessage << "Total H after trajectory  = " << H1
              << "  dH = " << H1 - H0 << "\n";
    std::cout.precision(current_precision);

    return (H1 - H0);
  }

 public:
  /////////////////////////////////////////
  // Constructor
  /////////////////////////////////////////
  HybridMonteCarlo(HMCparameters Pams, IntegratorType &_Int,
                   GridSerialRNG &_sRNG, GridParallelRNG &_pRNG, Field &_U)
      : Params(Pams), TheIntegrator(_Int), sRNG(_sRNG), pRNG(_pRNG), Ucur(_U) {}
  ~HybridMonteCarlo(){};

  void AddObservable(HmcObservable<Field> *obs) {
    Observables.push_back(obs);
  }

  void evolve(void) {
    Real DeltaH;

    Field Ucopy(Ucur._grid);

    Params.print_parameters();
    TheIntegrator.print_parameters();
    TheIntegrator.print_actions();

    // Actual updates (evolve a copy Ucopy then copy back eventually)
    unsigned int FinalTrajectory = Params.Trajectories + Params.NoMetropolisUntil + Params.StartTrajectory;
    for (int traj = Params.StartTrajectory; traj < FinalTrajectory; ++traj) {
      std::cout << GridLogMessage << "-- # Trajectory = " << traj << "\n";
      if (traj < Params.StartTrajectory + Params.NoMetropolisUntil) {
      	std::cout << GridLogMessage << "-- Thermalization" << std::endl;
    	}

    	double t0=usecond();
      Ucopy = Ucur;

      DeltaH = evolve_step(Ucopy);

      bool accept = true;
      if (traj >= Params.StartTrajectory + Params.NoMetropolisUntil) {
        accept = metropolis_test(DeltaH);
      } else {
      	std::cout << GridLogMessage << "Skipping Metropolis test" << std::endl;
      }

      if (accept) {
        Ucur = Ucopy;
      }
      double t1=usecond();
      std::cout << GridLogMessage << "Total time for trajectory (s): " << (t1-t0)/1e6 << std::endl;


      for (int obs = 0; obs < Observables.size(); obs++) {
      	std::cout << GridLogDebug << "Observables # " << obs << std::endl;
      	std::cout << GridLogDebug << "Observables total " << Observables.size() << std::endl;
      	std::cout << GridLogDebug << "Observables pointer " << Observables[obs] << std::endl;
        Observables[obs]->TrajectoryComplete(traj + 1, Ucur, sRNG, pRNG);
      }
      std::cout << GridLogMessage << ":::::::::::::::::::::::::::::::::::::::::::" << std::endl;
    }
  }
};

}  // QCD
}  // Grid

#endif
