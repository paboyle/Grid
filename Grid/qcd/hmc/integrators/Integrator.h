/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/hmc/integrators/Integrator.h

Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Guido Cossu <cossu@post.kek.jp>

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
#ifndef INTEGRATOR_INCLUDED
#define INTEGRATOR_INCLUDED

#include <memory>

namespace Grid {
namespace QCD {

class IntegratorParameters: Serializable {
public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(IntegratorParameters,
                std::string, name,      // name of the integrator
    unsigned int, MDsteps,  // number of outer steps
    RealD, trajL,           // trajectory length
  )

  IntegratorParameters(int MDsteps_ = 10, RealD trajL_ = 1.0)
      : MDsteps(MDsteps_),
        trajL(trajL_){
        // empty body constructor
  };


  template <class ReaderClass, typename std::enable_if<isReader<ReaderClass>::value, int >::type = 0 >
  IntegratorParameters(ReaderClass & Reader){
    std::cout << GridLogMessage << "Reading integrator\n";
    read(Reader, "Integrator", *this);
  }

  void print_parameters() const {
    std::cout << GridLogMessage << "[Integrator] Type               : " << name << std::endl;
    std::cout << GridLogMessage << "[Integrator] Trajectory length  : " << trajL << std::endl;
    std::cout << GridLogMessage << "[Integrator] Number of MD steps : " << MDsteps << std::endl;
    std::cout << GridLogMessage << "[Integrator] Step size          : " << trajL/MDsteps << std::endl;
  }
};

/*! @brief Class for Molecular Dynamics management */
template <class FieldImplementation, class SmearingPolicy, class RepresentationPolicy>
class Integrator {
 protected:
  typedef typename FieldImplementation::Field MomentaField;  //for readability
  typedef typename FieldImplementation::Field Field;

  int levels;  // number of integration levels
  double t_U;  // Track time passing on each level and for U and for P
  std::vector<double> t_P;  

  MomentaField P;
  SmearingPolicy& Smearer;
  RepresentationPolicy Representations;
  IntegratorParameters Params;

  const ActionSet<Field, RepresentationPolicy> as;

  void update_P(Field& U, int level, double ep) {
    t_P[level] += ep;
    update_P(P, U, level, ep);

    std::cout << GridLogIntegrator << "[" << level << "] P " << " dt " << ep << " : t_P " << t_P[level] << std::endl;
  }

  // to be used by the actionlevel class to iterate
  // over the representations
  struct _updateP {
    template <class FieldType, class GF, class Repr>
    void operator()(std::vector<Action<FieldType>*> repr_set, Repr& Rep,
                    GF& Mom, GF& U, double ep) {
      for (int a = 0; a < repr_set.size(); ++a) {
        FieldType forceR(U._grid);
        // Implement smearing only for the fundamental representation now
        repr_set.at(a)->deriv(Rep.U, forceR);
        GF force = Rep.RtoFundamentalProject(forceR);  // Ta for the fundamental rep
        Real force_abs = std::sqrt(norm2(force)/(U._grid->gSites()));
        std::cout << GridLogIntegrator << "Hirep Force average: " << force_abs << std::endl;
	Mom -= force * ep* HMC_MOMENTUM_DENOMINATOR;; 
      }
    }
  } update_P_hireps{};

  void update_P(MomentaField& Mom, Field& U, int level, double ep) {
    // input U actually not used in the fundamental case
    // Fundamental updates, include smearing

    for (int a = 0; a < as[level].actions.size(); ++a) {
      double start_full = usecond();
      Field force(U._grid);
      conformable(U._grid, Mom._grid);

      Field& Us = Smearer.get_U(as[level].actions.at(a)->is_smeared);
      double start_force = usecond();
      as[level].actions.at(a)->deriv(Us, force);  // deriv should NOT include Ta

      std::cout << GridLogIntegrator << "Smearing (on/off): " << as[level].actions.at(a)->is_smeared << std::endl;
      if (as[level].actions.at(a)->is_smeared) Smearer.smeared_force(force);
      force = FieldImplementation::projectForce(force); // Ta for gauge fields
      double end_force = usecond();
      Real force_abs = std::sqrt(norm2(force)/U._grid->gSites());
      std::cout << GridLogIntegrator << "["<<level<<"]["<<a<<"] Force average: " << force_abs << std::endl;
      Mom -= force * ep* HMC_MOMENTUM_DENOMINATOR;; 
      double end_full = usecond();
      double time_full  = (end_full - start_full) / 1e3;
      double time_force = (end_force - start_force) / 1e3;
      std::cout << GridLogMessage << "["<<level<<"]["<<a<<"] P update elapsed time: " << time_full << " ms (force: " << time_force << " ms)"  << std::endl;
    }

    // Force from the other representations
    as[level].apply(update_P_hireps, Representations, Mom, U, ep);
  }

  void update_U(Field& U, double ep) {
    update_U(P, U, ep);

    t_U += ep;
    int fl = levels - 1;
    std::cout << GridLogIntegrator << "   " << "[" << fl << "] U " << " dt " << ep << " : t_U " << t_U << std::endl;
  }
  
  void update_U(MomentaField& Mom, Field& U, double ep) {
    // exponential of Mom*U in the gauge fields case
    FieldImplementation::update_field(Mom, U, ep);

    // Update the smeared fields, can be implemented as observer
    Smearer.set_Field(U);

    // Update the higher representations fields
    Representations.update(U);  // void functions if fundamental representation
  }

  virtual void step(Field& U, int level, int first, int last) = 0;

 public:
  Integrator(GridBase* grid, IntegratorParameters Par,
             ActionSet<Field, RepresentationPolicy>& Aset,
             SmearingPolicy& Sm)
      : Params(Par),
        as(Aset),
        P(grid),
        levels(Aset.size()),
        Smearer(Sm),
        Representations(grid) {
    t_P.resize(levels, 0.0);
    t_U = 0.0;
    // initialization of smearer delegated outside of Integrator
  };

  virtual ~Integrator() {}

  virtual std::string integrator_name() = 0;

  void print_parameters(){
    std::cout << GridLogMessage << "[Integrator] Name : "<< integrator_name() << std::endl;
    Params.print_parameters();
  }

  void print_actions(){
        std::cout << GridLogMessage << ":::::::::::::::::::::::::::::::::::::::::" << std::endl;
        std::cout << GridLogMessage << "[Integrator] Action summary: "<<std::endl;
        for (int level = 0; level < as.size(); ++level) {
                std::cout << GridLogMessage << "[Integrator] ---- Level: "<< level << std::endl;
                for (int actionID = 0; actionID < as[level].actions.size(); ++actionID) {
                        std::cout << GridLogMessage << "["<< as[level].actions.at(actionID)->action_name() << "] ID: " << actionID << std::endl;
                        std::cout << as[level].actions.at(actionID)->LogParameters();
                }
        }
        std::cout << GridLogMessage << ":::::::::::::::::::::::::::::::::::::::::"<< std::endl;

  }

  void reverse_momenta(){
    P *= -1.0;
  }

  // to be used by the actionlevel class to iterate
  // over the representations
  struct _refresh {
    template <class FieldType, class Repr>
    void operator()(std::vector<Action<FieldType>*> repr_set, Repr& Rep,
                    GridParallelRNG& pRNG) {
      for (int a = 0; a < repr_set.size(); ++a){
        repr_set.at(a)->refresh(Rep.U, pRNG);
      
      std::cout << GridLogDebug << "Hirep refreshing pseudofermions" << std::endl;
    }
    }
  } refresh_hireps{};

  // Initialization of momenta and actions
  void refresh(Field& U, GridParallelRNG& pRNG) {
    assert(P._grid == U._grid);
    std::cout << GridLogIntegrator << "Integrator refresh\n";

    FieldImplementation::generate_momenta(P, pRNG);

    // Update the smeared fields, can be implemented as observer
    // necessary to keep the fields updated even after a reject
    // of the Metropolis
    Smearer.set_Field(U);
    // Set the (eventual) representations gauge fields
    Representations.update(U);

    // The Smearer is attached to a pointer of the gauge field
    // automatically gets the correct field
    // whether or not has been accepted in the previous sweep
    for (int level = 0; level < as.size(); ++level) {
      for (int actionID = 0; actionID < as[level].actions.size(); ++actionID) {
        // get gauge field from the SmearingPolicy and
        // based on the boolean is_smeared in actionID
        Field& Us = Smearer.get_U(as[level].actions.at(actionID)->is_smeared);
        as[level].actions.at(actionID)->refresh(Us, pRNG);
      }

      // Refresh the higher representation actions
      as[level].apply(refresh_hireps, Representations, pRNG);
    }
  }

  // to be used by the actionlevel class to iterate
  // over the representations
  struct _S {
    template <class FieldType, class Repr>
    void operator()(std::vector<Action<FieldType>*> repr_set, Repr& Rep, int level, RealD& H) {
      
      for (int a = 0; a < repr_set.size(); ++a) {
        RealD Hterm = repr_set.at(a)->S(Rep.U);
        std::cout << GridLogMessage << "S Level " << level << " term " << a << " H Hirep = " << Hterm << std::endl;
        H += Hterm;

      }
    }
  } S_hireps{};

  // Calculate action
  RealD S(Field& U) {  // here also U not used

    std::cout << GridLogIntegrator << "Integrator action\n";

    RealD H = - FieldImplementation::FieldSquareNorm(P)/HMC_MOMENTUM_DENOMINATOR; // - trace (P*P)/denom

    RealD Hterm;

    // Actions
    for (int level = 0; level < as.size(); ++level) {
      for (int actionID = 0; actionID < as[level].actions.size(); ++actionID) {
        // get gauge field from the SmearingPolicy and
        // based on the boolean is_smeared in actionID
        Field& Us = Smearer.get_U(as[level].actions.at(actionID)->is_smeared);
        std::cout << GridLogMessage << "S [" << level << "][" << actionID << "] action eval " << std::endl;
        Hterm = as[level].actions.at(actionID)->S(Us);
        std::cout << GridLogMessage << "S [" << level << "][" << actionID << "] H = " << Hterm << std::endl;
        H += Hterm;
      }
      as[level].apply(S_hireps, Representations, level, H);
    }

    return H;
  }

  void integrate(Field& U) {
    // reset the clocks
    t_U = 0;
    for (int level = 0; level < as.size(); ++level) {
      t_P[level] = 0;
    }

    for (int step = 0; step < Params.MDsteps; ++step) {  // MD step
      int first_step = (step == 0);
      int last_step = (step == Params.MDsteps - 1);
      this->step(U, 0, first_step, last_step);
    }

    // Check the clocks all match on all levels
    for (int level = 0; level < as.size(); ++level) {
      assert(fabs(t_U - t_P[level]) < 1.0e-6);  // must be the same
      std::cout << GridLogIntegrator << " times[" << level << "]= " << t_P[level] << " " << t_U << std::endl;
    }

    // and that we indeed got to the end of the trajectory
    assert(fabs(t_U - Params.trajL) < 1.0e-6);

  }


  

};
}
}
#endif  // INTEGRATOR_INCLUDED

