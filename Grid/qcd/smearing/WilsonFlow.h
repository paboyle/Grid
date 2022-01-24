/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/modules/plaquette.h

Copyright (C) 2017

Author: Guido Cossu <guido.cossu@ed.ac.uk>
Author: Christopher Kelly <ckelly@bnl.gov>

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

#pragma once

NAMESPACE_BEGIN(Grid);

template <class Gimpl>
class WilsonFlow: public Smear<Gimpl>{
public:
  //Store generic measurements to take during smearing process using std::function
  typedef std::function<void(int, RealD, const typename Gimpl::GaugeField &)> FunctionType;  //int: step,  RealD: flow time,  GaugeField : the gauge field
  
private:
  unsigned int Nstep;
  RealD epsilon; //for regular smearing this is the time step, for adaptive it is the initial time step
 
  std::vector< std::pair<int, FunctionType> > functions; //The int maps to the measurement frequency

  mutable WilsonGaugeAction<Gimpl> SG;

  //Evolve the gauge field by 1 step and update tau
  void evolve_step(typename Gimpl::GaugeField &U, RealD &tau) const;
  //Evolve the gauge field by 1 step and update tau and the current time step eps
  void evolve_step_adaptive(typename Gimpl::GaugeField&U, RealD &tau, RealD &eps, RealD maxTau) const;

public:
  INHERIT_GIMPL_TYPES(Gimpl)

  void resetActions(){ functions.clear(); }

  void addMeasurement(int meas_interval, FunctionType meas){ functions.push_back({meas_interval, meas}); }

  //Set the class to perform the default measurements: 
  //the plaquette energy density every step
  //the plaquette topological charge every 'topq_meas_interval' steps
  //and output to stdout
  void setDefaultMeasurements(int topq_meas_interval = 1);

  explicit WilsonFlow(unsigned int Nstep, RealD epsilon, unsigned int interval = 1):
  Nstep(Nstep),
    epsilon(epsilon),
    SG(WilsonGaugeAction<Gimpl>(3.0)) {
    // WilsonGaugeAction with beta 3.0
    assert(epsilon > 0.0);
    LogMessage();
    setDefaultMeasurements(interval);
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

  void smear_adaptive(GaugeField&, const GaugeField&, RealD maxTau) const;

  //Compute t^2 <E(t)> for time t from the plaquette
  static RealD energyDensityPlaquette(const RealD t, const GaugeField& U);

  //Compute t^2 <E(t)> for time t from the 1x1 cloverleaf form
  //t is the Wilson flow time
  static RealD energyDensityCloverleaf(const RealD t, const GaugeField& U);
  
  //Evolve the gauge field by Nstep steps of epsilon and return the energy density computed every interval steps
  //The smeared field is output as V
  std::vector<RealD> flowMeasureEnergyDensityPlaquette(GaugeField &V, const GaugeField& U, int measure_interval = 1);

  //Version that does not return the smeared field
  std::vector<RealD> flowMeasureEnergyDensityPlaquette(const GaugeField& U, int measure_interval = 1);


  //Evolve the gauge field by Nstep steps of epsilon and return the Cloverleaf energy density computed every interval steps
  //The smeared field is output as V
  std::vector<RealD> flowMeasureEnergyDensityCloverleaf(GaugeField &V, const GaugeField& U, int measure_interval = 1);

  //Version that does not return the smeared field
  std::vector<RealD> flowMeasureEnergyDensityCloverleaf(const GaugeField& U, int measure_interval = 1);
};


////////////////////////////////////////////////////////////////////////////////
// Implementations
////////////////////////////////////////////////////////////////////////////////
template <class Gimpl>
void WilsonFlow<Gimpl>::evolve_step(typename Gimpl::GaugeField &U, RealD &tau) const{
  GaugeField Z(U.Grid());
  GaugeField tmp(U.Grid());
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
  tau += epsilon;
}

template <class Gimpl>
void WilsonFlow<Gimpl>::evolve_step_adaptive(typename Gimpl::GaugeField &U, RealD &tau, RealD &eps, RealD maxTau) const{
  if (maxTau - tau < eps){
    eps = maxTau-tau;
  }
  //std::cout << GridLogMessage << "Integration epsilon : " << epsilon << std::endl;
  GaugeField Z(U.Grid());
  GaugeField Zprime(U.Grid());
  GaugeField tmp(U.Grid()), Uprime(U.Grid());
  Uprime = U;
  SG.deriv(U, Z);
  Zprime = -Z;
  Z *= 0.25;                                  // Z0 = 1/4 * F(U)
  Gimpl::update_field(Z, U, -2.0*eps);    // U = W1 = exp(ep*Z0)*W0

  Z *= -17.0/8.0;
  SG.deriv(U, tmp); Z += tmp;                 // -17/32*Z0 +Z1
  Zprime += 2.0*tmp;
  Z *= 8.0/9.0;                               // Z = -17/36*Z0 +8/9*Z1
  Gimpl::update_field(Z, U, -2.0*eps);    // U_= W2 = exp(ep*Z)*W1
    

  Z *= -4.0/3.0;
  SG.deriv(U, tmp); Z += tmp;                 // 4/3*(17/36*Z0 -8/9*Z1) +Z2
  Z *= 3.0/4.0;                               // Z = 17/36*Z0 -8/9*Z1 +3/4*Z2
  Gimpl::update_field(Z, U, -2.0*eps);    // V(t+e) = exp(ep*Z)*W2

  // Ramos 
  Gimpl::update_field(Zprime, Uprime, -2.0*eps); // V'(t+e) = exp(ep*Z')*W0
  // Compute distance as norm^2 of the difference
  GaugeField diffU = U - Uprime;
  RealD diff = norm2(diffU);
  // adjust integration step
    
  tau += eps;
  //std::cout << GridLogMessage << "Adjusting integration step with distance: " << diff << std::endl;
    
  eps = eps*0.95*std::pow(1e-4/diff,1./3.);
  //std::cout << GridLogMessage << "New epsilon : " << epsilon << std::endl;

}


template <class Gimpl>
RealD WilsonFlow<Gimpl>::energyDensityPlaquette(const RealD t, const GaugeField& U){
  static WilsonGaugeAction<Gimpl> SG(3.0);
  return 2.0 * t * t * SG.S(U)/U.Grid()->gSites();
}

//Compute t^2 <E(t)> for time from the 1x1 cloverleaf form
template <class Gimpl>
RealD WilsonFlow<Gimpl>::energyDensityCloverleaf(const RealD t, const GaugeField& U){
  typedef typename Gimpl::GaugeLinkField GaugeMat;
  typedef typename Gimpl::GaugeField GaugeLorentz;

  assert(Nd == 4);
  //E = 1/2 tr( F_munu F_munu )
  //However as  F_numu = -F_munu, only need to sum the trace of the squares of the following 6 field strengths:
  //F_01 F_02 F_03   F_12 F_13  F_23
  GaugeMat F(U.Grid());
  LatticeComplexD R(U.Grid());
  R = Zero();
  
  for(int mu=0;mu<3;mu++){
    for(int nu=mu+1;nu<4;nu++){
      WilsonLoops<Gimpl>::FieldStrength(F, U, mu, nu);
      R = R + trace(F*F);
    }
  }
  ComplexD out = sum(R);
  out = t*t*out / RealD(U.Grid()->gSites());
  return -real(out); //minus sign necessary for +ve energy
}


template <class Gimpl>
std::vector<RealD> WilsonFlow<Gimpl>::flowMeasureEnergyDensityPlaquette(GaugeField &V, const GaugeField& U, int measure_interval){
  std::vector<RealD> out;
  resetActions();
  addMeasurement(measure_interval, [&out](int step, RealD t, const typename Gimpl::GaugeField &U){ 
      std::cout << GridLogMessage << "[WilsonFlow] Computing plaquette energy density for step " << step << std::endl;
      out.push_back( energyDensityPlaquette(t,U) );
    });      
  smear(V,U);
  return out;
}

template <class Gimpl>
std::vector<RealD> WilsonFlow<Gimpl>::flowMeasureEnergyDensityPlaquette(const GaugeField& U, int measure_interval){
  GaugeField V(U);
  return flowMeasureEnergyDensityPlaquette(V,U, measure_interval);
}

template <class Gimpl>
std::vector<RealD> WilsonFlow<Gimpl>::flowMeasureEnergyDensityCloverleaf(GaugeField &V, const GaugeField& U, int measure_interval){
  std::vector<RealD> out;
  resetActions();
  addMeasurement(measure_interval, [&out](int step, RealD t, const typename Gimpl::GaugeField &U){ 
      std::cout << GridLogMessage << "[WilsonFlow] Computing Cloverleaf energy density for step " << step << std::endl;
      out.push_back( energyDensityCloverleaf(t,U) );
    });      
  smear(V,U);
  return out;
}

template <class Gimpl>
std::vector<RealD> WilsonFlow<Gimpl>::flowMeasureEnergyDensityCloverleaf(const GaugeField& U, int measure_interval){
  GaugeField V(U);
  return flowMeasureEnergyDensityCloverleaf(V,U, measure_interval);
}



//#define WF_TIMING 
template <class Gimpl>
void WilsonFlow<Gimpl>::smear(GaugeField& out, const GaugeField& in) const{
  out = in;
  RealD taus = 0.;
  for (unsigned int step = 1; step <= Nstep; step++) { //step indicates the number of smearing steps applied at the time of measurement
    auto start = std::chrono::high_resolution_clock::now();
    evolve_step(out, taus);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
#ifdef WF_TIMING
    std::cout << "Time to evolve " << diff.count() << " s\n";
#endif
    //Perform measurements
    for(auto const &meas : functions)
      if( step % meas.first == 0 ) meas.second(step,taus,out);
  }
}

template <class Gimpl>
void WilsonFlow<Gimpl>::smear_adaptive(GaugeField& out, const GaugeField& in, RealD maxTau) const{
  out = in;
  RealD taus = 0.;
  RealD eps = epsilon;
  unsigned int step = 0;
  do{
    step++;
    //std::cout << GridLogMessage << "Evolution time :"<< taus << std::endl;
    evolve_step_adaptive(out, taus, eps, maxTau);
    //Perform measurements
    for(auto const &meas : functions)
      if( step % meas.first == 0 ) meas.second(step,taus,out);
  } while (taus < maxTau);
}

template <class Gimpl>
void WilsonFlow<Gimpl>::setDefaultMeasurements(int topq_meas_interval){
  addMeasurement(1, [](int step, RealD t, const typename Gimpl::GaugeField &U){
      std::cout << GridLogMessage << "[WilsonFlow] Energy density (plaq) : "  << step << "  " << t << "  " << energyDensityPlaquette(t,U) << std::endl;
    });
  addMeasurement(topq_meas_interval, [](int step, RealD t, const typename Gimpl::GaugeField &U){
      std::cout << GridLogMessage << "[WilsonFlow] Top. charge           : "  << step << "  " << WilsonLoops<Gimpl>::TopologicalCharge(U) << std::endl;
    });
}


NAMESPACE_END(Grid);

