/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/hmc/integrators/Integrator_algorithm.h

Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>

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


/*! @file Integrator_algorithm.h
 * @brief Declaration of classes for the Molecular Dynamics algorithms
 *
 */
//--------------------------------------------------------------------

#ifndef INTEGRATOR_ALG_INCLUDED
#define INTEGRATOR_ALG_INCLUDED

NAMESPACE_BEGIN(Grid);

/* PAB:
 *
 * Recursive leapfrog; explanation of nested stepping
 *
 * Nested 1:4; units in dt for top level integrator
 *
 * CHROMA                           IroIro
 *   0        1                      0
 *  P 1/2                           P 1/2
 *          P 1/16                                  P1/16
 *                 U 1/8                                   U1/8
 *          P 1/8                                   P1/8
 *                 U 1/8                                   U1/8
 *          P 1/8                                   P1/8
 *                 U 1/8                                 U1/8
 *          P 1/8                                   P1/8
 *                 U 1/8                                 U1/8
 *          P 1/16                                  P1/8
 *  P 1                             P 1
 *          P 1/16                    * skipped --- avoids revaluating force
 *                 U 1/8                                 U1/8
 *          P 1/8                                   P1/8
 *                 U 1/8                                 U1/8
 *          P 1/8                                   P1/8
 *                 U 1/8                                 U1/8
 *          P 1/8                                   P1/8
 *                 U 1/8                                 U1/8
 *          P 1/16                                  P1/8
 *  P 1                             P 1
 *          P 1/16                    * skipped
 *                 U 1/8                                 U1/8
 *          P 1/8                                   P1/8
 *                 U 1/8                                 U1/8
 *          P 1/8                                   P1/8
 *                 U 1/8                                 U1/8
 *          P 1/8                                   P1/8
 *                 U 1/8                                 U1/8
 *          P 1/16                    * skipped
 *  P 1                             P 1
 *          P 1/16                                  P1/8
 *                 U 1/8                                 U1/8
 *          P 1/8                                   P1/8
 *                 U 1/8                                 U1/8
 *          P 1/8                                   P1/8
 *                 U 1/8                                 U1/8
 *          P 1/8                                   P1/8
 *                 U 1/8                                 U1/8
 *          P 1/16                                  P1/16
 *  P 1/2                            P 1/2
 */

template <class FieldImplementation_, class SmearingPolicy, class RepresentationPolicy = Representations<FundamentalRepresentation> >
class LeapFrog : public Integrator<FieldImplementation_, SmearingPolicy, RepresentationPolicy> 
{
public:
  typedef FieldImplementation_ FieldImplementation;
  typedef LeapFrog<FieldImplementation, SmearingPolicy, RepresentationPolicy> Algorithm;
  INHERIT_FIELD_TYPES(FieldImplementation);

  std::string integrator_name(){return "LeapFrog";}

  LeapFrog(GridBase* grid, IntegratorParameters Par, ActionSet<Field, RepresentationPolicy>& Aset, SmearingPolicy& Sm, Metric<Field>& M)
    : Integrator<FieldImplementation, SmearingPolicy, RepresentationPolicy>(grid, Par, Aset, Sm,M){};

  void step(Field& U, int level, int _first, int _last) {
    int fl = this->as.size() - 1;
    // level  : current level
    // fl     : final level
    // eps    : current step size

    // Get current level step size
    RealD eps = this->Params.trajL/this->Params.MDsteps;
    for (int l = 0; l <= level; ++l) eps /= this->as[l].multiplier;

    int multiplier = this->as[level].multiplier;
    for (int e = 0; e < multiplier; ++e) {
      int first_step = _first && (e == 0);
      int last_step = _last && (e == multiplier - 1);

      if (first_step) {  // initial half step
        this->update_P(U, level, eps / 2.0);
      }

      if (level == fl) {  // lowest level
        this->update_U(U, eps);
      } else {  // recursive function call
        this->step(U, level + 1, first_step, last_step);
      }

      int mm = last_step ? 1 : 2;
      this->update_P(U, level, mm * eps / 2.0);
    }
  }
};

template <class FieldImplementation_, class SmearingPolicy, class RepresentationPolicy = Representations<FundamentalRepresentation> >
class MinimumNorm2 : public Integrator<FieldImplementation_, SmearingPolicy, RepresentationPolicy> 
{
private:
//  const RealD lambda = 0.1931833275037836;

public:
  typedef FieldImplementation_ FieldImplementation;
  INHERIT_FIELD_TYPES(FieldImplementation);

  MinimumNorm2(GridBase* grid, IntegratorParameters Par, ActionSet<Field, RepresentationPolicy>& Aset, SmearingPolicy& Sm, Metric<Field>& M)
    : Integrator<FieldImplementation, SmearingPolicy, RepresentationPolicy>(grid, Par, Aset, Sm,M){};

  std::string integrator_name(){return "MininumNorm2";}

  void step(Field& U, int level, int _first, int _last) {
    // level  : current level
    // fl     : final level
    // eps    : current step size
    assert(level<3);
    RealD lambda= this->Params.lambda0;
    if (level>0) lambda= this->Params.lambda1;
    if (level>1) lambda= this->Params.lambda2;
    std::cout << GridLogMessage << "level: "<<level<< "lambda: "<<lambda<<std::endl;

    int fl = this->as.size() - 1;

    RealD eps = this->Params.trajL/this->Params.MDsteps * 2.0;
    for (int l = 0; l <= level; ++l) eps /= 2.0 * this->as[l].multiplier;

    // Nesting:  2xupdate_U of size eps/2
    // Next level is eps/2/multiplier

    int multiplier = this->as[level].multiplier;
    for (int e = 0; e < multiplier; ++e) {  // steps per step

      int first_step = _first && (e == 0);
      int last_step = _last && (e == multiplier - 1);

      if (first_step) {  // initial half step
        this->update_P(U, level, lambda * eps);
      }

      if (level == fl) {  // lowest level
        this->update_U(U, 0.5 * eps);
      } else {  // recursive function call
        this->step(U, level + 1, first_step, 0);
      }

      this->update_P(U, level, (1.0 - 2.0 * lambda) * eps);

      if (level == fl) {  // lowest level
        this->update_U(U, 0.5 * eps);
      } else {  // recursive function call
        this->step(U, level + 1, 0, last_step);
      }

      int mm = (last_step) ? 1 : 2;
      this->update_P(U, level, lambda * eps * mm);
    }
  }
};

template <class FieldImplementation_, class SmearingPolicy, class RepresentationPolicy = Representations<FundamentalRepresentation> >
class ForceGradient : public Integrator<FieldImplementation_, SmearingPolicy, RepresentationPolicy> 
{
private:
  const RealD lambda = 1.0 / 6.0;
  const RealD chi = 1.0 / 72.0;
  const RealD xi = 0.0;
  const RealD theta = 0.0;

public:
  typedef FieldImplementation_ FieldImplementation;
  INHERIT_FIELD_TYPES(FieldImplementation);

  // Looks like dH scales as dt^4. tested wilson/wilson 2 level.
  ForceGradient(GridBase* grid, IntegratorParameters Par,
                ActionSet<Field, RepresentationPolicy>& Aset,
                SmearingPolicy& Sm, Metric<Field>& M)
    : Integrator<FieldImplementation, SmearingPolicy, RepresentationPolicy>(
									    grid, Par, Aset, Sm,M){};

  std::string integrator_name(){return "ForceGradient";}
  
  void FG_update_P(Field& U, int level, double fg_dt, double ep) {
    Field Ufg(U.Grid());
    Field Pfg(U.Grid());
    Ufg = U;
    Pfg = Zero();
    std::cout << GridLogIntegrator << "FG update " << fg_dt << " " << ep << std::endl;
    // prepare_fg; no prediction/result cache for now
    // could relax CG stopping conditions for the
    // derivatives in the small step since the force gets multiplied by
    // a tiny dt^2 term relative to main force.
    //
    // Presently 4 force evals, and should have 3, so 1.33x too expensive.
    // could reduce this with sloppy CG to perhaps 1.15x too expensive
    // even without prediction.
    this->update_P(Pfg, Ufg, level, fg_dt);
    Pfg = Pfg*(1.0/fg_dt);
    this->update_U(Pfg, Ufg, fg_dt);
    this->update_P(Ufg, level, ep);
  }

  void step(Field& U, int level, int _first, int _last) {
    RealD eps = this->Params.trajL/this->Params.MDsteps * 2.0;
    for (int l = 0; l <= level; ++l) eps /= 2.0 * this->as[l].multiplier;

    RealD Chi = chi * eps * eps * eps;

    int fl = this->as.size() - 1;

    int multiplier = this->as[level].multiplier;

    for (int e = 0; e < multiplier; ++e) {  // steps per step

      int first_step = _first && (e == 0);
      int last_step = _last && (e == multiplier - 1);

      if (first_step) {  // initial half step
        this->update_P(U, level, lambda * eps);
      }

      if (level == fl) {  // lowest level
        this->update_U(U, 0.5 * eps);
      } else {  // recursive function call
        this->step(U, level + 1, first_step, 0);
      }

      this->FG_update_P(U, level, 2 * Chi / ((1.0 - 2.0 * lambda) * eps), (1.0 - 2.0 * lambda) * eps);

      if (level == fl) {  // lowest level
        this->update_U(U, 0.5 * eps);
      } else {  // recursive function call
        this->step(U, level + 1, 0, last_step);
      }

      int mm = (last_step) ? 1 : 2;
      this->update_P(U, level, lambda * eps * mm);
    }
  }
};

////////////////////////////////
// Riemannian Manifold HMC
// Girolami et al
////////////////////////////////



// correct
template <class FieldImplementation, class SmearingPolicy,
          class RepresentationPolicy =
              Representations<FundamentalRepresentation> >
class ImplicitLeapFrog : public Integrator<FieldImplementation, SmearingPolicy,
                                           RepresentationPolicy> {
 public:
  typedef ImplicitLeapFrog<FieldImplementation, SmearingPolicy, RepresentationPolicy>
      Algorithm;
  INHERIT_FIELD_TYPES(FieldImplementation);

  // Riemannian manifold metric operator
  // Hermitian operator Fisher

  std::string integrator_name(){return "ImplicitLeapFrog";}

  ImplicitLeapFrog(GridBase* grid, IntegratorParameters Par,
           ActionSet<Field, RepresentationPolicy>& Aset, SmearingPolicy& Sm, Metric<Field>& M)
      : Integrator<FieldImplementation, SmearingPolicy, RepresentationPolicy>(
            grid, Par, Aset, Sm, M){};

  void step(Field& U, int level, int _first, int _last) {
    int fl = this->as.size() - 1;
    // level  : current level
    // fl     : final level
    // eps    : current step size

    // Get current level step size
    RealD eps = this->Params.trajL/this->Params.MDsteps;
    for (int l = 0; l <= level; ++l) eps /= this->as[l].multiplier;

    int multiplier = this->as[level].multiplier;
    for (int e = 0; e < multiplier; ++e) {
      int first_step = _first && (e == 0);
      int last_step = _last && (e == multiplier - 1);

      if (first_step) {  // initial half step
       this->implicit_update_P(U, level, eps / 2.0);
      }

      if (level == fl) {  // lowest level
        this->implicit_update_U(U, eps,eps/2.);
      } else {  // recursive function call
        this->step(U, level + 1, first_step, last_step);
      }

      //int mm = last_step ? 1 : 2;
      if (last_step){
        this->update_P2(U, level, eps / 2.0);
      } else {
      this->implicit_update_P(U, level, eps, true);// works intermediate step
      }
    }
  }
};


template <class FieldImplementation, class SmearingPolicy,
          class RepresentationPolicy =
              Representations<FundamentalRepresentation> >
class ImplicitMinimumNorm2 : public Integrator<FieldImplementation, SmearingPolicy,
                                       RepresentationPolicy> {
 private:
//  const RealD lambda = 0.1931833275037836;

 public:
  INHERIT_FIELD_TYPES(FieldImplementation);

  ImplicitMinimumNorm2(GridBase* grid, IntegratorParameters Par,
               ActionSet<Field, RepresentationPolicy>& Aset, SmearingPolicy& Sm, Metric<Field>& M)
      : Integrator<FieldImplementation, SmearingPolicy, RepresentationPolicy>(
            grid, Par, Aset, Sm, M){};

  std::string integrator_name(){return "ImplicitMininumNorm2";}

  void step(Field& U, int level, int _first, int _last) {
    // level  : current level
    // fl     : final level
    // eps    : current step size

    int fl = this->as.size() - 1;
//    assert(Params.lambda.size()>level);
//    RealD lambda= Params.lambda[level];
    assert(level<3);
    RealD lambda= this->Params.lambda0;
    if (level>0) lambda= this->Params.lambda1;
    if (level>1) lambda= this->Params.lambda2;
    std::cout << GridLogMessage << "level: "<<level<< "lambda: "<<lambda<<std::endl;

  if(level<fl){

    RealD eps = this->Params.trajL/this->Params.MDsteps * 2.0;
    for (int l = 0; l <= level; ++l) eps /= 2.0 * this->as[l].multiplier;

    // Nesting:  2xupdate_U of size eps/2
    // Next level is eps/2/multiplier

    int multiplier = this->as[level].multiplier;
    for (int e = 0; e < multiplier; ++e) {  // steps per step

      int first_step = _first && (e == 0);
      int last_step = _last && (e == multiplier - 1);

      if (first_step) {  // initial half step
        this->update_P(U, level, lambda * eps);
      }

        this->step(U, level + 1, first_step, 0);

      this->update_P(U, level, (1.0 - 2.0 * lambda) * eps);

        this->step(U, level + 1, 0, last_step);

      int mm = (last_step) ? 1 : 2;
      this->update_P(U, level, lambda * eps * mm);
    }
  } 
  else 
  { // last level
    RealD eps = this->Params.trajL/this->Params.MDsteps * 2.0;
    for (int l = 0; l <= level; ++l) eps /= 2.0 * this->as[l].multiplier;

    // Nesting:  2xupdate_U of size eps/2
    // Next level is eps/2/multiplier

    int multiplier = this->as[level].multiplier;
    for (int e = 0; e < multiplier; ++e) {  // steps per step

      int first_step = _first && (e == 0);
      int last_step = _last && (e == multiplier - 1);

      if (first_step) {  // initial half step
        this->implicit_update_P(U, level, lambda * eps);
      }

      this->implicit_update_U(U, 0.5 * eps,lambda*eps);

      this->implicit_update_P(U, level, (1.0 - 2.0 * lambda) * eps, true);

      this->implicit_update_U(U, 0.5 * eps, (0.5-lambda)*eps);

      if (last_step) {
        this->update_P2(U, level, eps * lambda);
      } else {
        this->implicit_update_P(U, level, lambda * eps*2.0, true);
      }
    }
  }

  }
};

template <class FieldImplementation, class SmearingPolicy,
          class RepresentationPolicy =
              Representations<FundamentalRepresentation> >
class ImplicitCampostrini : public Integrator<FieldImplementation, SmearingPolicy,
                                       RepresentationPolicy> {
 private:
//  const RealD lambda = 0.1931833275037836;

 public:
  INHERIT_FIELD_TYPES(FieldImplementation);

  ImplicitCampostrini(GridBase* grid, IntegratorParameters Par,
               ActionSet<Field, RepresentationPolicy>& Aset, SmearingPolicy& Sm, Metric<Field>& M)
      : Integrator<FieldImplementation, SmearingPolicy, RepresentationPolicy>(
            grid, Par, Aset, Sm, M){};

  std::string integrator_name(){return "ImplicitCampostrini";}

  void step(Field& U, int level, int _first, int _last) {
    // level  : current level
    // fl     : final level
    // eps    : current step size

    int fl = this->as.size() - 1;
//    assert(Params.lambda.size()>level);
//    RealD lambda= Params.lambda[level];
    assert(level<3);
    RealD lambda= this->Params.lambda0;
    if (level>0) lambda= this->Params.lambda1;
    if (level>1) lambda= this->Params.lambda2;
    std::cout << GridLogMessage << "level: "<<level<< "lambda: "<<lambda<<std::endl;
    
    RealD sigma=pow(2.0,1./3.);

  if(level<fl){
//Still Omelyan. Needs to change step() to accept variable stepsize
    RealD eps = this->Params.trajL/this->Params.MDsteps * 2.0;
    for (int l = 0; l <= level; ++l) eps /= 2.0 * this->as[l].multiplier;

    // Nesting:  2xupdate_U of size eps/2
    // Next level is eps/2/multiplier

    int multiplier = this->as[level].multiplier;
    for (int e = 0; e < multiplier; ++e) {  // steps per step

      int first_step = _first && (e == 0);
      int last_step = _last && (e == multiplier - 1);

      if (first_step) {  // initial half step
        this->update_P(U, level, lambda * eps);
      }

        this->step(U, level + 1, first_step, 0);

      this->update_P(U, level, (1.0 - 2.0 * lambda) * eps);

        this->step(U, level + 1, 0, last_step);

      int mm = (last_step) ? 1 : 2;
      this->update_P(U, level, lambda * eps * mm);
    }
  } 
  else 
  { // last level
    RealD dt = this->Params.trajL/this->Params.MDsteps * 2.0;
    for (int l = 0; l <= level; ++l) dt /= 2.0 * this->as[l].multiplier;

    RealD epsilon = dt/(2.0 - sigma);

    int multiplier = this->as[level].multiplier;
    for (int e = 0; e < multiplier; ++e) {  // steps per step

      int first_step = _first && (e == 0);
      int last_step = _last && (e == multiplier - 1);
      // initial half step
      if (first_step) {  this->implicit_update_P(U, level, epsilon*0.5); }
      this->implicit_update_U(U, epsilon,epsilon*0.5);
      this->implicit_update_P(U, level, (1.0 - sigma) * epsilon *0.5, epsilon*0.5, true);
      this->implicit_update_U(U, -epsilon*sigma, -epsilon*sigma*0.5);
      this->implicit_update_P(U, level, (1.0 - sigma) * epsilon *0.5, -epsilon*sigma*0.5, true);
      this->implicit_update_U(U, epsilon,epsilon*0.5);
      if (last_step) { this->update_P2(U, level, epsilon*0.5 ); } 
      else
      this->implicit_update_P(U, level, epsilon,epsilon*0.5);
    }
  }

  }
};

NAMESPACE_END(Grid);

#endif  // INTEGRATOR_INCLUDED
