    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/Actions.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#ifndef GRID_QCD_ACTIONS_H
#define GRID_QCD_ACTIONS_H

// * Linear operators             (Hermitian and non-hermitian)  .. my LinearOperator
// * System solvers               (Hermitian and non-hermitian)  .. my OperatorFunction
// * MultiShift System solvers    (Hermitian and non-hermitian)  .. my OperatorFunction

////////////////////////////////////////////
// Abstract base interface
////////////////////////////////////////////
#include <qcd/action/ActionBase.h>
#include <qcd/action/ActionParams.h>

////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////
#include <qcd/action/gauge/GaugeImpl.h>
#include <qcd/utils/WilsonLoops.h>

#include <qcd/action/fermion/WilsonCompressor.h>     //used by all wilson type fermions
#include <qcd/action/fermion/FermionOperatorImpl.h>
#include <qcd/action/fermion/FermionOperator.h>
#include <qcd/action/fermion/WilsonKernels.h>        //used by all wilson type fermions

////////////////////////////////////////////
// Gauge Actions
////////////////////////////////////////////
#include <qcd/action/gauge/WilsonGaugeAction.h>
#include <qcd/action/gauge/PlaqPlusRectangleAction.h>

namespace Grid {
namespace QCD {

typedef WilsonGaugeAction<PeriodicGimplR>          WilsonGaugeActionR;
typedef WilsonGaugeAction<PeriodicGimplF>          WilsonGaugeActionF;
typedef WilsonGaugeAction<PeriodicGimplD>          WilsonGaugeActionD;
typedef PlaqPlusRectangleAction<PeriodicGimplR>    PlaqPlusRectangleActionR;
typedef PlaqPlusRectangleAction<PeriodicGimplF>    PlaqPlusRectangleActionF;
typedef PlaqPlusRectangleAction<PeriodicGimplD>    PlaqPlusRectangleActionD;
typedef IwasakiGaugeAction<PeriodicGimplR>         IwasakiGaugeActionR;
typedef IwasakiGaugeAction<PeriodicGimplF>         IwasakiGaugeActionF;
typedef IwasakiGaugeAction<PeriodicGimplD>         IwasakiGaugeActionD;
typedef SymanzikGaugeAction<PeriodicGimplR>        SymanzikGaugeActionR;
typedef SymanzikGaugeAction<PeriodicGimplF>        SymanzikGaugeActionF;
typedef SymanzikGaugeAction<PeriodicGimplD>        SymanzikGaugeActionD;


typedef WilsonGaugeAction<ConjugateGimplR>          ConjugateWilsonGaugeActionR;
typedef WilsonGaugeAction<ConjugateGimplF>          ConjugateWilsonGaugeActionF;
typedef WilsonGaugeAction<ConjugateGimplD>          ConjugateWilsonGaugeActionD;
typedef PlaqPlusRectangleAction<ConjugateGimplR>    ConjugatePlaqPlusRectangleActionR;
typedef PlaqPlusRectangleAction<ConjugateGimplF>    ConjugatePlaqPlusRectangleActionF;
typedef PlaqPlusRectangleAction<ConjugateGimplD>    ConjugatePlaqPlusRectangleActionD;
typedef IwasakiGaugeAction<ConjugateGimplR>         ConjugateIwasakiGaugeActionR;
typedef IwasakiGaugeAction<ConjugateGimplF>         ConjugateIwasakiGaugeActionF;
typedef IwasakiGaugeAction<ConjugateGimplD>         ConjugateIwasakiGaugeActionD;
typedef SymanzikGaugeAction<ConjugateGimplR>        ConjugateSymanzikGaugeActionR;
typedef SymanzikGaugeAction<ConjugateGimplF>        ConjugateSymanzikGaugeActionF;
typedef SymanzikGaugeAction<ConjugateGimplD>        ConjugateSymanzikGaugeActionD;

}}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Explicit explicit template instantiation is still required in the .cc files
//
// - CayleyFermion5D.cc
// - PartialFractionFermion5D.cc
// - WilsonFermion5D.cc
// - WilsonKernelsHand.cc
// - ContinuedFractionFermion5D.cc
// - WilsonFermion.cc
// - WilsonKernels.cc
//
// The explicit instantiation is only avoidable if we move this source to headers and end up with include/parse/recompile
// for EVERY .cc file. This define centralises the list and restores global push of impl cases
////////////////////////////////////////////////////////////////////////////////////////////////////

#define FermOpTemplateInstantiate(A) \
  template class A<WilsonImplF>;		\
  template class A<WilsonImplD>;  \
  template class A<GparityWilsonImplF>;		\
  template class A<GparityWilsonImplD>;		

////////////////////////////////////////////
// Fermion operators / actions
////////////////////////////////////////////

#include <qcd/action/fermion/WilsonFermion.h>       // 4d wilson like
#include <qcd/action/fermion/WilsonTMFermion.h>       // 4d wilson like
#include <qcd/action/fermion/WilsonFermion5D.h>     // 5d base used by all 5d overlap types

//#include <qcd/action/fermion/CloverFermion.h>

#include <qcd/action/fermion/CayleyFermion5D.h>     // Cayley types
#include <qcd/action/fermion/DomainWallFermion.h>
#include <qcd/action/fermion/DomainWallFermion.h>
#include <qcd/action/fermion/MobiusFermion.h>
#include <qcd/action/fermion/ScaledShamirFermion.h>
#include <qcd/action/fermion/MobiusZolotarevFermion.h>
#include <qcd/action/fermion/ShamirZolotarevFermion.h>
#include <qcd/action/fermion/OverlapWilsonCayleyTanhFermion.h>
#include <qcd/action/fermion/OverlapWilsonCayleyZolotarevFermion.h>

#include <qcd/action/fermion/ContinuedFractionFermion5D.h>               // Continued fraction
#include <qcd/action/fermion/OverlapWilsonContfracTanhFermion.h>
#include <qcd/action/fermion/OverlapWilsonContfracZolotarevFermion.h>

#include <qcd/action/fermion/PartialFractionFermion5D.h>                 // Partial fraction
#include <qcd/action/fermion/OverlapWilsonPartialFractionTanhFermion.h>
#include <qcd/action/fermion/OverlapWilsonPartialFractionZolotarevFermion.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// More maintainable to maintain the following typedef list centrally, as more "impl" targets
// are added, (e.g. extension for gparity, half precision project in comms etc..)
////////////////////////////////////////////////////////////////////////////////////////////////////


// Cayley 5d
namespace Grid {
  namespace QCD {

typedef WilsonFermion<WilsonImplR> WilsonFermionR;
typedef WilsonFermion<WilsonImplF> WilsonFermionF;
typedef WilsonFermion<WilsonImplD> WilsonFermionD;

typedef WilsonTMFermion<WilsonImplR> WilsonTMFermionR;
typedef WilsonTMFermion<WilsonImplF> WilsonTMFermionF;
typedef WilsonTMFermion<WilsonImplD> WilsonTMFermionD;

typedef DomainWallFermion<WilsonImplR> DomainWallFermionR;
typedef DomainWallFermion<WilsonImplF> DomainWallFermionF;
typedef DomainWallFermion<WilsonImplD> DomainWallFermionD;
typedef MobiusFermion<WilsonImplR> MobiusFermionR;
typedef MobiusFermion<WilsonImplF> MobiusFermionF;
typedef MobiusFermion<WilsonImplD> MobiusFermionD;
typedef ScaledShamirFermion<WilsonImplR> ScaledShamirFermionR;
typedef ScaledShamirFermion<WilsonImplF> ScaledShamirFermionF;
typedef ScaledShamirFermion<WilsonImplD> ScaledShamirFermionD;

typedef MobiusZolotarevFermion<WilsonImplR> MobiusZolotarevFermionR;
typedef MobiusZolotarevFermion<WilsonImplF> MobiusZolotarevFermionF;
typedef MobiusZolotarevFermion<WilsonImplD> MobiusZolotarevFermionD;
typedef ShamirZolotarevFermion<WilsonImplR> ShamirZolotarevFermionR;
typedef ShamirZolotarevFermion<WilsonImplF> ShamirZolotarevFermionF;
typedef ShamirZolotarevFermion<WilsonImplD> ShamirZolotarevFermionD;

typedef OverlapWilsonCayleyTanhFermion<WilsonImplR> OverlapWilsonCayleyTanhFermionR;
typedef OverlapWilsonCayleyTanhFermion<WilsonImplF> OverlapWilsonCayleyTanhFermionF;
typedef OverlapWilsonCayleyTanhFermion<WilsonImplD> OverlapWilsonCayleyTanhFermionD;
typedef OverlapWilsonCayleyZolotarevFermion<WilsonImplR> OverlapWilsonCayleyZolotarevFermionR;
typedef OverlapWilsonCayleyZolotarevFermion<WilsonImplF> OverlapWilsonCayleyZolotarevFermionF;
typedef OverlapWilsonCayleyZolotarevFermion<WilsonImplD> OverlapWilsonCayleyZolotarevFermionD;

// Continued fraction
typedef OverlapWilsonContFracTanhFermion<WilsonImplR> OverlapWilsonContFracTanhFermionR;
typedef OverlapWilsonContFracTanhFermion<WilsonImplF> OverlapWilsonContFracTanhFermionF;
typedef OverlapWilsonContFracTanhFermion<WilsonImplD> OverlapWilsonContFracTanhFermionD;
typedef OverlapWilsonContFracZolotarevFermion<WilsonImplR> OverlapWilsonContFracZolotarevFermionR;
typedef OverlapWilsonContFracZolotarevFermion<WilsonImplF> OverlapWilsonContFracZolotarevFermionF;
typedef OverlapWilsonContFracZolotarevFermion<WilsonImplD> OverlapWilsonContFracZolotarevFermionD;

// Partial fraction
typedef OverlapWilsonPartialFractionTanhFermion<WilsonImplR> OverlapWilsonPartialFractionTanhFermionR;
typedef OverlapWilsonPartialFractionTanhFermion<WilsonImplF> OverlapWilsonPartialFractionTanhFermionF;
typedef OverlapWilsonPartialFractionTanhFermion<WilsonImplD> OverlapWilsonPartialFractionTanhFermionD;

typedef OverlapWilsonPartialFractionZolotarevFermion<WilsonImplR> OverlapWilsonPartialFractionZolotarevFermionR;
typedef OverlapWilsonPartialFractionZolotarevFermion<WilsonImplF> OverlapWilsonPartialFractionZolotarevFermionF;
typedef OverlapWilsonPartialFractionZolotarevFermion<WilsonImplD> OverlapWilsonPartialFractionZolotarevFermionD;

// Gparity cases; partial list until tested
typedef WilsonFermion<GparityWilsonImplR>     GparityWilsonFermionR;
typedef WilsonFermion<GparityWilsonImplF>     GparityWilsonFermionF;
typedef WilsonFermion<GparityWilsonImplD>     GparityWilsonFermionD;
typedef DomainWallFermion<GparityWilsonImplR> GparityDomainWallFermionR;
typedef DomainWallFermion<GparityWilsonImplF> GparityDomainWallFermionF;
typedef DomainWallFermion<GparityWilsonImplD> GparityDomainWallFermionD;

  }}
///////////////////////////////////////////////////////////////////////////////
// G5 herm -- this has to live in QCD since dirac matrix is not in the broader sector of code
///////////////////////////////////////////////////////////////////////////////
#include <qcd/action/fermion/g5HermitianLinop.h>

////////////////////////////////////////
// Pseudo fermion combinations for HMC
////////////////////////////////////////
#include <qcd/action/pseudofermion/EvenOddSchurDifferentiable.h>

#include <qcd/action/pseudofermion/TwoFlavour.h>
#include <qcd/action/pseudofermion/TwoFlavourRatio.h>
#include <qcd/action/pseudofermion/TwoFlavourEvenOdd.h>
#include <qcd/action/pseudofermion/TwoFlavourEvenOddRatio.h>

#include <qcd/action/pseudofermion/OneFlavourRational.h>
#include <qcd/action/pseudofermion/OneFlavourRationalRatio.h>
#include <qcd/action/pseudofermion/OneFlavourEvenOddRational.h>
#include <qcd/action/pseudofermion/OneFlavourEvenOddRationalRatio.h>

#endif
