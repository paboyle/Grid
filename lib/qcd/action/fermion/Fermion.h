    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/Fermion_base_aggregate.h

    Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>

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
#ifndef  GRID_QCD_FERMION_H
#define  GRID_QCD_FERMION_H

// * Linear operators             (Hermitian and non-hermitian)  .. my LinearOperator
// * System solvers               (Hermitian and non-hermitian)  .. my OperatorFunction
// * MultiShift System solvers    (Hermitian and non-hermitian)  .. my OperatorFunction

////////////////////////////////////////////
// Abstract base interface
////////////////////////////////////////////
#include <Grid/qcd/action/ActionBase.h>
#include <Grid/qcd/action/ActionParams.h>

////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////
#include <Grid/qcd/action/gauge/GaugeImpl.h>
#include <Grid/qcd/utils/WilsonLoops.h>

#include <Grid/qcd/action/fermion/WilsonCompressor.h>     //used by all wilson type fermions
#include <Grid/qcd/action/fermion/FermionOperatorImpl.h>
#include <Grid/qcd/action/fermion/FermionOperator.h>
#include <Grid/qcd/action/fermion/WilsonKernels.h>        //used by all wilson type fermions

////////////////////////////////////////////
// Gauge Actions
////////////////////////////////////////////
#include <Grid/qcd/action/gauge/Photon.h>
#include <Grid/qcd/action/gauge/WilsonGaugeAction.h>
#include <Grid/qcd/action/gauge/PlaqPlusRectangleAction.h>

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

////////////////////////////////////////////
// Fermion operators / actions
////////////////////////////////////////////

#include <Grid/qcd/action/fermion/WilsonFermion.h>       // 4d wilson like
#include <Grid/qcd/action/fermion/WilsonTMFermion.h>       // 4d wilson like
#include <Grid/qcd/action/fermion/WilsonFermion5D.h>     // 5d base used by all 5d overlap types
//#include <Grid/qcd/action/fermion/CloverFermion.h>
#include <Grid/qcd/action/fermion/ImprovedStaggeredFermion.h>
#include <Grid/qcd/action/fermion/ImprovedStaggeredFermion5D.h>
#include <Grid/qcd/action/fermion/CayleyFermion5D.h>     // Cayley types
#include <Grid/qcd/action/fermion/DomainWallFermion.h>
#include <Grid/qcd/action/fermion/DomainWallFermion.h>
#include <Grid/qcd/action/fermion/MobiusFermion.h>
#include <Grid/qcd/action/fermion/ZMobiusFermion.h>
#include <Grid/qcd/action/fermion/ScaledShamirFermion.h>
#include <Grid/qcd/action/fermion/MobiusZolotarevFermion.h>
#include <Grid/qcd/action/fermion/ShamirZolotarevFermion.h>
#include <Grid/qcd/action/fermion/OverlapWilsonCayleyTanhFermion.h>
#include <Grid/qcd/action/fermion/OverlapWilsonCayleyZolotarevFermion.h>
#include <Grid/qcd/action/fermion/ContinuedFractionFermion5D.h>               // Continued fraction
#include <Grid/qcd/action/fermion/OverlapWilsonContfracTanhFermion.h>
#include <Grid/qcd/action/fermion/OverlapWilsonContfracZolotarevFermion.h>
#include <Grid/qcd/action/fermion/PartialFractionFermion5D.h>                 // Partial fraction
#include <Grid/qcd/action/fermion/OverlapWilsonPartialFractionTanhFermion.h>
#include <Grid/qcd/action/fermion/OverlapWilsonPartialFractionZolotarevFermion.h>
///////////////////////////////////////////////////////////////////////////////
// G5 herm -- this has to live in QCD since dirac matrix is not in the broader sector of code
///////////////////////////////////////////////////////////////////////////////
#include <Grid/qcd/action/fermion/g5HermitianLinop.h>

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

typedef WilsonFermion<WilsonAdjImplR> WilsonAdjFermionR;
typedef WilsonFermion<WilsonAdjImplF> WilsonAdjFermionF;
typedef WilsonFermion<WilsonAdjImplD> WilsonAdjFermionD;

typedef WilsonFermion<WilsonTwoIndexSymmetricImplR> WilsonTwoIndexSymmetricFermionR;
typedef WilsonFermion<WilsonTwoIndexSymmetricImplF> WilsonTwoIndexSymmetricFermionF;
typedef WilsonFermion<WilsonTwoIndexSymmetricImplD> WilsonTwoIndexSymmetricFermionD;

typedef WilsonTMFermion<WilsonImplR> WilsonTMFermionR;
typedef WilsonTMFermion<WilsonImplF> WilsonTMFermionF;
typedef WilsonTMFermion<WilsonImplD> WilsonTMFermionD;

typedef DomainWallFermion<WilsonImplR> DomainWallFermionR;
typedef DomainWallFermion<WilsonImplF> DomainWallFermionF;
typedef DomainWallFermion<WilsonImplD> DomainWallFermionD;

typedef MobiusFermion<WilsonImplR> MobiusFermionR;
typedef MobiusFermion<WilsonImplF> MobiusFermionF;
typedef MobiusFermion<WilsonImplD> MobiusFermionD;

typedef ZMobiusFermion<ZWilsonImplR> ZMobiusFermionR;
typedef ZMobiusFermion<ZWilsonImplF> ZMobiusFermionF;
typedef ZMobiusFermion<ZWilsonImplD> ZMobiusFermionD;

// Ls vectorised 
typedef DomainWallFermion<DomainWallVec5dImplR> DomainWallFermionVec5dR;
typedef DomainWallFermion<DomainWallVec5dImplF> DomainWallFermionVec5dF;
typedef DomainWallFermion<DomainWallVec5dImplD> DomainWallFermionVec5dD;

typedef MobiusFermion<DomainWallVec5dImplR> MobiusFermionVec5dR;
typedef MobiusFermion<DomainWallVec5dImplF> MobiusFermionVec5dF;
typedef MobiusFermion<DomainWallVec5dImplD> MobiusFermionVec5dD;

typedef ZMobiusFermion<ZDomainWallVec5dImplR> ZMobiusFermionVec5dR;
typedef ZMobiusFermion<ZDomainWallVec5dImplF> ZMobiusFermionVec5dF;
typedef ZMobiusFermion<ZDomainWallVec5dImplD> ZMobiusFermionVec5dD;


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

typedef WilsonTMFermion<GparityWilsonImplR> GparityWilsonTMFermionR;
typedef WilsonTMFermion<GparityWilsonImplF> GparityWilsonTMFermionF;
typedef WilsonTMFermion<GparityWilsonImplD> GparityWilsonTMFermionD;
typedef MobiusFermion<GparityWilsonImplR> GparityMobiusFermionR;
typedef MobiusFermion<GparityWilsonImplF> GparityMobiusFermionF;
typedef MobiusFermion<GparityWilsonImplD> GparityMobiusFermionD;

typedef ImprovedStaggeredFermion<StaggeredImplR> ImprovedStaggeredFermionR;
typedef ImprovedStaggeredFermion<StaggeredImplF> ImprovedStaggeredFermionF;
typedef ImprovedStaggeredFermion<StaggeredImplD> ImprovedStaggeredFermionD;

typedef ImprovedStaggeredFermion5D<StaggeredImplR> ImprovedStaggeredFermion5DR;
typedef ImprovedStaggeredFermion5D<StaggeredImplF> ImprovedStaggeredFermion5DF;
typedef ImprovedStaggeredFermion5D<StaggeredImplD> ImprovedStaggeredFermion5DD;

typedef ImprovedStaggeredFermion5D<StaggeredVec5dImplR> ImprovedStaggeredFermionVec5dR;
typedef ImprovedStaggeredFermion5D<StaggeredVec5dImplF> ImprovedStaggeredFermionVec5dF;
typedef ImprovedStaggeredFermion5D<StaggeredVec5dImplD> ImprovedStaggeredFermionVec5dD;


  }}

#include <Grid/qcd/action/pseudofermion/OneFlavourRational.h>
#include <Grid/qcd/action/pseudofermion/OneFlavourRationalRatio.h>
#include <Grid/qcd/action/pseudofermion/OneFlavourEvenOddRational.h>
#include <Grid/qcd/action/pseudofermion/OneFlavourEvenOddRationalRatio.h>

////////////////////
// Scalar actions
////////////////////
#include <Grid/qcd/action/scalar/Scalar.h>

#endif
