#ifndef GRID_QCD_ACTIONS_H
#define GRID_QCD_ACTIONS_H


// Some reorganisation likely required as both Chroma and IroIro
// are separating the concept of the operator from that of action.
//
// The FermAction contains methods to create 
//
// * Linear operators             (Hermitian and non-hermitian)  .. my LinearOperator
// * System solvers               (Hermitian and non-hermitian)  .. my OperatorFunction
// * MultiShift System solvers    (Hermitian and non-hermitian)  .. my OperatorFunction


////////////////////////////////////////////
// Abstract base interface
////////////////////////////////////////////
#include <qcd/action/ActionBase.h>

////////////////////////////////////////////
// Gauge Actions
////////////////////////////////////////////
#include <qcd/action/gauge/WilsonGaugeAction.h>


////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////
#include <qcd/action/fermion/WilsonCompressor.h>     //used by all wilson type fermions
#include <qcd/action/fermion/FermionOperatorImpl.h>
#include <qcd/action/fermion/FermionOperator.h>
#include <qcd/action/fermion/WilsonKernels.h>        //used by all wilson type fermions


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
  template class A<GparityWilsonImplF>;		\
  template class A<GparityWilsonImplD>;		\
  template class A<WilsonImplF>;		\
  template class A<WilsonImplD>;

////////////////////////////////////////////
// Fermion operators / actions
////////////////////////////////////////////

#include <qcd/action/fermion/WilsonFermion.h>       // 4d wilson like
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
#include <qcd/action/fermion/OverlapWilsonContFracTanhFermion.h>
#include <qcd/action/fermion/OverlapWilsonContFracZolotarevFermion.h>

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
#include <qcd/action/pseudofermion/TwoFlavour.h>
#include <qcd/action/pseudofermion/TwoFlavourEvenOdd.h>

#endif
