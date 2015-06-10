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
#include <qcd/action/fermion/FermionOperator.h>

////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////
#include <qcd/action/fermion/WilsonCompressor.h>     //used by all wilson type fermions
#include <qcd/action/fermion/WilsonKernels.h>        //used by all wilson type fermions

////////////////////////////////////////////
// 4D formulations
////////////////////////////////////////////
#include <qcd/action/fermion/WilsonFermion.h>
//#include <qcd/action/fermion/CloverFermion.h>

////////////////////////////////////////////
// 5D formulations...
////////////////////////////////////////////

#include <qcd/action/fermion/WilsonFermion5D.h> // used by all 5d overlap types

//////////
// Cayley
//////////
#include <qcd/action/fermion/CayleyFermion5D.h>

#include <qcd/action/fermion/DomainWallFermion.h>
#include <qcd/action/fermion/DomainWallFermion.h>

#include <qcd/action/fermion/MobiusFermion.h>
#include <qcd/action/fermion/ScaledShamirFermion.h>
#include <qcd/action/fermion/OverlapWilsonCayleyTanhFermion.h>

#include <qcd/action/fermion/MobiusZolotarevFermion.h>
#include <qcd/action/fermion/ShamirZolotarevFermion.h>
#include <qcd/action/fermion/OverlapWilsonCayleyZolotarevFermion.h>

//////////////////////
// Continued fraction
//////////////////////
#include <qcd/action/fermion/ContinuedFractionFermion5D.h>
#include <qcd/action/fermion/OverlapWilsonContfracTanhFermion.h>
#include <qcd/action/fermion/OverlapWilsonContfracZolotarevFermion.h>

//////////////////////
// Partial fraction
//////////////////////
#include <qcd/action/fermion/PartialFractionFermion5D.h>
#include <qcd/action/fermion/OverlapWilsonPartialFractionTanhFermion.h>
#include <qcd/action/fermion/OverlapWilsonPartialFractionZolotarevFermion.h>

///////////////////////////////////////////////////////////////////////////////
// G5 herm -- this has to live in QCD since dirac matrix is not in the broader sector of code
///////////////////////////////////////////////////////////////////////////////
#include <qcd/action/fermion/g5HermitianLinop.h>

#endif
