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
#pragma once

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
// - DomainWallEOFAFermion.cc
// - MobiusEOFAFermion.cc
//
// The explicit instantiation is only avoidable if we move this source to headers and end up with include/parse/recompile
// for EVERY .cc file. This define centralises the list and restores global push of impl cases
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////
// Fermion operators / actions
////////////////////////////////////////////

#include <Grid/qcd/action/fermion/WilsonFermion.h>       // 4d wilson like
NAMESPACE_CHECK(Wilson);
#include <Grid/qcd/action/fermion/WilsonTMFermion.h>       // 4d wilson like
NAMESPACE_CHECK(WilsonTM);
#include <Grid/qcd/action/fermion/WilsonCloverFermion.h> // 4d wilson clover fermions
NAMESPACE_CHECK(WilsonClover);
#include <Grid/qcd/action/fermion/WilsonFermion5D.h>     // 5d base used by all 5d overlap types
NAMESPACE_CHECK(Wilson5D);

#include <Grid/qcd/action/fermion/NaiveStaggeredFermion.h>
#include <Grid/qcd/action/fermion/ImprovedStaggeredFermion.h>
#include <Grid/qcd/action/fermion/ImprovedStaggeredFermion5D.h>
NAMESPACE_CHECK(Staggered);

#include <Grid/qcd/action/fermion/CayleyFermion5D.h>     // Cayley types
#include <Grid/qcd/action/fermion/DomainWallFermion.h>
#include <Grid/qcd/action/fermion/DomainWallEOFAFermion.h>
#include <Grid/qcd/action/fermion/MobiusFermion.h>
#include <Grid/qcd/action/fermion/MobiusEOFAFermion.h>
#include <Grid/qcd/action/fermion/ZMobiusFermion.h>
NAMESPACE_CHECK(DomainWall);

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
NAMESPACE_CHECK(Overlap);
///////////////////////////////////////////////////////////////////////////////
// G5 herm -- this has to live in QCD since dirac matrix is not in the broader sector of code
///////////////////////////////////////////////////////////////////////////////
#include <Grid/qcd/action/fermion/g5HermitianLinop.h>

///////////////////////////////////////////////////////////////////////////////
// Fourier accelerated Pauli Villars inverse support
///////////////////////////////////////////////////////////////////////////////
#include <Grid/qcd/action/fermion/WilsonTMFermion5D.h>   
NAMESPACE_CHECK(WilsonTM5);

////////////////////////////////////////////////////////////////////////////////
// Move this group to a DWF specific tools/algorithms subdir? 
////////////////////////////////////////////////////////////////////////////////
#include <Grid/qcd/action/fermion/SchurDiagTwoKappa.h>
#include <Grid/qcd/action/fermion/FourierAcceleratedPV.h>
#include <Grid/qcd/action/fermion/PauliVillarsInverters.h>
#include <Grid/qcd/action/fermion/Reconstruct5Dprop.h>
#include <Grid/qcd/action/fermion/MADWF.h>
NAMESPACE_CHECK(DWFutils);

////////////////////////////////////////////////////////////////////////////////////////////////////
// More maintainable to maintain the following typedef list centrally, as more "impl" targets
// are added, (e.g. extension for gparity, half precision project in comms etc..)
////////////////////////////////////////////////////////////////////////////////////////////////////

// Cayley 5d
NAMESPACE_BEGIN(Grid);

typedef WilsonFermion<WilsonImplR> WilsonFermionR;
typedef WilsonFermion<WilsonImplF> WilsonFermionF;
typedef WilsonFermion<WilsonImplD> WilsonFermionD;

typedef WilsonFermion<WilsonImplRL> WilsonFermionRL;
typedef WilsonFermion<WilsonImplFH> WilsonFermionFH;
typedef WilsonFermion<WilsonImplDF> WilsonFermionDF;

typedef WilsonFermion<WilsonAdjImplR> WilsonAdjFermionR;
typedef WilsonFermion<WilsonAdjImplF> WilsonAdjFermionF;
typedef WilsonFermion<WilsonAdjImplD> WilsonAdjFermionD;

typedef WilsonFermion<WilsonTwoIndexSymmetricImplR> WilsonTwoIndexSymmetricFermionR;
typedef WilsonFermion<WilsonTwoIndexSymmetricImplF> WilsonTwoIndexSymmetricFermionF;
typedef WilsonFermion<WilsonTwoIndexSymmetricImplD> WilsonTwoIndexSymmetricFermionD;

typedef WilsonFermion<WilsonTwoIndexAntiSymmetricImplR> WilsonTwoIndexAntiSymmetricFermionR;
typedef WilsonFermion<WilsonTwoIndexAntiSymmetricImplF> WilsonTwoIndexAntiSymmetricFermionF;
typedef WilsonFermion<WilsonTwoIndexAntiSymmetricImplD> WilsonTwoIndexAntiSymmetricFermionD;

// Twisted mass fermion
typedef WilsonTMFermion<WilsonImplR> WilsonTMFermionR;
typedef WilsonTMFermion<WilsonImplF> WilsonTMFermionF;
typedef WilsonTMFermion<WilsonImplD> WilsonTMFermionD;

// Clover fermions
typedef WilsonCloverFermion<WilsonImplR> WilsonCloverFermionR;
typedef WilsonCloverFermion<WilsonImplF> WilsonCloverFermionF;
typedef WilsonCloverFermion<WilsonImplD> WilsonCloverFermionD;

typedef WilsonCloverFermion<WilsonAdjImplR> WilsonCloverAdjFermionR;
typedef WilsonCloverFermion<WilsonAdjImplF> WilsonCloverAdjFermionF;
typedef WilsonCloverFermion<WilsonAdjImplD> WilsonCloverAdjFermionD;

typedef WilsonCloverFermion<WilsonTwoIndexSymmetricImplR> WilsonCloverTwoIndexSymmetricFermionR;
typedef WilsonCloverFermion<WilsonTwoIndexSymmetricImplF> WilsonCloverTwoIndexSymmetricFermionF;
typedef WilsonCloverFermion<WilsonTwoIndexSymmetricImplD> WilsonCloverTwoIndexSymmetricFermionD;

typedef WilsonCloverFermion<WilsonTwoIndexAntiSymmetricImplR> WilsonCloverTwoIndexAntiSymmetricFermionR;
typedef WilsonCloverFermion<WilsonTwoIndexAntiSymmetricImplF> WilsonCloverTwoIndexAntiSymmetricFermionF;
typedef WilsonCloverFermion<WilsonTwoIndexAntiSymmetricImplD> WilsonCloverTwoIndexAntiSymmetricFermionD;

// Domain Wall fermions
typedef DomainWallFermion<WilsonImplR> DomainWallFermionR;
typedef DomainWallFermion<WilsonImplF> DomainWallFermionF;
typedef DomainWallFermion<WilsonImplD> DomainWallFermionD;

typedef DomainWallFermion<WilsonImplRL> DomainWallFermionRL;
typedef DomainWallFermion<WilsonImplFH> DomainWallFermionFH;
typedef DomainWallFermion<WilsonImplDF> DomainWallFermionDF;

typedef DomainWallEOFAFermion<WilsonImplR> DomainWallEOFAFermionR;
typedef DomainWallEOFAFermion<WilsonImplF> DomainWallEOFAFermionF;
typedef DomainWallEOFAFermion<WilsonImplD> DomainWallEOFAFermionD;

typedef DomainWallEOFAFermion<WilsonImplRL> DomainWallEOFAFermionRL;
typedef DomainWallEOFAFermion<WilsonImplFH> DomainWallEOFAFermionFH;
typedef DomainWallEOFAFermion<WilsonImplDF> DomainWallEOFAFermionDF;

typedef MobiusFermion<WilsonImplR> MobiusFermionR;
typedef MobiusFermion<WilsonImplF> MobiusFermionF;
typedef MobiusFermion<WilsonImplD> MobiusFermionD;

typedef MobiusFermion<WilsonImplRL> MobiusFermionRL;
typedef MobiusFermion<WilsonImplFH> MobiusFermionFH;
typedef MobiusFermion<WilsonImplDF> MobiusFermionDF;

typedef MobiusEOFAFermion<WilsonImplR> MobiusEOFAFermionR;
typedef MobiusEOFAFermion<WilsonImplF> MobiusEOFAFermionF;
typedef MobiusEOFAFermion<WilsonImplD> MobiusEOFAFermionD;

typedef MobiusEOFAFermion<WilsonImplRL> MobiusEOFAFermionRL;
typedef MobiusEOFAFermion<WilsonImplFH> MobiusEOFAFermionFH;
typedef MobiusEOFAFermion<WilsonImplDF> MobiusEOFAFermionDF;

typedef ZMobiusFermion<ZWilsonImplR> ZMobiusFermionR;
typedef ZMobiusFermion<ZWilsonImplF> ZMobiusFermionF;
typedef ZMobiusFermion<ZWilsonImplD> ZMobiusFermionD;

typedef ZMobiusFermion<ZWilsonImplRL> ZMobiusFermionRL;
typedef ZMobiusFermion<ZWilsonImplFH> ZMobiusFermionFH;
typedef ZMobiusFermion<ZWilsonImplDF> ZMobiusFermionDF;

// Ls vectorised
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

typedef WilsonFermion<GparityWilsonImplRL>     GparityWilsonFermionRL;
typedef WilsonFermion<GparityWilsonImplFH>     GparityWilsonFermionFH;
typedef WilsonFermion<GparityWilsonImplDF>     GparityWilsonFermionDF;

typedef DomainWallFermion<GparityWilsonImplR> GparityDomainWallFermionR;
typedef DomainWallFermion<GparityWilsonImplF> GparityDomainWallFermionF;
typedef DomainWallFermion<GparityWilsonImplD> GparityDomainWallFermionD;

typedef DomainWallFermion<GparityWilsonImplRL> GparityDomainWallFermionRL;
typedef DomainWallFermion<GparityWilsonImplFH> GparityDomainWallFermionFH;
typedef DomainWallFermion<GparityWilsonImplDF> GparityDomainWallFermionDF;

typedef DomainWallEOFAFermion<GparityWilsonImplR> GparityDomainWallEOFAFermionR;
typedef DomainWallEOFAFermion<GparityWilsonImplF> GparityDomainWallEOFAFermionF;
typedef DomainWallEOFAFermion<GparityWilsonImplD> GparityDomainWallEOFAFermionD;

typedef DomainWallEOFAFermion<GparityWilsonImplRL> GparityDomainWallEOFAFermionRL;
typedef DomainWallEOFAFermion<GparityWilsonImplFH> GparityDomainWallEOFAFermionFH;
typedef DomainWallEOFAFermion<GparityWilsonImplDF> GparityDomainWallEOFAFermionDF;

typedef WilsonTMFermion<GparityWilsonImplR> GparityWilsonTMFermionR;
typedef WilsonTMFermion<GparityWilsonImplF> GparityWilsonTMFermionF;
typedef WilsonTMFermion<GparityWilsonImplD> GparityWilsonTMFermionD;

typedef WilsonTMFermion<GparityWilsonImplRL> GparityWilsonTMFermionRL;
typedef WilsonTMFermion<GparityWilsonImplFH> GparityWilsonTMFermionFH;
typedef WilsonTMFermion<GparityWilsonImplDF> GparityWilsonTMFermionDF;

typedef MobiusFermion<GparityWilsonImplR> GparityMobiusFermionR;
typedef MobiusFermion<GparityWilsonImplF> GparityMobiusFermionF;
typedef MobiusFermion<GparityWilsonImplD> GparityMobiusFermionD;

typedef MobiusFermion<GparityWilsonImplRL> GparityMobiusFermionRL;
typedef MobiusFermion<GparityWilsonImplFH> GparityMobiusFermionFH;
typedef MobiusFermion<GparityWilsonImplDF> GparityMobiusFermionDF;

typedef MobiusEOFAFermion<GparityWilsonImplR> GparityMobiusEOFAFermionR;
typedef MobiusEOFAFermion<GparityWilsonImplF> GparityMobiusEOFAFermionF;
typedef MobiusEOFAFermion<GparityWilsonImplD> GparityMobiusEOFAFermionD;

typedef MobiusEOFAFermion<GparityWilsonImplRL> GparityMobiusEOFAFermionRL;
typedef MobiusEOFAFermion<GparityWilsonImplFH> GparityMobiusEOFAFermionFH;
typedef MobiusEOFAFermion<GparityWilsonImplDF> GparityMobiusEOFAFermionDF;

typedef ImprovedStaggeredFermion<StaggeredImplR> ImprovedStaggeredFermionR;
typedef ImprovedStaggeredFermion<StaggeredImplF> ImprovedStaggeredFermionF;
typedef ImprovedStaggeredFermion<StaggeredImplD> ImprovedStaggeredFermionD;

typedef NaiveStaggeredFermion<StaggeredImplR> NaiveStaggeredFermionR;
typedef NaiveStaggeredFermion<StaggeredImplF> NaiveStaggeredFermionF;
typedef NaiveStaggeredFermion<StaggeredImplD> NaiveStaggeredFermionD;

typedef ImprovedStaggeredFermion5D<StaggeredImplR> ImprovedStaggeredFermion5DR;
typedef ImprovedStaggeredFermion5D<StaggeredImplF> ImprovedStaggeredFermion5DF;
typedef ImprovedStaggeredFermion5D<StaggeredImplD> ImprovedStaggeredFermion5DD;

NAMESPACE_END(Grid);

////////////////////
// Scalar QED actions
// TODO: this needs to move to another header after rename to Fermion.h
////////////////////
#include <Grid/qcd/action/scalar/Scalar.h>
#include <Grid/qcd/action/gauge/Photon.h>


