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
#include <Grid/qcd/action/fermion/DWFSlow.h>       // Slow DWF

#include <Grid/qcd/action/fermion/WilsonFermion.h>       // 4d wilson like
NAMESPACE_CHECK(Wilson);
#include <Grid/qcd/action/fermion/WilsonTMFermion.h>       // 4d wilson like
NAMESPACE_CHECK(WilsonTM);
#include <Grid/qcd/action/fermion/WilsonCloverFermion.h> // 4d wilson clover fermions
#include <Grid/qcd/action/fermion/CompactWilsonCloverFermion.h> // 4d compact wilson clover fermions
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

typedef WilsonFermion<WilsonImplD2> WilsonFermionD2;
typedef WilsonFermion<WilsonImplF> WilsonFermionF;
typedef WilsonFermion<WilsonImplD> WilsonFermionD;

typedef WilsonFermion<WilsonAdjImplF> WilsonAdjFermionF;
typedef WilsonFermion<WilsonAdjImplD> WilsonAdjFermionD;

typedef WilsonFermion<WilsonTwoIndexSymmetricImplF> WilsonTwoIndexSymmetricFermionF;
typedef WilsonFermion<WilsonTwoIndexSymmetricImplD> WilsonTwoIndexSymmetricFermionD;

typedef WilsonFermion<WilsonTwoIndexAntiSymmetricImplF> WilsonTwoIndexAntiSymmetricFermionF;
typedef WilsonFermion<WilsonTwoIndexAntiSymmetricImplD> WilsonTwoIndexAntiSymmetricFermionD;

// Twisted mass fermion
typedef WilsonTMFermion<WilsonImplD2> WilsonTMFermionD2;
typedef WilsonTMFermion<WilsonImplF> WilsonTMFermionF;
typedef WilsonTMFermion<WilsonImplD> WilsonTMFermionD;

// Clover fermions
template <typename WImpl> using WilsonClover = WilsonCloverFermion<WImpl, CloverHelpers<WImpl>>;
template <typename WImpl> using WilsonExpClover = WilsonCloverFermion<WImpl, ExpCloverHelpers<WImpl>>;

typedef WilsonClover<WilsonImplD2> WilsonCloverFermionD2;
typedef WilsonClover<WilsonImplF> WilsonCloverFermionF;
typedef WilsonClover<WilsonImplD> WilsonCloverFermionD;

typedef WilsonExpClover<WilsonImplD2> WilsonExpCloverFermionD2;
typedef WilsonExpClover<WilsonImplF> WilsonExpCloverFermionF;
typedef WilsonExpClover<WilsonImplD> WilsonExpCloverFermionD;

typedef WilsonClover<WilsonAdjImplF> WilsonCloverAdjFermionF;
typedef WilsonClover<WilsonAdjImplD> WilsonCloverAdjFermionD;

typedef WilsonClover<WilsonTwoIndexSymmetricImplF> WilsonCloverTwoIndexSymmetricFermionF;
typedef WilsonClover<WilsonTwoIndexSymmetricImplD> WilsonCloverTwoIndexSymmetricFermionD;

typedef WilsonClover<WilsonTwoIndexAntiSymmetricImplF> WilsonCloverTwoIndexAntiSymmetricFermionF;
typedef WilsonClover<WilsonTwoIndexAntiSymmetricImplD> WilsonCloverTwoIndexAntiSymmetricFermionD;

// Compact Clover fermions
template <typename WImpl> using CompactWilsonClover = CompactWilsonCloverFermion<WImpl, CompactCloverHelpers<WImpl>>;
template <typename WImpl> using CompactWilsonExpClover = CompactWilsonCloverFermion<WImpl, CompactExpCloverHelpers<WImpl>>;

typedef CompactWilsonClover<WilsonImplD2> CompactWilsonCloverFermionD2;
typedef CompactWilsonClover<WilsonImplF> CompactWilsonCloverFermionF;
typedef CompactWilsonClover<WilsonImplD> CompactWilsonCloverFermionD;

typedef CompactWilsonExpClover<WilsonImplD2> CompactWilsonExpCloverFermionD2;
typedef CompactWilsonExpClover<WilsonImplF> CompactWilsonExpCloverFermionF;
typedef CompactWilsonExpClover<WilsonImplD> CompactWilsonExpCloverFermionD;

typedef CompactWilsonClover<WilsonAdjImplF> CompactWilsonCloverAdjFermionF;
typedef CompactWilsonClover<WilsonAdjImplD> CompactWilsonCloverAdjFermionD;

typedef CompactWilsonClover<WilsonTwoIndexSymmetricImplF> CompactWilsonCloverTwoIndexSymmetricFermionF;
typedef CompactWilsonClover<WilsonTwoIndexSymmetricImplD> CompactWilsonCloverTwoIndexSymmetricFermionD;

typedef CompactWilsonClover<WilsonTwoIndexAntiSymmetricImplF> CompactWilsonCloverTwoIndexAntiSymmetricFermionF;
typedef CompactWilsonClover<WilsonTwoIndexAntiSymmetricImplD> CompactWilsonCloverTwoIndexAntiSymmetricFermionD;

// Domain Wall fermions
typedef DomainWallFermion<WilsonImplF> DomainWallFermionF;
typedef DomainWallFermion<WilsonImplD> DomainWallFermionD;
typedef DomainWallFermion<WilsonImplD2> DomainWallFermionD2;

typedef DomainWallEOFAFermion<WilsonImplD2> DomainWallEOFAFermionD2;
typedef DomainWallEOFAFermion<WilsonImplF> DomainWallEOFAFermionF;
typedef DomainWallEOFAFermion<WilsonImplD> DomainWallEOFAFermionD;

typedef MobiusFermion<WilsonImplD2> MobiusFermionD2;
typedef MobiusFermion<WilsonImplF> MobiusFermionF;
typedef MobiusFermion<WilsonImplD> MobiusFermionD;

typedef MobiusEOFAFermion<WilsonImplD2> MobiusEOFAFermionD2;
typedef MobiusEOFAFermion<WilsonImplF> MobiusEOFAFermionF;
typedef MobiusEOFAFermion<WilsonImplD> MobiusEOFAFermionD;

typedef ZMobiusFermion<ZWilsonImplD2> ZMobiusFermionD2;
typedef ZMobiusFermion<ZWilsonImplF> ZMobiusFermionF;
typedef ZMobiusFermion<ZWilsonImplD> ZMobiusFermionD;

typedef ScaledShamirFermion<WilsonImplD2> ScaledShamirFermionD2;
typedef ScaledShamirFermion<WilsonImplF> ScaledShamirFermionF;
typedef ScaledShamirFermion<WilsonImplD> ScaledShamirFermionD;

typedef MobiusZolotarevFermion<WilsonImplD2> MobiusZolotarevFermionD2;
typedef MobiusZolotarevFermion<WilsonImplF> MobiusZolotarevFermionF;
typedef MobiusZolotarevFermion<WilsonImplD> MobiusZolotarevFermionD;
typedef ShamirZolotarevFermion<WilsonImplD2> ShamirZolotarevFermionD2;
typedef ShamirZolotarevFermion<WilsonImplF> ShamirZolotarevFermionF;
typedef ShamirZolotarevFermion<WilsonImplD> ShamirZolotarevFermionD;

typedef OverlapWilsonCayleyTanhFermion<WilsonImplD2> OverlapWilsonCayleyTanhFermionD2;
typedef OverlapWilsonCayleyTanhFermion<WilsonImplF> OverlapWilsonCayleyTanhFermionF;
typedef OverlapWilsonCayleyTanhFermion<WilsonImplD> OverlapWilsonCayleyTanhFermionD;
typedef OverlapWilsonCayleyZolotarevFermion<WilsonImplD2> OverlapWilsonCayleyZolotarevFermionD2;
typedef OverlapWilsonCayleyZolotarevFermion<WilsonImplF> OverlapWilsonCayleyZolotarevFermionF;
typedef OverlapWilsonCayleyZolotarevFermion<WilsonImplD> OverlapWilsonCayleyZolotarevFermionD;

// Continued fraction
typedef OverlapWilsonContFracTanhFermion<WilsonImplD2> OverlapWilsonContFracTanhFermionD2;
typedef OverlapWilsonContFracTanhFermion<WilsonImplF> OverlapWilsonContFracTanhFermionF;
typedef OverlapWilsonContFracTanhFermion<WilsonImplD> OverlapWilsonContFracTanhFermionD;
typedef OverlapWilsonContFracZolotarevFermion<WilsonImplD2> OverlapWilsonContFracZolotarevFermionD2;
typedef OverlapWilsonContFracZolotarevFermion<WilsonImplF> OverlapWilsonContFracZolotarevFermionF;
typedef OverlapWilsonContFracZolotarevFermion<WilsonImplD> OverlapWilsonContFracZolotarevFermionD;

// Partial fraction
typedef OverlapWilsonPartialFractionTanhFermion<WilsonImplD2> OverlapWilsonPartialFractionTanhFermionD2;
typedef OverlapWilsonPartialFractionTanhFermion<WilsonImplF> OverlapWilsonPartialFractionTanhFermionF;
typedef OverlapWilsonPartialFractionTanhFermion<WilsonImplD> OverlapWilsonPartialFractionTanhFermionD;

typedef OverlapWilsonPartialFractionZolotarevFermion<WilsonImplD2> OverlapWilsonPartialFractionZolotarevFermionD2;
typedef OverlapWilsonPartialFractionZolotarevFermion<WilsonImplF> OverlapWilsonPartialFractionZolotarevFermionF;
typedef OverlapWilsonPartialFractionZolotarevFermion<WilsonImplD> OverlapWilsonPartialFractionZolotarevFermionD;

// Gparity cases; partial list until tested
typedef WilsonFermion<GparityWilsonImplF>     GparityWilsonFermionF;
typedef WilsonFermion<GparityWilsonImplD>     GparityWilsonFermionD;

typedef DomainWallFermion<GparityWilsonImplF> GparityDomainWallFermionF;
typedef DomainWallFermion<GparityWilsonImplD> GparityDomainWallFermionD;

typedef DomainWallEOFAFermion<GparityWilsonImplR> GparityDomainWallEOFAFermionD2;
typedef DomainWallEOFAFermion<GparityWilsonImplF> GparityDomainWallEOFAFermionF;
typedef DomainWallEOFAFermion<GparityWilsonImplD> GparityDomainWallEOFAFermionD;

typedef WilsonTMFermion<GparityWilsonImplR> GparityWilsonTMFermionD2;
typedef WilsonTMFermion<GparityWilsonImplF> GparityWilsonTMFermionF;
typedef WilsonTMFermion<GparityWilsonImplD> GparityWilsonTMFermionD;

typedef MobiusFermion<GparityWilsonImplR> GparityMobiusFermionD2;
typedef MobiusFermion<GparityWilsonImplF> GparityMobiusFermionF;
typedef MobiusFermion<GparityWilsonImplD> GparityMobiusFermionD;

typedef MobiusEOFAFermion<GparityWilsonImplR> GparityMobiusEOFAFermionD2;
typedef MobiusEOFAFermion<GparityWilsonImplF> GparityMobiusEOFAFermionF;
typedef MobiusEOFAFermion<GparityWilsonImplD> GparityMobiusEOFAFermionD;

typedef ImprovedStaggeredFermion<StaggeredImplF> ImprovedStaggeredFermionF;
typedef ImprovedStaggeredFermion<StaggeredImplD> ImprovedStaggeredFermionD;

typedef NaiveStaggeredFermion<StaggeredImplF> NaiveStaggeredFermionF;
typedef NaiveStaggeredFermion<StaggeredImplD> NaiveStaggeredFermionD;

typedef ImprovedStaggeredFermion5D<StaggeredImplF> ImprovedStaggeredFermion5DF;
typedef ImprovedStaggeredFermion5D<StaggeredImplD> ImprovedStaggeredFermion5DD;

NAMESPACE_END(Grid);

////////////////////
// Scalar QED actions
// TODO: this needs to move to another header after rename to Fermion.h
////////////////////
#include <Grid/qcd/action/scalar/Scalar.h>
#include <Grid/qcd/action/gauge/Photon.h>


