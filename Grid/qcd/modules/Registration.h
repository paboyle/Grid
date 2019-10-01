/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/modules/Registration.h

Copyright (C) 2016

Author: Guido Cossu <guido.cossu@ed.ac.uk>

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

#ifndef MODULES_REGISTRATION_H
#define MODULES_REGISTRATION_H

			   // simplify with macros

			   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			   // Actions
			   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef WilsonGModule<ImplementationPolicy> WilsonGMod;
typedef SymanzikGModule<ImplementationPolicy> SymanzikGMod;
typedef IwasakiGModule<ImplementationPolicy> IwasakiGMod;
typedef DBW2GModule<ImplementationPolicy> DBW2GMod;
typedef RBCGModule<ImplementationPolicy> RBCGMod;
typedef PlaqPlusRectangleGModule<ImplementationPolicy> PlaqPlusRectangleGMod;

static Registrar<WilsonGMod,            HMC_ActionModuleFactory<gauge_string, typename ImplementationPolicy::Field, Serialiser> > __WGmodXMLInit("Wilson"); 
static Registrar<SymanzikGMod,          HMC_ActionModuleFactory<gauge_string, typename ImplementationPolicy::Field, Serialiser> > __SymGmodXMLInit("Symanzik"); 
static Registrar<IwasakiGMod,           HMC_ActionModuleFactory<gauge_string, typename ImplementationPolicy::Field, Serialiser> > __IwGmodXMLInit("Iwasaki"); 
static Registrar<DBW2GMod,              HMC_ActionModuleFactory<gauge_string, typename ImplementationPolicy::Field, Serialiser> > __DBW2GmodXMLInit("DBW2"); 
static Registrar<RBCGMod,               HMC_ActionModuleFactory<gauge_string, typename ImplementationPolicy::Field, Serialiser> > __RBCGmodXMLInit("RBC"); 
static Registrar<PlaqPlusRectangleGMod, HMC_ActionModuleFactory<gauge_string, typename ImplementationPolicy::Field, Serialiser> > __PPRectGmodXMLInit("PlaqPlusRect"); 


// FIXME more general implementation
static Registrar<TwoFlavourFModule<FermionImplementationPolicy> ,      
		 HMC_ActionModuleFactory<gauge_string, typename ImplementationPolicy::Field, Serialiser> > __TwoFlavourFmodXMLInit("TwoFlavours"); 
static Registrar<TwoFlavourRatioFModule<FermionImplementationPolicy> , 
		 HMC_ActionModuleFactory<gauge_string, typename ImplementationPolicy::Field, Serialiser> > __TwoFlavourRatioFmodXMLInit("TwoFlavoursRatio"); 
static Registrar<TwoFlavourEOFModule<FermionImplementationPolicy> ,    
		 HMC_ActionModuleFactory<gauge_string, typename ImplementationPolicy::Field, Serialiser> > __TwoFlavourEOFmodXMLInit("TwoFlavoursEvenOdd"); 
static Registrar<TwoFlavourRatioEOFModule<FermionImplementationPolicy>,
		 HMC_ActionModuleFactory<gauge_string, typename ImplementationPolicy::Field, Serialiser> > __TwoFlavourRatioEOFmodXMLInit("TwoFlavoursEvenOddRatio"); 
static Registrar<OneFlavourFModule<FermionImplementationPolicy> ,      
		 HMC_ActionModuleFactory<gauge_string, typename ImplementationPolicy::Field, Serialiser> > __OneFlavourFmodXMLInit("OneFlavour"); 
static Registrar<OneFlavourEOFModule<FermionImplementationPolicy> ,    
		 HMC_ActionModuleFactory<gauge_string, typename ImplementationPolicy::Field, Serialiser> > __OneFlavourEOFmodXMLInit("OneFlavourEvenOdd"); 
static Registrar<OneFlavourRatioFModule<FermionImplementationPolicy> , 
		 HMC_ActionModuleFactory<gauge_string, typename ImplementationPolicy::Field, Serialiser> > __OneFlavourRatioFmodXMLInit("OneFlavourRatio"); 
static Registrar<OneFlavourRatioEOFModule<FermionImplementationPolicy>,
		 HMC_ActionModuleFactory<gauge_string, typename ImplementationPolicy::Field, Serialiser> > __OneFlavourRatioEOFmodXMLInit("OneFlavourEvenOddRatio"); 



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Solvers
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Now a specific registration with a fermion field
// here must instantiate CG and CR for every new fermion field type (macro!!)

static Registrar< ConjugateGradientModule<WilsonFermionR::FermionField>,   
                  HMC_SolverModuleFactory<solver_string, WilsonFermionR::FermionField, Serialiser> > __CGWFmodXMLInit("ConjugateGradient"); 
static Registrar< ConjugateResidualModule<WilsonFermionR::FermionField>,   
                  HMC_SolverModuleFactory<solver_string, WilsonFermionR::FermionField, Serialiser> > __CRWFmodXMLInit("ConjugateResidual"); 

// add the staggered, scalar versions here


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fermion operators
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static Registrar< WilsonFermionModule<FermionImplementationPolicy>,   
                  HMC_FermionOperatorModuleFactory<fermionop_string, FermionImplementationPolicy, Serialiser> > __WilsonFOPmodXMLInit("Wilson"); 
static Registrar< MobiusFermionModule<FermionImplementationPolicy>,   
                  HMC_FermionOperatorModuleFactory<fermionop_string, FermionImplementationPolicy, Serialiser> > __MobiusFOPmodXMLInit("Mobius");
static Registrar< DomainWallFermionModule<FermionImplementationPolicy>,   
                  HMC_FermionOperatorModuleFactory<fermionop_string, FermionImplementationPolicy, Serialiser> > __DWFOPmodXMLInit("DomainWall");


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Observables
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static Registrar<PlaquetteMod<ImplementationPolicy>, HMC_ObservablesModuleFactory<observable_string, typename ImplementationPolicy::Field, Serialiser> > __OBSPLmodXMLInit("Plaquette"); 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Checkpointers
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static Registrar<BinaryCPModule<ImplementationPolicy>, HMC_CPModuleFactory<cp_string, ImplementationPolicy, Serialiser> > __CPBinarymodXMLInit("Binary");
static Registrar<NerscCPModule<ImplementationPolicy> , HMC_CPModuleFactory<cp_string, ImplementationPolicy, Serialiser> > __CPNerscmodXMLInit("Nersc");

#ifdef HAVE_LIME
static Registrar<ILDGCPModule<ImplementationPolicy>  , HMC_CPModuleFactory<cp_string, ImplementationPolicy, Serialiser> > __CPILDGmodXMLInit("ILDG");
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Integrators
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static Registrar< HMCLeapFrog<ImplementationPolicy, RepresentationPolicy, Serialiser>      , HMCRunnerModuleFactory<hmc_string, Serialiser> > __HMCLFmodXMLInit("LeapFrog");
static Registrar< HMCMinimumNorm2<ImplementationPolicy, RepresentationPolicy, Serialiser>  , HMCRunnerModuleFactory<hmc_string, Serialiser> > __HMCMN2modXMLInit("MinimumNorm2");
static Registrar< HMCForceGradient<ImplementationPolicy, RepresentationPolicy, Serialiser> , HMCRunnerModuleFactory<hmc_string, Serialiser> > __HMCFGmodXMLInit("ForceGradient");

typedef HMCRunnerModuleFactory<hmc_string, Serialiser > HMCModuleFactory;

#endif
