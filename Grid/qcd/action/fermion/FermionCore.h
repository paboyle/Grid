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
#ifndef  GRID_QCD_FERMION_CORE_H
#define  GRID_QCD_FERMION_CORE_H

#include <Grid/GridCore.h>
#include <Grid/GridQCDcore.h>
#include <Grid/qcd/action/ActionCore.h>

////////////////////////////////////////////
// Fermion prereqs
////////////////////////////////////////////
#include <Grid/qcd/action/fermion/WilsonCompressor.h>     //used by all wilson type fermions
#include <Grid/qcd/action/fermion/FermionOperatorImpl.h>
#include <Grid/qcd/action/fermion/FermionOperator.h>
#include <Grid/qcd/action/fermion/WilsonKernels.h>        //used by all wilson type fermions
#include <Grid/qcd/action/fermion/StaggeredKernels.h>        //used by all wilson type fermions

#define FermOpStaggeredTemplateInstantiate(A) \
  template class A<StaggeredImplF>; \
  template class A<StaggeredImplD>; 

#define FermOpStaggeredVec5dTemplateInstantiate(A) \
  template class A<StaggeredVec5dImplF>; \
  template class A<StaggeredVec5dImplD>; 

#define FermOp4dVecTemplateInstantiate(A) \
  template class A<WilsonImplF>;		\
  template class A<WilsonImplD>;		\
  template class A<ZWilsonImplF>;		\
  template class A<ZWilsonImplD>;		\
  template class A<GparityWilsonImplF>;		\
  template class A<GparityWilsonImplD>;		\
  template class A<WilsonImplFH>;		\
  template class A<WilsonImplDF>;		\
  template class A<ZWilsonImplFH>;		\
  template class A<ZWilsonImplDF>;		\
  template class A<GparityWilsonImplFH>;		\
  template class A<GparityWilsonImplDF>;		


#define AdjointFermOpTemplateInstantiate(A) \
  template class A<WilsonAdjImplF>; \
  template class A<WilsonAdjImplD>; 

#define TwoIndexFermOpTemplateInstantiate(A) \
  template class A<WilsonTwoIndexSymmetricImplF>; \
  template class A<WilsonTwoIndexSymmetricImplD>; \
  template class A<WilsonTwoIndexAntiSymmetricImplF>; \
  template class A<WilsonTwoIndexAntiSymmetricImplD>;

#define FermOp5dVecTemplateInstantiate(A) \
  template class A<DomainWallVec5dImplF>;	\
  template class A<DomainWallVec5dImplD>;	\
  template class A<ZDomainWallVec5dImplF>;	\
  template class A<ZDomainWallVec5dImplD>;	\
  template class A<DomainWallVec5dImplFH>;	\
  template class A<DomainWallVec5dImplDF>;	\
  template class A<ZDomainWallVec5dImplFH>;	\
  template class A<ZDomainWallVec5dImplDF>;	

#define FermOpTemplateInstantiate(A) \
 FermOp4dVecTemplateInstantiate(A) \
 FermOp5dVecTemplateInstantiate(A) 

#define GparityFermOpTemplateInstantiate(A) 

#endif
