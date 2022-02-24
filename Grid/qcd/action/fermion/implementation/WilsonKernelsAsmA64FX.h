/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid



    Source file: ./lib/qcd/action/fermion/WilsonKernelsAsmA64FX.h

    Copyright (C) 2020

Author: Nils Meyer  <nils.meyer@ur.de>  Regensburg University

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

//#if defined(A64FXASM)
#if defined(A64FX)

// safety include
#include <arm_sve.h>

// undefine everything related to kernels
#include <simd/Fujitsu_A64FX_undef.h>


    ///////////////////////////////////////////////////////////
    // If we are A64FX specialise the single precision routine
    ///////////////////////////////////////////////////////////
#if defined(DSLASHINTRIN)
//#pragma message ("A64FX Dslash: intrin")
#include <simd/Fujitsu_A64FX_intrin_single.h>
#else
#pragma message ("A64FX Dslash: asm")
#include <simd/Fujitsu_A64FX_asm_single.h>
#endif

/// Switch off the 5d vectorised code optimisations
#undef DWFVEC5D

/////////////////////////////////////////////////////////////////
// XYZT vectorised, undag Kernel, single
/////////////////////////////////////////////////////////////////
#undef KERNEL_DAG
#define INTERIOR_AND_EXTERIOR
#undef INTERIOR
#undef EXTERIOR

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<WilsonImplF>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<ZWilsonImplF>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

//#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
//template<> void
//WilsonKernels<WilsonImplFH>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
//						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
//#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

//#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
//template<> void
//WilsonKernels<ZWilsonImplFH>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
//						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
//#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>


#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<WilsonImplF>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<ZWilsonImplF>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

//#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
//template<> void
//WilsonKernels<WilsonImplFH>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
//						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
//#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

//#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
//template<> void
//WilsonKernels<ZWilsonImplFH>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
//						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
//#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>


#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<WilsonImplF>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<ZWilsonImplF>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

//#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
//template<> void
//WilsonKernels<WilsonImplFH>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
//						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
//#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

//#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
//template<> void
//WilsonKernels<ZWilsonImplFH>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
//						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
//#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>



/////////////////////////////////////////////////////////////////
// XYZT vectorised, dag Kernel, single
/////////////////////////////////////////////////////////////////
#define KERNEL_DAG
#define INTERIOR_AND_EXTERIOR
#undef INTERIOR
#undef EXTERIOR

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<WilsonImplF>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<ZWilsonImplF>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

//#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
//template<> void
//WilsonKernels<WilsonImplFH>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
//						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
//#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

//#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
//template<> void
//WilsonKernels<ZWilsonImplFH>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
//						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
//#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>


#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<WilsonImplF>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<ZWilsonImplF>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

//#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
//template<> void
//WilsonKernels<WilsonImplFH>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
//						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
//#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

//#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
//template<> void
//WilsonKernels<ZWilsonImplFH>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
//						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
//#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>


#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<WilsonImplF>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<ZWilsonImplF>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

//#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
//template<> void
//WilsonKernels<WilsonImplFH>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
//						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
//#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

//#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
//template<> void
//WilsonKernels<ZWilsonImplFH>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
//						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
//#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>



// undefine
#include <simd/Fujitsu_A64FX_undef.h>

///////////////////////////////////////////////////////////
// If we are A64FX specialise the double precision routine
///////////////////////////////////////////////////////////

#if defined(DSLASHINTRIN)
#include <simd/Fujitsu_A64FX_intrin_double.h>
#else
#include <simd/Fujitsu_A64FX_asm_double.h>
#endif

// former KNL
//#define MAYBEPERM(A,perm) if (perm) { A ; }
//#define MULT_2SPIN(ptr,pf) MULT_ADDSUB_2SPIN(ptr,pf)
//#define COMPLEX_SIGNS(isigns) vComplexD *isigns = &signsD[0];


#define INTERIOR_AND_EXTERIOR
#undef  INTERIOR
#undef  EXTERIOR

/////////////////////////////////////////////////////////////////
// XYZT vectorised, undag Kernel, double
/////////////////////////////////////////////////////////////////
#undef KERNEL_DAG
#define INTERIOR_AND_EXTERIOR
#undef INTERIOR
#undef EXTERIOR

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<WilsonImplD>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<ZWilsonImplD>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

// #pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
// template<> void
// WilsonKernels<WilsonImplDF>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
// 						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
// #include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

// #pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
// template<> void
// WilsonKernels<ZWilsonImplDF>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
// 						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
// #include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>


#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<WilsonImplD>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<ZWilsonImplD>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

// #pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
// template<> void
// WilsonKernels<WilsonImplDF>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
// 						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
// #include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

// #pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
// template<> void
// WilsonKernels<ZWilsonImplDF>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
// 						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
// #include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>


#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<WilsonImplD>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<ZWilsonImplD>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

// #pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
// template<> void
// WilsonKernels<WilsonImplDF>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
// 						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
// #include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

// #pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
// template<> void
// WilsonKernels<ZWilsonImplDF>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
// 						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
// #include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>


/////////////////////////////////////////////////////////////////
// XYZT vectorised, dag Kernel, double
/////////////////////////////////////////////////////////////////
#define KERNEL_DAG
#define INTERIOR_AND_EXTERIOR
#undef INTERIOR
#undef EXTERIOR

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<WilsonImplD>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<ZWilsonImplD>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

// #pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
// template<> void
// WilsonKernels<WilsonImplDF>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
// 						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
// #include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

// #pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
// template<> void
// WilsonKernels<ZWilsonImplDF>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
// 						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
// #include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>


#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<WilsonImplD>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<ZWilsonImplD>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

// #pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
// template<> void
// WilsonKernels<WilsonImplDF>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
// 						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
// #include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

// #pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
// template<> void
// WilsonKernels<ZWilsonImplDF>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
// 						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
// #include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>


#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<WilsonImplD>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

#pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
template<> void
WilsonKernels<ZWilsonImplD>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

// #pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
// template<> void
// WilsonKernels<WilsonImplDF>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
// 						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
// #include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>

// #pragma GCC optimize ("-O3", "-fno-schedule-insns", "-fno-schedule-insns2")
// template<> void
// WilsonKernels<ZWilsonImplDF>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
// 						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
// #include <qcd/action/fermion/implementation/WilsonKernelsAsmBodyA64FX.h>




// undefs
#include <simd/Fujitsu_A64FX_undef.h>

#endif //A64FXASM
