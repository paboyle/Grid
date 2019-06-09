/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/WilsonKernelsGpu.cc

Copyright (C) 2018

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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

#pragma once

#include <Grid/qcd/action/fermion/FermionCore.h>

NAMESPACE_BEGIN(Grid);


#define KERNEL_CALL(A) \
      const uint64_t nsimd = Simd::Nsimd(); \
      const uint64_t    NN = Nsite*Ls*nsimd;\
      accelerator_loopN( sss, NN, {         \
	  uint64_t cur  = sss;              \
	  cur = cur / nsimd;                \
	  uint64_t   s  = cur%Ls;           \
	  cur = cur / Ls;                   \
	  uint64_t   sU = cur;              \
	  WilsonKernels<Impl>::A(st_v,U_v[sU],buf,Ls,s,sU,in_v,out_v);\
      });
 
#define HOST_CALL(A) \
  const uint64_t nsimd = Simd::Nsimd();					\
  const uint64_t    NN = Nsite*Ls;					\
  SIMT_loop( ss, NN, nsimd, {						\
      int sF = ss;							\
      int sU = ss/Ls;							\
      WilsonKernels<Impl>::A(st_v,U_v,buf,sF,sU,in_v,out_v);	\
  });

#define ASM_CALL(A) \
  SIMT_loop( ss, Nsite, {						\
      int sU = ss;							\
      int sF = ss*Ls;							\
      WilsonKernels<Impl>::A(st_v,U_v,buf,sF,sU,Ls,1,in_v,out_v);	\
  });

template <class Impl>
void WilsonKernels<Impl>::DhopKernel(int Opt,StencilImpl &st,  DoubledGaugeField &U, SiteHalfSpinor * buf,
				     int Ls, int Nsite, const FermionField &in, FermionField &out,
				     int interior,int exterior) 
{
    auto U_v   =   U.View();
    auto in_v  =  in.View();
    auto out_v = out.View();
    auto st_v  =  st.View();

   if( interior && exterior ) { 
     if (Opt == WilsonKernelsStatic::OptGeneric    ) { HOST_CALL(GenericDhopSite); return;}
#ifndef GRID_NVCC
     if (Opt == WilsonKernelsStatic::OptHandUnroll ) { HOST_CALL(HandDhopSite);    return;}
     if (Opt == WilsonKernelsStatic::OptInlineAsm  ) {  ASM_CALL(AsmDhopSite);     return;}
#endif
   } else if( interior ) {
     if (Opt == WilsonKernelsStatic::OptGeneric    ) { HOST_CALL(GenericDhopSiteInt); return;}
#ifndef GRID_NVCC
     if (Opt == WilsonKernelsStatic::OptHandUnroll ) { HOST_CALL(HandDhopSiteInt);    return;}
     if (Opt == WilsonKernelsStatic::OptInlineAsm  ) {  ASM_CALL(AsmDhopSiteInt);     return;}
#endif
   } else if( exterior ) { 
     if (Opt == WilsonKernelsStatic::OptGeneric    ) { HOST_CALL(GenericDhopSiteExt); return;}
#ifndef GRID_NVCC
     if (Opt == WilsonKernelsStatic::OptHandUnroll ) { HOST_CALL(HandDhopSiteExt);    return;}
     if (Opt == WilsonKernelsStatic::OptInlineAsm  ) {  ASM_CALL(AsmDhopSiteExt);     return;}
#endif
   }
   assert(0 && " Kernel optimisation case not covered ");
  }
  template <class Impl>
  void WilsonKernels<Impl>::DhopDagKernel(int Opt,StencilImpl &st,  DoubledGaugeField &U, SiteHalfSpinor * buf,
					  int Ls, int Nsite, const FermionField &in, FermionField &out,
					  int interior,int exterior) 
  {
    auto U_v   = U.View();
    auto in_v  = in.View();
    auto out_v = out.View();
    auto st_v  = st.View();

   if( interior && exterior ) { 
     if (Opt == WilsonKernelsStatic::OptGeneric    ) { HOST_CALL(GenericDhopSiteDag); return;}
#ifndef GRID_NVCC
     if (Opt == WilsonKernelsStatic::OptHandUnroll ) { HOST_CALL(HandDhopSiteDag);    return;}
     if (Opt == WilsonKernelsStatic::OptInlineAsm  ) {  ASM_CALL(AsmDhopSiteDag);     return;}
#endif
   } else if( interior ) {
     if (Opt == WilsonKernelsStatic::OptGeneric    ) { HOST_CALL(GenericDhopSiteDagInt); return;}
#ifndef GRID_NVCC
     if (Opt == WilsonKernelsStatic::OptHandUnroll ) { HOST_CALL(HandDhopSiteDagInt);    return;}
     if (Opt == WilsonKernelsStatic::OptInlineAsm  ) {  ASM_CALL(AsmDhopSiteDagInt);     return;}
#endif
   } else if( exterior ) { 
     if (Opt == WilsonKernelsStatic::OptGeneric    ) { HOST_CALL(GenericDhopSiteDagExt); return;}
#ifndef GRID_NVCC
     if (Opt == WilsonKernelsStatic::OptHandUnroll ) { HOST_CALL(HandDhopSiteDagExt);    return;}
     if (Opt == WilsonKernelsStatic::OptInlineAsm  ) {  ASM_CALL(AsmDhopSiteDagExt);     return;}
#endif
   }
   assert(0 && " Kernel optimisation case not covered ");
  }

NAMESPACE_END(Grid);

