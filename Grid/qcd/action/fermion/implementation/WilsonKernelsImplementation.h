/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/WilsonKernels.cc

Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#pragma once

#include <Grid/qcd/action/fermion/FermionCore.h>

NAMESPACE_BEGIN(Grid);


////////////////////////////////////////////
// Generic implementation; move to different file?
////////////////////////////////////////////

accelerator_inline void get_stencil(StencilEntry * mem, StencilEntry &chip)
{
#ifdef __CUDA_ARCH__
  static_assert(sizeof(StencilEntry)==sizeof(uint4),"Unexpected Stencil Entry Size"); 
  uint4 * mem_pun  = (uint4 *)mem; // force 128 bit loads
  uint4 * chip_pun = (uint4 *)&chip;
  * chip_pun = * mem_pun;
#else 
  chip = *mem;
#endif
  return;
}
  
#define GENERIC_STENCIL_LEG(Dir,spProj,Recon)			\
  SE = st.GetEntry(ptype, Dir, sF);				\
  if (SE->_is_local) {						\
    int perm= SE->_permute;					\
    auto tmp = coalescedReadPermute(in[SE->_offset],ptype,perm,lane);	\
    spProj(chi,tmp);						\
  } else {							\
    chi = coalescedRead(buf[SE->_offset],lane);			\
  }								\
  synchronise();						\
  Impl::multLink(Uchi, U[sU], chi, Dir, SE, st);		\
  Recon(result, Uchi);
  
#define GENERIC_STENCIL_LEG_INT(Dir,spProj,Recon)		\
  SE = st.GetEntry(ptype, Dir, sF);				\
  if (SE->_is_local) {						\
    int perm= SE->_permute;					\
    auto tmp = coalescedReadPermute(in[SE->_offset],ptype,perm,lane);	\
    spProj(chi,tmp);						\
  } else if ( st.same_node[Dir] ) {				\
    chi = coalescedRead(buf[SE->_offset],lane);			\
  }								\
  synchronise();						\
  if (SE->_is_local || st.same_node[Dir] ) {			\
    Impl::multLink(Uchi, U[sU], chi, Dir, SE, st);		\
    Recon(result, Uchi);					\
  }								\
  synchronise();						

#define GENERIC_STENCIL_LEG_EXT(Dir,spProj,Recon)		\
  SE = st.GetEntry(ptype, Dir, sF);				\
  if ((!SE->_is_local) && (!st.same_node[Dir]) ) {		\
    auto chi = coalescedRead(buf[SE->_offset],lane);		\
    Impl::multLink(Uchi, U[sU], chi, Dir, SE, st);		\
    Recon(result, Uchi);					\
    nmu++;							\
  }								\
  synchronise();						

#define GENERIC_DHOPDIR_LEG_BODY(Dir,spProj,Recon)		\
    if (SE->_is_local ) {					\
      int perm= SE->_permute;					\
      auto tmp = coalescedReadPermute(in[SE->_offset],ptype,perm,lane);	\
      spProj(chi,tmp);						\
    } else {							\
      chi = coalescedRead(buf[SE->_offset],lane);		\
    }								\
    synchronise();						\
    Impl::multLink(Uchi, U[sU], chi, dir, SE, st);		\
    Recon(result, Uchi);					

#define GENERIC_DHOPDIR_LEG(Dir,spProj,Recon)			\
  if (gamma == Dir) {						\
    GENERIC_DHOPDIR_LEG_BODY(Dir,spProj,Recon);			\
  }


  ////////////////////////////////////////////////////////////////////
  // All legs kernels ; comms then compute
  ////////////////////////////////////////////////////////////////////
template <class Impl>
void WilsonKernels<Impl>::GenericDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U,
					     SiteHalfSpinor *buf, int sF,
					     int sU, const FermionFieldView &in, FermionFieldView &out)
{
  typedef decltype(coalescedRead(buf[0]))   calcHalfSpinor;
  typedef decltype(coalescedRead(in[0])) calcSpinor;
  calcHalfSpinor chi;
  //  calcHalfSpinor *chi_p;
  calcHalfSpinor Uchi;
  calcSpinor result;
  StencilEntry *SE;
  int ptype;
  const int Nsimd = SiteHalfSpinor::Nsimd();
  const int lane=SIMTlane(Nsimd);
  GENERIC_STENCIL_LEG(Xp,spProjXp,spReconXp);
  GENERIC_STENCIL_LEG(Yp,spProjYp,accumReconYp);
  GENERIC_STENCIL_LEG(Zp,spProjZp,accumReconZp);
  GENERIC_STENCIL_LEG(Tp,spProjTp,accumReconTp);
  GENERIC_STENCIL_LEG(Xm,spProjXm,accumReconXm);
  GENERIC_STENCIL_LEG(Ym,spProjYm,accumReconYm);
  GENERIC_STENCIL_LEG(Zm,spProjZm,accumReconZm);
  GENERIC_STENCIL_LEG(Tm,spProjTm,accumReconTm);
  coalescedWrite(out[sF],result,lane);
};

template <class Impl>
void WilsonKernels<Impl>::GenericDhopSite(StencilView &st, DoubledGaugeFieldView &U,
					  SiteHalfSpinor *buf, int sF,
					  int sU, const FermionFieldView &in, FermionFieldView &out) 
{
  typedef decltype(coalescedRead(buf[0])) calcHalfSpinor;
  typedef decltype(coalescedRead(in[0]))  calcSpinor;
  calcHalfSpinor chi;
  //  calcHalfSpinor *chi_p;
  calcHalfSpinor Uchi;
  calcSpinor result;
  StencilEntry *SE;
  int ptype;

  const int Nsimd = SiteHalfSpinor::Nsimd();
  const int lane=SIMTlane(Nsimd);
  GENERIC_STENCIL_LEG(Xm,spProjXp,spReconXp);
  GENERIC_STENCIL_LEG(Ym,spProjYp,accumReconYp);
  GENERIC_STENCIL_LEG(Zm,spProjZp,accumReconZp);
  GENERIC_STENCIL_LEG(Tm,spProjTp,accumReconTp);
  GENERIC_STENCIL_LEG(Xp,spProjXm,accumReconXm);
  GENERIC_STENCIL_LEG(Yp,spProjYm,accumReconYm);
  GENERIC_STENCIL_LEG(Zp,spProjZm,accumReconZm);
  GENERIC_STENCIL_LEG(Tp,spProjTm,accumReconTm);
  coalescedWrite(out[sF], result,lane);
};
  ////////////////////////////////////////////////////////////////////
  // Interior kernels
  ////////////////////////////////////////////////////////////////////
template <class Impl>
void WilsonKernels<Impl>::GenericDhopSiteDagInt(StencilView &st,  DoubledGaugeFieldView &U,
						SiteHalfSpinor *buf, int sF,
						int sU, const FermionFieldView &in, FermionFieldView &out)
{
  typedef decltype(coalescedRead(buf[0])) calcHalfSpinor;
  typedef decltype(coalescedRead(in[0]))  calcSpinor;
  calcHalfSpinor chi;
  //  calcHalfSpinor *chi_p;
  calcHalfSpinor Uchi;
  calcSpinor result;
  StencilEntry *SE;
  int ptype;
  const int Nsimd = SiteHalfSpinor::Nsimd();
  const int lane=SIMTlane(Nsimd);

  result=Zero();
  GENERIC_STENCIL_LEG_INT(Xp,spProjXp,accumReconXp);
  GENERIC_STENCIL_LEG_INT(Yp,spProjYp,accumReconYp);
  GENERIC_STENCIL_LEG_INT(Zp,spProjZp,accumReconZp);
  GENERIC_STENCIL_LEG_INT(Tp,spProjTp,accumReconTp);
  GENERIC_STENCIL_LEG_INT(Xm,spProjXm,accumReconXm);
  GENERIC_STENCIL_LEG_INT(Ym,spProjYm,accumReconYm);
  GENERIC_STENCIL_LEG_INT(Zm,spProjZm,accumReconZm);
  GENERIC_STENCIL_LEG_INT(Tm,spProjTm,accumReconTm);
  coalescedWrite(out[sF], result,lane);
};

template <class Impl>
void WilsonKernels<Impl>::GenericDhopSiteInt(StencilView &st,  DoubledGaugeFieldView &U,
							 SiteHalfSpinor *buf, int sF,
							 int sU, const FermionFieldView &in, FermionFieldView &out) 
{
  typedef decltype(coalescedRead(buf[0])) calcHalfSpinor;
  typedef decltype(coalescedRead(in[0]))  calcSpinor;
  const int Nsimd = SiteHalfSpinor::Nsimd();
  const int lane=SIMTlane(Nsimd);

  calcHalfSpinor chi;
  //  calcHalfSpinor *chi_p;
  calcHalfSpinor Uchi;
  calcSpinor result;
  StencilEntry *SE;
  int ptype;
  result=Zero();
  GENERIC_STENCIL_LEG_INT(Xm,spProjXp,accumReconXp);
  GENERIC_STENCIL_LEG_INT(Ym,spProjYp,accumReconYp);
  GENERIC_STENCIL_LEG_INT(Zm,spProjZp,accumReconZp);
  GENERIC_STENCIL_LEG_INT(Tm,spProjTp,accumReconTp);
  GENERIC_STENCIL_LEG_INT(Xp,spProjXm,accumReconXm);
  GENERIC_STENCIL_LEG_INT(Yp,spProjYm,accumReconYm);
  GENERIC_STENCIL_LEG_INT(Zp,spProjZm,accumReconZm);
  GENERIC_STENCIL_LEG_INT(Tp,spProjTm,accumReconTm);
  coalescedWrite(out[sF], result,lane);
};
////////////////////////////////////////////////////////////////////
// Exterior kernels
////////////////////////////////////////////////////////////////////
template <class Impl>
void WilsonKernels<Impl>::GenericDhopSiteDagExt(StencilView &st,  DoubledGaugeFieldView &U,
						SiteHalfSpinor *buf, int sF,
						int sU, const FermionFieldView &in, FermionFieldView &out)
{
  typedef decltype(coalescedRead(buf[0])) calcHalfSpinor;
  typedef decltype(coalescedRead(in[0]))  calcSpinor;
  //  calcHalfSpinor *chi_p;
  calcHalfSpinor Uchi;
  calcSpinor result;
  StencilEntry *SE;
  int ptype;
  int nmu=0;
  const int Nsimd = SiteHalfSpinor::Nsimd();
  const int lane=SIMTlane(Nsimd);
  result=Zero();
  GENERIC_STENCIL_LEG_EXT(Xp,spProjXp,accumReconXp);
  GENERIC_STENCIL_LEG_EXT(Yp,spProjYp,accumReconYp);
  GENERIC_STENCIL_LEG_EXT(Zp,spProjZp,accumReconZp);
  GENERIC_STENCIL_LEG_EXT(Tp,spProjTp,accumReconTp);
  GENERIC_STENCIL_LEG_EXT(Xm,spProjXm,accumReconXm);
  GENERIC_STENCIL_LEG_EXT(Ym,spProjYm,accumReconYm);
  GENERIC_STENCIL_LEG_EXT(Zm,spProjZm,accumReconZm);
  GENERIC_STENCIL_LEG_EXT(Tm,spProjTm,accumReconTm);
  if ( nmu ) { 
    auto out_t = coalescedRead(out[sF],lane);
    out_t = out_t + result;
    coalescedWrite(out[sF],out_t,lane);
  }
};

template <class Impl>
void WilsonKernels<Impl>::GenericDhopSiteExt(StencilView &st,  DoubledGaugeFieldView &U,
					     SiteHalfSpinor *buf, int sF,
					     int sU, const FermionFieldView &in, FermionFieldView &out) 
{
  typedef decltype(coalescedRead(buf[0])) calcHalfSpinor;
  typedef decltype(coalescedRead(in[0]))  calcSpinor;
  //  calcHalfSpinor *chi_p;
  calcHalfSpinor Uchi;
  calcSpinor result;
  StencilEntry *SE;
  int ptype;
  int nmu=0;
  const int Nsimd = SiteHalfSpinor::Nsimd();
  const int lane=SIMTlane(Nsimd);
  result=Zero();
  GENERIC_STENCIL_LEG_EXT(Xm,spProjXp,accumReconXp);
  GENERIC_STENCIL_LEG_EXT(Ym,spProjYp,accumReconYp);
  GENERIC_STENCIL_LEG_EXT(Zm,spProjZp,accumReconZp);
  GENERIC_STENCIL_LEG_EXT(Tm,spProjTp,accumReconTp);
  GENERIC_STENCIL_LEG_EXT(Xp,spProjXm,accumReconXm);
  GENERIC_STENCIL_LEG_EXT(Yp,spProjYm,accumReconYm);
  GENERIC_STENCIL_LEG_EXT(Zp,spProjZm,accumReconZm);
  GENERIC_STENCIL_LEG_EXT(Tp,spProjTm,accumReconTm);
  if ( nmu ) { 
    auto out_t = coalescedRead(out[sF],lane);
    out_t = out_t + result;
    coalescedWrite(out[sF],out_t,lane);
  }
};

#define DhopDirMacro(Dir,spProj,spRecon)	\
  template <class Impl>							\
  void WilsonKernels<Impl>::DhopDir##Dir(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf, int sF, \
					 int sU, const FermionFieldView &in, FermionFieldView &out, int dir) \
  {									\
  typedef decltype(coalescedRead(buf[0])) calcHalfSpinor;		\
  typedef decltype(coalescedRead(in[0]))  calcSpinor;			\
  calcHalfSpinor chi;							\
  calcSpinor result;							\
  calcHalfSpinor Uchi;							\
  StencilEntry *SE;							\
  int ptype;								\
  const int Nsimd = SiteHalfSpinor::Nsimd();				\
  const int lane=SIMTlane(Nsimd);					\
									\
  SE = st.GetEntry(ptype, dir, sF);					\
  GENERIC_DHOPDIR_LEG_BODY(Dir,spProj,spRecon);				\
  coalescedWrite(out[sF], result,lane);					\
  }									

DhopDirMacro(Xp,spProjXp,spReconXp);
DhopDirMacro(Yp,spProjYp,spReconYp);
DhopDirMacro(Zp,spProjZp,spReconZp);
DhopDirMacro(Tp,spProjTp,spReconTp);
DhopDirMacro(Xm,spProjXm,spReconXm);
DhopDirMacro(Ym,spProjYm,spReconYm);
DhopDirMacro(Zm,spProjZm,spReconZm);
DhopDirMacro(Tm,spProjTm,spReconTm);

template <class Impl> 
void WilsonKernels<Impl>::DhopDirK( StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf, int sF,
				    int sU, const FermionFieldView &in, FermionFieldView &out, int dir, int gamma) 
{
  typedef decltype(coalescedRead(buf[0])) calcHalfSpinor;
  typedef decltype(coalescedRead(in[0]))  calcSpinor;
  calcHalfSpinor chi;
  calcSpinor result;
  calcHalfSpinor Uchi;
  StencilEntry *SE;
  int ptype;
  const int Nsimd = SiteHalfSpinor::Nsimd();
  const int lane=SIMTlane(Nsimd);

  SE = st.GetEntry(ptype, dir, sF);
  GENERIC_DHOPDIR_LEG(Xp,spProjXp,spReconXp);
  GENERIC_DHOPDIR_LEG(Yp,spProjYp,spReconYp);
  GENERIC_DHOPDIR_LEG(Zp,spProjZp,spReconZp);
  GENERIC_DHOPDIR_LEG(Tp,spProjTp,spReconTp);
  GENERIC_DHOPDIR_LEG(Xm,spProjXm,spReconXm);
  GENERIC_DHOPDIR_LEG(Ym,spProjYm,spReconYm);
  GENERIC_DHOPDIR_LEG(Zm,spProjZm,spReconZm);
  GENERIC_DHOPDIR_LEG(Tm,spProjTm,spReconTm);
  coalescedWrite(out[sF], result,lane);
}

template <class Impl>
void WilsonKernels<Impl>::DhopDirAll( StencilImpl &st, DoubledGaugeField &U,SiteHalfSpinor *buf, int Ls,
				      int Nsite, const FermionField &in, std::vector<FermionField> &out) 
{
   auto U_v   = U.View();
   auto in_v  = in.View();
   auto st_v  = st.View();

   auto out_Xm = out[0].View();
   auto out_Ym = out[1].View();
   auto out_Zm = out[2].View();
   auto out_Tm = out[3].View();
   auto out_Xp = out[4].View();
   auto out_Yp = out[5].View();
   auto out_Zp = out[6].View();
   auto out_Tp = out[7].View();

   accelerator_forNB(sss,Nsite*Ls,Simd::Nsimd(),{
      int sU=sss/Ls;				
      int sF =sss;				
      DhopDirXm(st_v,U_v,st.CommBuf(),sF,sU,in_v,out_Xm,0);
      DhopDirYm(st_v,U_v,st.CommBuf(),sF,sU,in_v,out_Ym,1);
      DhopDirZm(st_v,U_v,st.CommBuf(),sF,sU,in_v,out_Zm,2);
      DhopDirTm(st_v,U_v,st.CommBuf(),sF,sU,in_v,out_Tm,3);
      DhopDirXp(st_v,U_v,st.CommBuf(),sF,sU,in_v,out_Xp,4);
      DhopDirYp(st_v,U_v,st.CommBuf(),sF,sU,in_v,out_Yp,5);
      DhopDirZp(st_v,U_v,st.CommBuf(),sF,sU,in_v,out_Zp,6);
      DhopDirTp(st_v,U_v,st.CommBuf(),sF,sU,in_v,out_Tp,7);
   });
}


template <class Impl>
void WilsonKernels<Impl>::DhopDirKernel( StencilImpl &st, DoubledGaugeField &U,SiteHalfSpinor *buf, int Ls,
					 int Nsite, const FermionField &in, FermionField &out, int dirdisp, int gamma) 
{
  assert(dirdisp<=7);
  assert(dirdisp>=0);

   auto U_v   = U.View();
   auto in_v  = in.View();
   auto out_v = out.View();
   auto st_v  = st.View();
#define LoopBody(Dir)				\
   case Dir :			\
     accelerator_forNB(ss,Nsite,Simd::Nsimd(),{	\
       for(int s=0;s<Ls;s++){			\
	 int sU=ss;				\
	 int sF = s+Ls*sU;						\
	 DhopDir##Dir(st_v,U_v,st.CommBuf(),sF,sU,in_v,out_v,dirdisp);\
       }							       \
       });							       \
     break;

   switch(gamma){
   LoopBody(Xp);
   LoopBody(Yp);
   LoopBody(Zp);
   LoopBody(Tp);

   LoopBody(Xm);
   LoopBody(Ym);
   LoopBody(Zm);
   LoopBody(Tm);
   default:
     assert(0);
     break;
   }
#undef LoopBody
} 

#define KERNEL_CALLNB(A) \
  const uint64_t    NN = Nsite*Ls;					\
  accelerator_forNB( ss, NN, Simd::Nsimd(), {				\
      int sF = ss;							\
      int sU = ss/Ls;							\
      WilsonKernels<Impl>::A(st_v,U_v,buf,sF,sU,in_v,out_v);		\
  });

#define KERNEL_CALL(A) KERNEL_CALLNB(A); accelerator_barrier(); 

#define ASM_CALL(A)							\
  thread_for( ss, Nsite, {						\
    int sU = ss;							\
    int sF = ss*Ls;							\
    WilsonKernels<Impl>::A(st_v,U_v,buf,sF,sU,Ls,1,in_v,out_v);		\
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
     if (Opt == WilsonKernelsStatic::OptGeneric    ) { KERNEL_CALL(GenericDhopSite); return;}
#ifndef GRID_NVCC
     if (Opt == WilsonKernelsStatic::OptHandUnroll ) { KERNEL_CALL(HandDhopSite);    return;}
     if (Opt == WilsonKernelsStatic::OptInlineAsm  ) {  ASM_CALL(AsmDhopSite);    return;}
#endif
   } else if( interior ) {
     if (Opt == WilsonKernelsStatic::OptGeneric    ) { KERNEL_CALLNB(GenericDhopSiteInt); return;}
#ifndef GRID_NVCC
     if (Opt == WilsonKernelsStatic::OptHandUnroll ) { KERNEL_CALLNB(HandDhopSiteInt);    return;}
     if (Opt == WilsonKernelsStatic::OptInlineAsm  ) {  ASM_CALL(AsmDhopSiteInt);    return;}
#endif
   } else if( exterior ) { 
     if (Opt == WilsonKernelsStatic::OptGeneric    ) { KERNEL_CALL(GenericDhopSiteExt); return;}
#ifndef GRID_NVCC
     if (Opt == WilsonKernelsStatic::OptHandUnroll ) { KERNEL_CALL(HandDhopSiteExt);    return;}
     if (Opt == WilsonKernelsStatic::OptInlineAsm  ) {  ASM_CALL(AsmDhopSiteExt);    return;}
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
     if (Opt == WilsonKernelsStatic::OptGeneric    ) { KERNEL_CALL(GenericDhopSiteDag); return;}
#ifndef GRID_NVCC
     if (Opt == WilsonKernelsStatic::OptHandUnroll ) { KERNEL_CALL(HandDhopSiteDag);    return;}
     if (Opt == WilsonKernelsStatic::OptInlineAsm  ) {  ASM_CALL(AsmDhopSiteDag);     return;}
#endif
   } else if( interior ) {
     if (Opt == WilsonKernelsStatic::OptGeneric    ) { KERNEL_CALL(GenericDhopSiteDagInt); return;}
#ifndef GRID_NVCC
     if (Opt == WilsonKernelsStatic::OptHandUnroll ) { KERNEL_CALL(HandDhopSiteDagInt);    return;}
     if (Opt == WilsonKernelsStatic::OptInlineAsm  ) {  ASM_CALL(AsmDhopSiteDagInt);     return;}
#endif
   } else if( exterior ) { 
     if (Opt == WilsonKernelsStatic::OptGeneric    ) { KERNEL_CALL(GenericDhopSiteDagExt); return;}
#ifndef GRID_NVCC
     if (Opt == WilsonKernelsStatic::OptHandUnroll ) { KERNEL_CALL(HandDhopSiteDagExt);    return;}
     if (Opt == WilsonKernelsStatic::OptInlineAsm  ) {  ASM_CALL(AsmDhopSiteDagExt);     return;}
#endif
   }
   assert(0 && " Kernel optimisation case not covered ");
  }

NAMESPACE_END(Grid);

