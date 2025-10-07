/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/TwoSpinWilsonKernels.cc

Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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

#define GENERIC_STENCIL_LEG(Dir,spProj,Recon)			\
  SE = st.GetEntry(ptype, Dir, sF);				\
  if (SE->_is_local) {						\
    int perm= SE->_permute;					\
    auto tmp = coalescedReadPermute(in[SE->_offset],ptype,perm,lane);	\
    spProj(chi,tmp);						\
  } else {							\
    chi = coalescedRead(buf[SE->_offset],lane);			\
  }								\
  acceleratorSynchronise();					\
  Impl::multLink(Uchi, U[sU], chi, Dir, SE, st);		\
  Recon(result, Uchi);

#define GENERIC_STENCIL_LEG_INT(Dir,spProj,Recon)		\
  SE = st.GetEntry(ptype, Dir, sF);				\
  if (SE->_is_local) {						\
    int perm= SE->_permute;					\
    auto tmp = coalescedReadPermute(in[SE->_offset],ptype,perm,lane);	\
    spProj(chi,tmp);							\
    Impl::multLink(Uchi, U[sU], chi, Dir, SE, st);			\
    Recon(result, Uchi);						\
  }									\
  acceleratorSynchronise();

#define GENERIC_STENCIL_LEG_EXT(Dir,spProj,Recon)		\
  SE = st.GetEntry(ptype, Dir, sF);				\
  if (!SE->_is_local ) {		\
    auto chi = coalescedRead(buf[SE->_offset],lane);		\
    Impl::multLink(Uchi, U[sU], chi, Dir, SE, st);		\
    Recon(result, Uchi);					\
    nmu++;							\
  }								\
  acceleratorSynchronise();

#define GENERIC_DHOPDIR_LEG_BODY(Dir,spProj,Recon)		\
    if (SE->_is_local ) {					\
      int perm= SE->_permute;					\
      auto tmp = coalescedReadPermute(in[SE->_offset],ptype,perm,lane);	\
      spProj(chi,tmp);						\
    } else {							\
      chi = coalescedRead(buf[SE->_offset],lane);		\
    }								\
    acceleratorSynchronise();					\
    Impl::multLink(Uchi, U[sU], chi, dir, SE, st);		\
    Recon(result, Uchi);

#define GENERIC_DHOPDIR_LEG(Dir,spProj,Recon)			\
  if (gamma == Dir) {						\
    GENERIC_DHOPDIR_LEG_BODY(Dir,spProj,Recon);			\
  }

////////////////////////////////////////////////////////////////////
// All legs kernels ; comms then compute
////////////////////////////////////////////////////////////////////
template <class Impl> accelerator_inline
void TwoSpinWilsonKernels<Impl>::DhopSiteDag(StencilView &st, DoubledGaugeFieldView &U,
					     SiteSpinor *buf, int sF,
					     int sU, const FermionFieldView &in, FermionFieldView &out)
{
  typedef decltype(coalescedRead(in[0])) calcSpinor;
  calcSpinor chi;
  calcSpinor Uchi;
  calcSpinor result;
  StencilEntry *SE;
  int ptype;
  const int Nsimd = SiteSpinor::Nsimd();
  const int lane=acceleratorSIMTlane(Nsimd);
  GENERIC_STENCIL_LEG(Xp,pauliProjXp,pauliAssign);
  GENERIC_STENCIL_LEG(Yp,pauliProjYp,pauliAdd);
  GENERIC_STENCIL_LEG(Zp,pauliProjZp,pauliAdd);
  GENERIC_STENCIL_LEG(Xm,pauliProjXm,pauliAdd);
  GENERIC_STENCIL_LEG(Ym,pauliProjYm,pauliAdd);
  GENERIC_STENCIL_LEG(Zm,pauliProjZm,pauliAdd);
  coalescedWrite(out[sF],result,lane);
};

template <class Impl> accelerator_inline
void TwoSpinWilsonKernels<Impl>::GenericDhopSite(StencilView &st, DoubledGaugeFieldView &U,
					  SiteSpinor *buf, int sF,
					  int sU, const FermionFieldView &in, FermionFieldView &out)
{
  typedef decltype(coalescedRead(in[0]))  calcSpinor;
  calcSpinor chi;
  //  calcSpinor *chi_p;
  calcSpinor Uchi;
  calcSpinor result;
  StencilEntry *SE;
  int ptype;

  const int Nsimd = SiteSpinor::Nsimd();
  const int lane=acceleratorSIMTlane(Nsimd);
  GENERIC_STENCIL_LEG(Xm,pauliProjXp,pauliAssign);
  GENERIC_STENCIL_LEG(Ym,pauliProjYp,pauliAdd);
  GENERIC_STENCIL_LEG(Zm,pauliProjZp,pauliAdd);
  GENERIC_STENCIL_LEG(Xp,pauliProjXm,pauliAdd);
  GENERIC_STENCIL_LEG(Yp,pauliProjYm,pauliAdd);
  GENERIC_STENCIL_LEG(Zp,pauliProjZm,pauliAdd);
  coalescedWrite(out[sF], result,lane);
};
  ////////////////////////////////////////////////////////////////////
  // Interior kernels
  ////////////////////////////////////////////////////////////////////
template <class Impl> accelerator_inline
void TwoSpinWilsonKernels<Impl>::GenericDhopSiteDagInt(StencilView &st,  DoubledGaugeFieldView &U,
						       SiteSpinor *buf, int sF,
						       int sU, const FermionFieldView &in, FermionFieldView &out)
{
  typedef decltype(coalescedRead(in[0]))  calcSpinor;
  calcSpinor chi;
  //  calcSpinor *chi_p;
  calcSpinor Uchi;
  calcSpinor result;
  StencilEntry *SE;
  int ptype;
  const int Nsimd = SiteSpinor::Nsimd();
  const int lane=acceleratorSIMTlane(Nsimd);

  result=Zero();
  GENERIC_STENCIL_LEG_INT(Xp,pauliProjXp,pauliAdd);
  GENERIC_STENCIL_LEG_INT(Yp,pauliProjYp,pauliAdd);
  GENERIC_STENCIL_LEG_INT(Zp,pauliProjZp,pauliAdd);
  GENERIC_STENCIL_LEG_INT(Xm,pauliProjXm,pauliAdd);
  GENERIC_STENCIL_LEG_INT(Ym,pauliProjYm,pauliAdd);
  GENERIC_STENCIL_LEG_INT(Zm,pauliProjZm,pauliAdd);
  coalescedWrite(out[sF], result,lane);
};

template <class Impl> accelerator_inline
void TwoSpinWilsonKernels<Impl>::GenericDhopSiteInt(StencilView &st,  DoubledGaugeFieldView &U,
						    SiteSpinor *buf, int sF,
						    int sU, const FermionFieldView &in, FermionFieldView &out)
{
  typedef decltype(coalescedRead(in[0]))  calcSpinor;
  const int Nsimd = SiteSpinor::Nsimd();
  const int lane=acceleratorSIMTlane(Nsimd);

  calcSpinor chi;
  //  calcSpinor *chi_p;
  calcSpinor Uchi;
  calcSpinor result;
  StencilEntry *SE;
  int ptype;
  result=Zero();
  GENERIC_STENCIL_LEG_INT(Xm,pauliProjXp,pauliAdd);
  GENERIC_STENCIL_LEG_INT(Ym,pauliProjYp,pauliAdd);
  GENERIC_STENCIL_LEG_INT(Zm,pauliProjZp,pauliAdd);
  GENERIC_STENCIL_LEG_INT(Xp,pauliProjXm,pauliAdd);
  GENERIC_STENCIL_LEG_INT(Yp,pauliProjYm,pauliAdd);
  GENERIC_STENCIL_LEG_INT(Zp,pauliProjZm,pauliAdd);
  coalescedWrite(out[sF], result,lane);
};
////////////////////////////////////////////////////////////////////
// Exterior kernels
////////////////////////////////////////////////////////////////////
template <class Impl> accelerator_inline
void TwoSpinWilsonKernels<Impl>::GenericDhopSiteDagExt(StencilView &st,  DoubledGaugeFieldView &U,
						SiteSpinor *buf, int sF,
						int sU, const FermionFieldView &in, FermionFieldView &out)
{
  typedef decltype(coalescedRead(in[0]))  calcSpinor;
  //  calcSpinor *chi_p;
  calcSpinor Uchi;
  calcSpinor result;
  StencilEntry *SE;
  int ptype;
  int nmu=0;
  const int Nsimd = SiteSpinor::Nsimd();
  const int lane=acceleratorSIMTlane(Nsimd);
  result=Zero();
  GENERIC_STENCIL_LEG_EXT(Xp,pauliProjXp,pauliAdd);
  GENERIC_STENCIL_LEG_EXT(Yp,pauliProjYp,pauliAdd);
  GENERIC_STENCIL_LEG_EXT(Zp,pauliProjZp,pauliAdd);
  GENERIC_STENCIL_LEG_EXT(Xm,pauliProjXm,pauliAdd);
  GENERIC_STENCIL_LEG_EXT(Ym,pauliProjYm,pauliAdd);
  GENERIC_STENCIL_LEG_EXT(Zm,pauliProjZm,pauliAdd);
  if ( nmu ) {
    auto out_t = coalescedRead(out[sF],lane);
    out_t = out_t + result;
    coalescedWrite(out[sF],out_t,lane);
  }
};

template <class Impl> accelerator_inline
void TwoSpinWilsonKernels<Impl>::GenericDhopSiteExt(StencilView &st,  DoubledGaugeFieldView &U,
					     SiteSpinor *buf, int sF,
					     int sU, const FermionFieldView &in, FermionFieldView &out)
{
  typedef decltype(coalescedRead(in[0]))  calcSpinor;
  //  calcSpinor *chi_p;
  calcSpinor Uchi;
  calcSpinor result;
  StencilEntry *SE;
  int ptype;
  int nmu=0;
  const int Nsimd = SiteSpinor::Nsimd();
  const int lane=acceleratorSIMTlane(Nsimd);
  result=Zero();
  GENERIC_STENCIL_LEG_EXT(Xm,pauliProjXp,pauliAdd);
  GENERIC_STENCIL_LEG_EXT(Ym,pauliProjYp,pauliAdd);
  GENERIC_STENCIL_LEG_EXT(Zm,pauliProjZp,pauliAdd);
  GENERIC_STENCIL_LEG_EXT(Xp,pauliProjXm,pauliAdd);
  GENERIC_STENCIL_LEG_EXT(Yp,pauliProjYm,pauliAdd);
  GENERIC_STENCIL_LEG_EXT(Zp,pauliProjZm,pauliAdd);
  if ( nmu ) {
    auto out_t = coalescedRead(out[sF],lane);
    out_t = out_t + result;
    coalescedWrite(out[sF],out_t,lane);
  }
};

#define DhopDirMacro(Dir,spProj,spRecon)	\
  template <class Impl> accelerator_inline				\
  void TwoSpinWilsonKernels<Impl>::DhopDir##Dir(StencilView &st, DoubledGaugeFieldView &U,SiteSpinor *buf, int sF, \
					 int sU, const FermionFieldView &in, FermionFieldView &out, int dir) \
  {									\
  typedef decltype(coalescedRead(in[0]))  calcSpinor;			\
  calcSpinor chi;							\
  calcSpinor result;							\
  calcSpinor Uchi;							\
  StencilEntry *SE;							\
  int ptype;								\
  const int Nsimd = SiteSpinor::Nsimd();				\
  const int lane=acceleratorSIMTlane(Nsimd);					\
									\
  SE = st.GetEntry(ptype, dir, sF);					\
  GENERIC_DHOPDIR_LEG_BODY(Dir,spProj,spRecon);				\
  coalescedWrite(out[sF], result,lane);					\
  }

DhopDirMacro(Xp,pauliProjXp,pauliAssign);
DhopDirMacro(Yp,pauliProjYp,pauliAssign);
DhopDirMacro(Zp,pauliProjZp,pauliAssign);
DhopDirMacro(Xm,pauliProjXm,pauliAssign);
DhopDirMacro(Ym,pauliProjYm,pauliAssign);
DhopDirMacro(Zm,pauliProjZm,pauliAssign);

template <class Impl> accelerator_inline
void TwoSpinWilsonKernels<Impl>::DhopDirK( StencilView &st, DoubledGaugeFieldView &U,SiteSpinor *buf, int sF,
				    int sU, const FermionFieldView &in, FermionFieldView &out, int dir, int gamma)
{
  typedef decltype(coalescedRead(in[0]))  calcSpinor;
  calcSpinor chi;
  calcSpinor result;
  calcSpinor Uchi;
  StencilEntry *SE;
  int ptype;
  const int Nsimd = SiteSpinor::Nsimd();
  const int lane=acceleratorSIMTlane(Nsimd);

  SE = st.GetEntry(ptype, dir, sF);
  GENERIC_DHOPDIR_LEG(Xp,pauliProjXp,pauliAssign);
  GENERIC_DHOPDIR_LEG(Yp,pauliProjYp,pauliAssign);
  GENERIC_DHOPDIR_LEG(Zp,pauliProjZp,pauliAssign);
  GENERIC_DHOPDIR_LEG(Xm,pauliProjXm,pauliAssign);
  GENERIC_DHOPDIR_LEG(Ym,pauliProjYm,pauliAssign);
  GENERIC_DHOPDIR_LEG(Zm,pauliProjZm,pauliAssign);
  coalescedWrite(out[sF], result,lane);
}

template <class Impl>
void TwoSpinWilsonKernels<Impl>::DhopDirAll( StencilImpl &st, DoubledGaugeField &U,SiteSpinor *buf, int Ls,
				      int Nsite, const FermionField &in, std::vector<FermionField> &out)
{
   autoView(U_v  ,U,AcceleratorRead);
   autoView(in_v ,in,AcceleratorRead);
   autoView(st_v ,st,AcceleratorRead);

   autoView(out_Xm,out[0],AcceleratorWrite);
   autoView(out_Ym,out[1],AcceleratorWrite);
   autoView(out_Zm,out[2],AcceleratorWrite);
   autoView(out_Xp,out[4],AcceleratorWrite);
   autoView(out_Yp,out[5],AcceleratorWrite);
   autoView(out_Zp,out[6],AcceleratorWrite);
   auto CBp=st.CommBuf();
   accelerator_for(sss,Nsite*Ls,Simd::Nsimd(),{
      int sU=sss/Ls;
      int sF =sss;
      DhopDirXm(st_v,U_v,CBp,sF,sU,in_v,out_Xm,0);
      DhopDirYm(st_v,U_v,CBp,sF,sU,in_v,out_Ym,1);
      DhopDirZm(st_v,U_v,CBp,sF,sU,in_v,out_Zm,2);
      DhopDirXp(st_v,U_v,CBp,sF,sU,in_v,out_Xp,3);
      DhopDirYp(st_v,U_v,CBp,sF,sU,in_v,out_Yp,4);
      DhopDirZp(st_v,U_v,CBp,sF,sU,in_v,out_Zp,5);
   });
}


template <class Impl>
void TwoSpinWilsonKernels<Impl>::DhopDirKernel( StencilImpl &st, DoubledGaugeField &U,SiteSpinor *buf, int Ls,
					 int Nsite, const FermionField &in, FermionField &out, int dirdisp, int gamma)
{
  assert(dirdisp<=5);
  assert(dirdisp>=0);

   autoView(U_v  ,U  ,AcceleratorRead);
   autoView(in_v ,in ,AcceleratorRead);
   autoView(out_v,out,AcceleratorWrite);
   autoView(st_v ,st ,AcceleratorRead);
   auto CBp=st.CommBuf();
#define LoopBody(Dir)				\
   case Dir :					\
     accelerator_for(ss,Nsite,Simd::Nsimd(),{	\
       for(int s=0;s<Ls;s++){			\
	 int sU=ss;				\
	 int sF = s+Ls*sU;						\
	 DhopDir##Dir(st_v,U_v,CBp,sF,sU,in_v,out_v,dirdisp);\
       }							       \
       });							       \
     break;

   switch(gamma){
   LoopBody(Xp);
   LoopBody(Yp);
   LoopBody(Zp);

   LoopBody(Xm);
   LoopBody(Ym);
   LoopBody(Zm);
   default:
     assert(0);
     break;
   }
#undef LoopBody
}


#define KERNEL_CALLNB(A)						\
  const uint64_t    NN = Nsite*Ls;					\
  accelerator_forNB( ss, NN, Simd::Nsimd(), {				\
      int sF = ss;							\
      int sU = ss/Ls;							\
      TwoSpinWilsonKernels<Impl>::A(st_v,U_v,buf,sF,sU,in_v,out_v);		\
    });

#define KERNEL_CALL(A) KERNEL_CALLNB(A); accelerator_barrier();

#define KERNEL_CALL_EXT(A)						\
  const uint64_t    sz = st.surface_list.size();			\
  auto ptr = &st.surface_list[0];					\
  accelerator_forNB( ss, sz, Simd::Nsimd(), {				\
      int sF = ptr[ss];							\
      int sU = sF/Ls;							\
      TwoSpinWilsonKernels<Impl>::A(st_v,U_v,buf,sF,sU,in_v,out_v);		\
    });									\
  accelerator_barrier();


template <class Impl>
void TwoSpinWilsonKernels<Impl>::DhopKernel(StencilImpl &st,  DoubledGaugeField &U, SiteSpinor * buf,
					    int Ls, int Nsite, const FermionField &in, FermionField &out,
					    int interior,int exterior)
{
  autoView(U_v  ,  U,AcceleratorRead);
  autoView(in_v , in,AcceleratorRead);
  autoView(out_v,out,AcceleratorWrite);
  autoView(st_v , st,AcceleratorRead);
  
  if( interior && exterior ) {
    acceleratorFenceComputeStream();
    KERNEL_CALL(GenericDhopSite);
    return;
  } else if( interior ) {
    KERNEL_CALLNB(GenericDhopSiteInt);
    return;
  } else if( exterior ) {
    //     // dependent on result of merge
    acceleratorFenceComputeStream();
    KERNEL_CALL_EXT(GenericDhopSiteExt);
    return;
  }
  assert(0 && " Kernel optimisation case not covered ");
}

template <class Impl>
void TwoSpinWilsonKernels<Impl>::DhopDagKernel(StencilImpl &st,  DoubledGaugeField &U, SiteSpinor * buf,
					       int Ls, int Nsite, const FermionField &in, FermionField &out,
					       int interior,int exterior)
{
  autoView(U_v  ,U,AcceleratorRead);
  autoView(in_v ,in,AcceleratorRead);
  autoView(out_v,out,AcceleratorWrite);
  autoView(st_v ,st,AcceleratorRead);
  
  if( interior && exterior ) {
    acceleratorFenceComputeStream();
    KERNEL_CALL(GenericDhopSiteDag);
    return;
  } else if( interior ) {
    KERNEL_CALLNB(GenericDhopSiteDagInt); return;
  } else if( exterior ) {
    // Dependent on result of merge
    acceleratorFenceComputeStream();
    KERNEL_CALL_EXT(GenericDhopSiteDagExt); return;
  }
  assert(0 && " Kernel optimisation case not covered ");
}

#undef KERNEL_CALLNB
#undef KERNEL_CALL

NAMESPACE_END(Grid);
