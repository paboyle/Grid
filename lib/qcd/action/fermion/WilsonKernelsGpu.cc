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
#include <Grid/qcd/action/fermion/FermionCore.h>

NAMESPACE_BEGIN(Grid);

//////////////////////////////////////////////////////////////
// Gpu implementation; thread loop is implicit ; move to header
//////////////////////////////////////////////////////////////
accelerator_inline void synchronise(void) 
{
#ifdef __CUDA_ARCH__
  __syncthreads();
#endif
  return;
}
accelerator_inline int get_my_lanes(int Nsimd) 
{
#ifdef __CUDA_ARCH__
  return 1;
#else 
  return Nsimd;
#endif
}
accelerator_inline int get_my_lane_offset(int Nsimd) 
{
#ifdef __CUDA_ARCH__
  return ( (threadIdx.x) % Nsimd);
#else
  return 0;
#endif
}

accelerator_inline void get_stencil(StencilEntry * mem, StencilEntry &chip)
{
#if 0
  chip = *mem;
#else 
  assert(sizeof(StencilEntry)==sizeof(uint4));
  uint4 * mem_pun  = (uint4 *)mem;
  uint4 * chip_pun = (uint4 *)&chip;
  * chip_pun = * mem_pun;
#endif
  return;
}

#ifdef GPU_VEC
#if 1
#define GPU_COALESCED_STENCIL_LEG_PROJ(Dir,spProj)			\
  synchronise();							\
  if (SE._is_local) {							\
    int mask = Nsimd >> (ptype + 1);					\
    int plane= SE._permute ? (lane ^ mask) : lane;			\
    auto in_l = extractLane(plane,in[SE._offset+s]);			\
    spProj(chi,in_l);							\
  } else {								\
    chi  = extractLane(lane,buf[SE._offset+s]);			\
  }									\
  synchronise();
#else 
#define GPU_COALESCED_STENCIL_LEG_PROJ(Dir,spProj)			\
  { int mask = Nsimd >> (ptype + 1);					\
  int plane= SE._permute ? (lane ^ mask) : lane;			\
  synchronise();							\
  auto in_l = extractLane(plane,in[SE._offset+s]);			\
  synchronise();							\
  spProj(chi,in_l); }							
#endif
#else 
#define GPU_COALESCED_STENCIL_LEG_PROJ(Dir,spProj)			\
  synchronise();							\
  if (SE._is_local) {							\
    auto in_t = in[SE._offset+s];					\
    if (SE._permute) {							\
      spProj(tmp, in_t);						\
      permute(chi, tmp, ptype);						\
    } else {								\
      spProj(chi, in_t);						\
    }									\
  } else {								\
    chi  = buf[SE._offset+s];						\
  }									\
  synchronise();
#endif

template <class Impl>
accelerator_inline void WilsonKernels<Impl>::GpuDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U,
						     SiteHalfSpinor *buf, int Ls, int s,
						     int sU, const FermionFieldView &in, FermionFieldView &out)
{
#ifdef GPU_VEC
  typename SiteHalfSpinor::scalar_object chi;
  typename SiteHalfSpinor::scalar_object Uchi;
  typename SiteSpinor::scalar_object   result;
#else 
  SiteHalfSpinor chi;
  SiteHalfSpinor Uchi;
  SiteHalfSpinor tmp;
  SiteSpinor   result;
#endif
  typedef typename SiteSpinor::scalar_type scalar_type;
  typedef typename SiteSpinor::vector_type vector_type;
  constexpr int Nsimd = sizeof(vector_type)/sizeof(scalar_type);

  uint64_t lane_offset= get_my_lane_offset(Nsimd);
  uint64_t lanes      = get_my_lanes(Nsimd);

  StencilEntry *SE_mem;
  StencilEntry SE; 

  int ptype;
  uint64_t ssF = Ls * sU;
  uint64_t sF  = ssF + s;
#ifndef __CUDA_ARCH__
  for(int lane = lane_offset;lane<lane_offset+lanes;lane++){
#else
  int lane = lane_offset; {
#endif
    SE_mem = st.GetEntry(ptype, Xp, ssF); get_stencil(SE_mem,SE);
    GPU_COALESCED_STENCIL_LEG_PROJ(Xp,spProjXp); 
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Xp);
    spReconXp(result, Uchi);

    SE_mem = st.GetEntry(ptype, Yp, ssF); get_stencil(SE_mem,SE);
    GPU_COALESCED_STENCIL_LEG_PROJ(Yp,spProjYp);
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Yp);
    accumReconYp(result, Uchi);
      
    SE_mem = st.GetEntry(ptype, Zp, ssF); get_stencil(SE_mem,SE);
    GPU_COALESCED_STENCIL_LEG_PROJ(Zp,spProjZp);
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Zp);
    accumReconZp(result, Uchi);

    SE_mem = st.GetEntry(ptype, Tp, ssF); get_stencil(SE_mem,SE);
    GPU_COALESCED_STENCIL_LEG_PROJ(Tp,spProjTp);
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Tp);
    accumReconTp(result, Uchi);

    SE_mem = st.GetEntry(ptype, Xm, ssF); get_stencil(SE_mem,SE);
    GPU_COALESCED_STENCIL_LEG_PROJ(Xm,spProjXm);
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Xm);
    accumReconXm(result, Uchi);

    SE_mem = st.GetEntry(ptype, Ym, ssF); get_stencil(SE_mem,SE);
    GPU_COALESCED_STENCIL_LEG_PROJ(Ym,spProjYm);
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Ym);
    accumReconYm(result, Uchi);

    SE_mem = st.GetEntry(ptype, Zm, ssF); get_stencil(SE_mem,SE);
    GPU_COALESCED_STENCIL_LEG_PROJ(Zm,spProjZm);
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Zm);
    accumReconZm(result, Uchi);

    SE_mem = st.GetEntry(ptype, Tm, ssF); get_stencil(SE_mem,SE);
    GPU_COALESCED_STENCIL_LEG_PROJ(Tm,spProjTm); 
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Tm);
    accumReconTm(result, Uchi);

    synchronise();
#ifdef GPU_VEC
    insertLane (lane,out[sF],result);
#else
  vstream(out[sF], result);
#endif
  }
}

template <class Impl>
accelerator_inline void WilsonKernels<Impl>::GpuDhopSite(StencilView &st, SiteDoubledGaugeField &U,
						  SiteHalfSpinor *buf,  int Ls, int s,
						  int sU, const FermionFieldView &in, FermionFieldView &out) 
{
#ifdef GPU_VEC
  typename SiteHalfSpinor::scalar_object chi;
  typename SiteHalfSpinor::scalar_object Uchi;
  typename SiteSpinor::scalar_object   result;
#else 
  SiteHalfSpinor chi;
  SiteHalfSpinor Uchi;
  SiteHalfSpinor tmp;
  SiteSpinor   result;
#endif
  typedef typename SiteSpinor::scalar_type scalar_type;
  typedef typename SiteSpinor::vector_type vector_type;
  constexpr int Nsimd = sizeof(vector_type)/sizeof(scalar_type);

  uint64_t lane_offset= get_my_lane_offset(Nsimd);
  uint64_t lanes      = get_my_lanes(Nsimd);

  //  printf (" sU %d s %d Nsimd %d lanes %ld lane_off %ld\n",sU, s, Nsimd, lanes, lane_offset);

  StencilEntry *SE_mem;
  StencilEntry SE;
  int ptype;
  // Forces some degree of coalesce on the table look ups
  // Could also use wide load instructions on the data structure
  uint64_t ssF = Ls * sU;
  uint64_t sF  = ssF + s;

#ifndef __CUDA_ARCH__
  for(int lane = lane_offset;lane<lane_offset+lanes;lane++){
#else
  int lane = lane_offset; {
#endif
    SE_mem = st.GetEntry(ptype, Xp, ssF); get_stencil(SE_mem,SE);
    GPU_COALESCED_STENCIL_LEG_PROJ(Xp,spProjXm); 
    Impl::multLinkGpu(lane,Uchi,U,chi,Xp);
    spReconXm(result, Uchi);

    SE_mem = st.GetEntry(ptype, Yp, ssF); get_stencil(SE_mem,SE);
    GPU_COALESCED_STENCIL_LEG_PROJ(Yp,spProjYm);
    Impl::multLinkGpu(lane,Uchi,U,chi,Yp);
    accumReconYm(result, Uchi);
      
    SE_mem = st.GetEntry(ptype, Zp, ssF); get_stencil(SE_mem,SE);
    GPU_COALESCED_STENCIL_LEG_PROJ(Zp,spProjZm);
    Impl::multLinkGpu(lane,Uchi,U,chi,Zp);
    accumReconZm(result, Uchi);

    SE_mem = st.GetEntry(ptype, Tp, ssF); get_stencil(SE_mem,SE);
    GPU_COALESCED_STENCIL_LEG_PROJ(Tp,spProjTm);
    Impl::multLinkGpu(lane,Uchi,U,chi,Tp);
    accumReconTm(result, Uchi);

    SE_mem = st.GetEntry(ptype, Xm, ssF); get_stencil(SE_mem,SE);
    GPU_COALESCED_STENCIL_LEG_PROJ(Xm,spProjXp);
    Impl::multLinkGpu(lane,Uchi,U,chi,Xm);
    accumReconXp(result, Uchi);

    SE_mem = st.GetEntry(ptype, Ym, ssF); get_stencil(SE_mem,SE);
    GPU_COALESCED_STENCIL_LEG_PROJ(Ym,spProjYp);
    Impl::multLinkGpu(lane,Uchi,U,chi,Ym);
    accumReconYp(result, Uchi);

    SE_mem = st.GetEntry(ptype, Zm, ssF); get_stencil(SE_mem,SE);
    GPU_COALESCED_STENCIL_LEG_PROJ(Zm,spProjZp);
    Impl::multLinkGpu(lane,Uchi,U,chi,Zm);
    accumReconZp(result, Uchi);

    SE_mem = st.GetEntry(ptype, Tm, ssF); get_stencil(SE_mem,SE);
    GPU_COALESCED_STENCIL_LEG_PROJ(Tm,spProjTp); 
    Impl::multLinkGpu(lane,Uchi,U,chi,Tm);
    accumReconTp(result, Uchi);

    synchronise();
#ifdef GPU_VEC
    insertLane (lane,out[sF],result);
#else
  vstream(out[sF], result);
#endif
  }

};

// Template specialise Gparity to empty for now
#define GPU_EMPTY(A)							\
  template <>								\
accelerator_inline void							\
WilsonKernels<A>::GpuDhopSite(StencilView &st,				\
			      SiteDoubledGaugeField &U,			\
			      SiteHalfSpinor *buf, int Ls, int sF,	\
			      int sU,					\
			      const FermionFieldView &in,		\
			      FermionFieldView &out) { assert(0);};	\
  template <>								\
  accelerator_inline void							\
  WilsonKernels<A>::GpuDhopSiteDag(StencilView &st,			\
				   DoubledGaugeFieldView &U,		\
				   SiteHalfSpinor *buf, int Ls,int sF,	\
				   int sU,				\
				   const FermionFieldView &in,		\
				   FermionFieldView &out) { assert(0);};

GPU_EMPTY(GparityWilsonImplF);
GPU_EMPTY(GparityWilsonImplFH);
GPU_EMPTY(GparityWilsonImplD);
GPU_EMPTY(GparityWilsonImplDF);

template <class Impl>
void WilsonKernels<Impl>::Dhop(int Opt,StencilImpl &st,  DoubledGaugeField &U, SiteHalfSpinor * buf,
			       int Ls, int Nsite, const FermionField &in, FermionField &out,
			       int interior,int exterior) 
{
    auto U_v   = U.View();
    auto in_v  = in.View();
    auto out_v = out.View();
    auto st_v  = st.View();
    if ( (Opt == WilsonKernelsStatic::OptGpu) && interior && exterior ) { 
      const uint64_t nsimd = Simd::Nsimd();
      const uint64_t    NN = Nsite*Ls*nsimd;
      accelerator_loopN( sss, NN, {
	  uint64_t cur  = sss;
	  //	  uint64_t lane = cur % nsimd;
	  cur = cur / nsimd;
	  uint64_t   s  = cur%Ls;
	  uint64_t   sF = cur;         cur = cur / Ls;
	  uint64_t   sU = cur;
	  WilsonKernels<Impl>::GpuDhopSite(st_v,U_v[sU],buf,Ls,s,sU,in_v,out_v);
      });
    } else { 
      accelerator_loop( ss, U_v, {
	int sU = ss;
        int sF = Ls * sU;
        DhopSite(Opt,st_v,U_v,st.CommBuf(),sF,sU,Ls,1,in_v,out_v);
      });
    }
  }
  template <class Impl>
  void WilsonKernels<Impl>::DhopDag(int Opt,StencilImpl &st,  DoubledGaugeField &U, SiteHalfSpinor * buf,
				    int Ls, int Nsite, const FermionField &in, FermionField &out,
				    int interior,int exterior) 
  {
    auto U_v   = U.View();
    auto in_v  = in.View();
    auto out_v = out.View();
    auto st_v  = st.View();

    if ( (Opt == WilsonKernelsStatic::OptGpu) && interior && exterior ) { 
      const uint64_t nsimd = Simd::Nsimd();
      const uint64_t    NN = Nsite*Ls*nsimd;
      accelerator_loopN( sss, NN, {
	  uint64_t cur  = sss;
	  // uint64_t lane = cur % nsimd;
	  cur = cur / nsimd;
	  uint64_t   s  = cur%Ls;
	  uint64_t   sF = cur;         cur = cur / Ls;
	  uint64_t   sU = cur;
	  WilsonKernels<Impl>::GpuDhopSiteDag(st_v,U_v,buf,Ls,s,sU,in_v,out_v);
      });
    } else { 
      accelerator_loop( ss, U_v, {
	int sU = ss;
        int sF = Ls * sU;
        DhopSiteDag(Opt,st,U_v,st.CommBuf(),sF,sU,Ls,1,in_v,out_v);
      });
    }
  }


/*
GPU_EMPTY(DomainWallVec5dImplF);
GPU_EMPTY(DomainWallVec5dImplFH);
GPU_EMPTY(DomainWallVec5dImplD);
GPU_EMPTY(DomainWallVec5dImplDF);
GPU_EMPTY(ZDomainWallVec5dImplF);
GPU_EMPTY(ZDomainWallVec5dImplFH);
GPU_EMPTY(ZDomainWallVec5dImplD);
GPU_EMPTY(ZDomainWallVec5dImplDF);
*/

FermOpTemplateInstantiate(WilsonKernels);
AdjointFermOpTemplateInstantiate(WilsonKernels);
TwoIndexFermOpTemplateInstantiate(WilsonKernels);

NAMESPACE_END(Grid);

