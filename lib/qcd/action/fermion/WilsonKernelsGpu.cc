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


#define GPU_COALESCED_STENCIL_LEG_PROJ(Dir,spProj)			\
  synchronise();							\
  if (SE->_is_local) {							\
    int mask = Nsimd >> (ptype + 1);					\
    int plane= SE->_permute ? (lane ^ mask) : lane;			\
    auto in_l = extractLane(plane,in[SE->_offset]);			\
    spProj(chi,in_l);							\
  } else {								\
    chi  = extractLane(lane,buf[SE->_offset]);				\
  }									\
  synchronise();

template <class Impl>
accelerator void WilsonKernels<Impl>::GpuDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U,
						     SiteHalfSpinor *buf, int sF,
						     int sU, const FermionFieldView &in, FermionFieldView &out)
{
  typename SiteHalfSpinor::scalar_object chi;
  typename SiteHalfSpinor::scalar_object Uchi;
  typename SiteSpinor::scalar_object   result;
  typedef typename SiteSpinor::scalar_type scalar_type;
  typedef typename SiteSpinor::vector_type vector_type;

  constexpr int Nsimd = sizeof(vector_type)/sizeof(scalar_type);

  uint64_t lane_offset= get_my_lane_offset(Nsimd);
  uint64_t lanes      = get_my_lanes(Nsimd);

  StencilEntry *SE;
  int ptype;

#ifndef __CUDA_ARCH__
  for(int lane = lane_offset;lane<lane_offset+lanes;lane++){
#else
  int lane = lane_offset; {
#endif
    SE = st.GetEntry(ptype, Xp, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Xp,spProjXp); 
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Xp);
    spReconXp(result, Uchi);

    SE = st.GetEntry(ptype, Yp, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Yp,spProjYp);
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Yp);
    accumReconYp(result, Uchi);
      
    SE = st.GetEntry(ptype, Zp, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Zp,spProjZp);
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Zp);
    accumReconZp(result, Uchi);

    SE = st.GetEntry(ptype, Tp, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Tp,spProjTp);
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Tp);
    accumReconTp(result, Uchi);

    SE = st.GetEntry(ptype, Xm, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Xm,spProjXm);
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Xm);
    accumReconXm(result, Uchi);

    SE = st.GetEntry(ptype, Ym, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Ym,spProjYm);
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Ym);
    accumReconYm(result, Uchi);


    SE = st.GetEntry(ptype, Zm, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Zm,spProjZm);
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Zm);
    accumReconZm(result, Uchi);

    SE = st.GetEntry(ptype, Tm, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Tm,spProjTm); 
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Tm);
    accumReconTm(result, Uchi);

    synchronise();
    insertLane (lane,out[sF],result);
  }
}

template <class Impl>
accelerator void WilsonKernels<Impl>::GpuDhopSite(StencilView &st, DoubledGaugeFieldView &U,
						  SiteHalfSpinor *buf, int sF,
						  int sU, const FermionFieldView &in, FermionFieldView &out) 
{
  typename SiteHalfSpinor::scalar_object chi;
  typename SiteHalfSpinor::scalar_object Uchi;
  typename SiteSpinor::scalar_object   result;
  typedef typename SiteSpinor::scalar_type scalar_type;
  typedef typename SiteSpinor::vector_type vector_type;

  constexpr int Nsimd = sizeof(vector_type)/sizeof(scalar_type);

  uint64_t lane_offset= get_my_lane_offset(Nsimd);
  uint64_t lanes      = get_my_lanes(Nsimd);

  StencilEntry *SE;
  int ptype;

#ifndef __CUDA_ARCH__
  for(int lane = lane_offset;lane<lane_offset+lanes;lane++){
#else
  int lane = lane_offset; {
#endif
    SE = st.GetEntry(ptype, Xp, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Xp,spProjXm); 
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Xp);
    spReconXm(result, Uchi);

    SE = st.GetEntry(ptype, Yp, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Yp,spProjYm);
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Yp);
    accumReconYm(result, Uchi);
      
    SE = st.GetEntry(ptype, Zp, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Zp,spProjZm);
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Zp);
    accumReconZm(result, Uchi);

    SE = st.GetEntry(ptype, Tp, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Tp,spProjTm);
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Tp);
    accumReconTm(result, Uchi);

    SE = st.GetEntry(ptype, Xm, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Xm,spProjXp);
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Xm);
    accumReconXp(result, Uchi);

    SE = st.GetEntry(ptype, Ym, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Ym,spProjYp);
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Ym);
    accumReconYp(result, Uchi);

    SE = st.GetEntry(ptype, Zm, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Zm,spProjZp);
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Zm);
    accumReconZp(result, Uchi);

    SE = st.GetEntry(ptype, Tm, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Tm,spProjTp); 
    Impl::multLinkGpu(lane,Uchi,U[sU],chi,Tm);
    accumReconTp(result, Uchi);

    synchronise();
    insertLane (lane,out[sF],result);
  }

};

// Template specialise Gparity to empty for now
#define GPU_EMPTY(A)							\
  template <>								\
accelerator void							\
WilsonKernels<A>::GpuDhopSite(StencilView &st,				\
			      DoubledGaugeFieldView &U,			\
			      SiteHalfSpinor *buf, int sF,		\
			      int sU,					\
			      const FermionFieldView &in,		\
			      FermionFieldView &out) { assert(0);};	\
  template <>								\
  accelerator void							\
  WilsonKernels<A>::GpuDhopSiteDag(StencilView &st,			\
				DoubledGaugeFieldView &U,		\
				   SiteHalfSpinor *buf, int sF,		\
				int sU,					\
				const FermionFieldView &in,		\
				   FermionFieldView &out) { assert(0);};

GPU_EMPTY(GparityWilsonImplF);
GPU_EMPTY(GparityWilsonImplFH);
GPU_EMPTY(GparityWilsonImplD);
GPU_EMPTY(GparityWilsonImplDF);

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

