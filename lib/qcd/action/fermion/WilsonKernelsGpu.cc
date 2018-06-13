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
// Gpu implementation; thread loop is implicit
//////////////////////////////////////////////////////////////
__host__ __device__ inline void synchronise(void) 
{
#ifdef __CUDA_ARCH__
  __syncthreads();
#endif
  return;
}

#define GPU_DSLASH_COALESCE  
#ifdef GPU_DSLASH_COALESCE

__host__ __device__ inline int get_my_lanes(int Nsimd) 
{
#ifdef __CUDA_ARCH__
  return 1;
#else 
  return Nsimd;
#endif
}
__host__ __device__ inline int get_my_lane_offset(int Nsimd) 
{
#ifdef __CUDA_ARCH__
  return ( (threadIdx.x) % Nsimd);
#else
  return 0;
#endif
}

////////////////////////////////////////////////////////////////////////
// Extract/Insert a single lane; do this locally in this file.
// Don't need a global version really.
////////////////////////////////////////////////////////////////////////
template<class vobj> accelerator_inline
typename vobj::scalar_object extractLaneGpu(int lane, const vobj & __restrict__ vec)
{
  typedef typename vobj::scalar_object scalar_object;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  constexpr int words=sizeof(vobj)/sizeof(vector_type);
  constexpr int Nsimd=vector_type::Nsimd();

  scalar_object extracted;
  scalar_type * __restrict__  sp = (scalar_type *)&extracted; // Type pun
  scalar_type * __restrict__  vp = (scalar_type *)&vec;
  for(int w=0;w<words;w++){
    sp[w]=vp[w*Nsimd+lane];
  }
  return extracted;
}

template<class vobj> accelerator_inline
void insertLaneFloat2(int lane, vobj & __restrict__ vec,const typename vobj::scalar_object & __restrict__ extracted)
{
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  constexpr int words=sizeof(vobj)/sizeof(vector_type);
  constexpr int Nsimd=vector_type::Nsimd();

  float2 * __restrict__ sp = (float2 *)&extracted;
  float2 * __restrict__ vp = (float2 *)&vec;
  for(int w=0;w<words;w++){
    vp[w*Nsimd+lane]=sp[w];
  }
}
template<class vobj> accelerator_inline
typename vobj::scalar_object extractLaneFloat2(int lane, const vobj & __restrict__ vec)
{
  typedef typename vobj::scalar_object scalar_object;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  constexpr int words=sizeof(vobj)/sizeof(vector_type);
  constexpr int Nsimd=vector_type::Nsimd();

  scalar_object extracted;
  float2 * __restrict__  sp = (float2 *)&extracted; // Type pun
  float2 * __restrict__  vp = (float2 *)&vec;
  for(int w=0;w<words;w++){
    sp[w]=vp[w*Nsimd+lane];
  }
  return extracted;
}


#define GPU_COALESCED_STENCIL_LEG_PROJ(Dir,spProj)		\
  if (SE->_is_local) {						\
    auto in_l = extractLaneGpu(lane,in[SE->_offset]);		\
    spProj(chi,in_l);						\
  } else {							\
    chi  = extractLaneGpu(lane,buf[SE->_offset]);		\
  }								

template <class Impl>
accelerator void WilsonKernels<Impl>::GpuDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U,
						     SiteHalfSpinor *buf, int sF,int LLs,
						     int sU, const FermionFieldView &in, FermionFieldView &out)
{
  typename SiteHalfSpinor::scalar_object tmp;
  typename SiteHalfSpinor::scalar_object chi;
  typename SiteHalfSpinor::scalar_object Uchi;
  typename SiteSpinor::scalar_object   result;

  typedef typename SiteSpinor::scalar_type scalar_type;
  typedef typename SiteSpinor::vector_type vector_type;
  constexpr int Nsimd = sizeof(vector_type)/sizeof(scalar_type);

  uint64_t lane_offset= get_my_lane_offset(Nsimd);
  uint64_t lanes      = get_my_lanes (Nsimd);

  StencilEntry *SE;
  int ptype;

  for(int lane = lane_offset;lane<lane_offset+lanes;lane++){
  for(int s=0;s<LLs;s++){
  for(int mu=0;mu<2*Nd;mu++) {

    SE = st.GetEntry(ptype, mu, sF);
    switch(mu){
    case Xp:
      GPU_COALESCED_STENCIL_LEG_PROJ(Xp,spProjXp); break;
    case Yp:
      GPU_COALESCED_STENCIL_LEG_PROJ(Yp,spProjYp); break;
    case Zp:
      GPU_COALESCED_STENCIL_LEG_PROJ(Zp,spProjZp); break;
    case Tp:
      GPU_COALESCED_STENCIL_LEG_PROJ(Tp,spProjTp); break;
    case Xm:
      GPU_COALESCED_STENCIL_LEG_PROJ(Xm,spProjXm); break;
    case Ym:
      GPU_COALESCED_STENCIL_LEG_PROJ(Ym,spProjYm); break;
    case Zm:
      GPU_COALESCED_STENCIL_LEG_PROJ(Zm,spProjZm); break;
    case Tm:
    default:
      GPU_COALESCED_STENCIL_LEG_PROJ(Tm,spProjTm); break;
    }
    synchronise();

    auto U_l = extractLaneGpu(lane,U[sU](mu));
    Uchi()=U_l*chi();

    switch(mu){
    case Xp:
      spReconXp(result, Uchi);    break;
    case Yp:
      accumReconYp(result, Uchi); break;
    case Zp:
      accumReconZp(result, Uchi); break;
    case Tp:
      accumReconTp(result, Uchi); break;
    case Xm:
      accumReconXm(result, Uchi); break;
    case Ym:
      accumReconYm(result, Uchi); break;
    case Zm:
      accumReconZm(result, Uchi); break;
    case Tm:
    default:
      accumReconTm(result, Uchi); break;
    }
    synchronise();
  }
  insertLane (lane,out[sF],result);
  sF++;
  }}
};

template <class Impl>
accelerator void WilsonKernels<Impl>::GpuDhopSite(StencilView &st, DoubledGaugeFieldView &U,
						  SiteHalfSpinor *buf, int sF,int LLs,
						  int sU, const FermionFieldView &in, FermionFieldView &out) 
{
  typename SiteHalfSpinor::scalar_object tmp;
  typename SiteHalfSpinor::scalar_object chi;
  typename SiteHalfSpinor::scalar_object Uchi;
  typename SiteSpinor::scalar_object   result;
  typedef typename SiteSpinor::scalar_type scalar_type;
  typedef typename SiteSpinor::vector_type vector_type;

  constexpr int Nsimd = sizeof(vector_type)/sizeof(scalar_type);

  uint64_t lane_offset= get_my_lane_offset(Nsimd);
  uint64_t lanes      = get_my_lanes(Nsimd);

  //  printf("Evaluating site %d  Nsimd %d : lanes %ld %ld - %ld\n",sF,Nsimd,lanes,lane_offset,lane_offset+lanes-1);

  StencilEntry *SE;
  int ptype;
  
  for(int lane = lane_offset;lane<lane_offset+lanes;lane++){
  for(int s=0;s<LLs;s++){
#if 1
    int mu=0;
    SE = st.GetEntry(ptype, mu, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Xp,spProjXm); 
    { auto U_l = extractLaneFloat2(lane,U[sU](mu)); Uchi() =  U_l * chi();}
    spReconXm(result, Uchi);

    mu++; SE = st.GetEntry(ptype, mu, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Yp,spProjYm);
    { auto U_l = extractLaneFloat2(lane,U[sU](mu)); Uchi() =  U_l * chi();}
    accumReconYm(result, Uchi);
      
    mu++; SE = st.GetEntry(ptype, mu, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Zp,spProjZm);
    { auto U_l = extractLaneFloat2(lane,U[sU](mu)); Uchi() =  U_l * chi();}
    accumReconZm(result, Uchi);

    mu++; SE = st.GetEntry(ptype, mu, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Tp,spProjTm);
    { auto U_l = extractLaneFloat2(lane,U[sU](mu)); Uchi() =  U_l * chi();}
    accumReconTm(result, Uchi);

    mu++; SE = st.GetEntry(ptype, mu, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Xm,spProjXp);
    { auto U_l = extractLaneFloat2(lane,U[sU](mu)); Uchi() =  U_l * chi();}
    accumReconXp(result, Uchi);

    mu++; SE = st.GetEntry(ptype, mu, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Ym,spProjYp);
    { auto U_l = extractLaneFloat2(lane,U[sU](mu)); Uchi() =  U_l * chi();}
    accumReconYp(result, Uchi);


    mu++; SE = st.GetEntry(ptype, mu, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Zm,spProjZp);
    { auto U_l = extractLaneFloat2(lane,U[sU](mu)); Uchi() =  U_l * chi();}
    accumReconZp(result, Uchi);

    mu++; SE = st.GetEntry(ptype, mu, sF);
    GPU_COALESCED_STENCIL_LEG_PROJ(Tm,spProjTp); 
    { auto U_l = extractLaneFloat2(lane,U[sU](mu)); Uchi() =  U_l * chi();}
    accumReconTp(result, Uchi);

#else 
  for(int mu=0;mu<2*Nd;mu++) {
 
    SE = st.GetEntry(ptype, mu, sF);

    switch(mu){
    case Xp:
      GPU_COALESCED_STENCIL_LEG_PROJ(Xp,spProjXm); break;
    case Yp:
      GPU_COALESCED_STENCIL_LEG_PROJ(Yp,spProjYm); break;
    case Zp:
      GPU_COALESCED_STENCIL_LEG_PROJ(Zp,spProjZm); break;
    case Tp:
      GPU_COALESCED_STENCIL_LEG_PROJ(Tp,spProjTm); break;
    case Xm:
      GPU_COALESCED_STENCIL_LEG_PROJ(Xm,spProjXp); break;
    case Ym:
      GPU_COALESCED_STENCIL_LEG_PROJ(Ym,spProjYp); break;
    case Zm:
      GPU_COALESCED_STENCIL_LEG_PROJ(Zm,spProjZp); break;
    case Tm:
    default:
      GPU_COALESCED_STENCIL_LEG_PROJ(Tm,spProjTp); break;
    }

    auto U_l = extractLaneGpu(lane,U[sU](mu));

    auto tmp =  U_l * chi();

    Uchi() = tmp;

    switch(mu){
    case Xp:
      spReconXm(result, Uchi); break;
    case Yp:
      accumReconYm(result, Uchi); break;
    case Zp:
      accumReconZm(result, Uchi); break;
    case Tp:
      accumReconTm(result, Uchi); break;
    case Xm:
      accumReconXp(result, Uchi); break;
    case Ym:
      accumReconYp(result, Uchi); break;
    case Zm:
      accumReconZp(result, Uchi); break;
    case Tm:
    default:
      accumReconTp(result, Uchi); break;
    }
  }
#endif
  insertLaneFloat2 (lane,out[sF],result);
  sF++;
  }}
};

// Template specialise Gparity to empty for now
#define GPU_EMPTY(A)							\
  template <>								\
accelerator void							\
WilsonKernels<A>::GpuDhopSite(StencilView &st,				\
			      DoubledGaugeFieldView &U,			\
			      SiteHalfSpinor *buf, int sF, int LLs,	\
			      int sU,					\
			      const FermionFieldView &in,		\
			      FermionFieldView &out) {};		\
  template <>								\
  accelerator void							\
  WilsonKernels<A>::GpuDhopSiteDag(StencilView &st,			\
				DoubledGaugeFieldView &U,		\
				   SiteHalfSpinor *buf, int sF,int LLs,	\
				int sU,					\
				const FermionFieldView &in,		\
				FermionFieldView &out) {};

GPU_EMPTY(GparityWilsonImplF);
GPU_EMPTY(GparityWilsonImplFH);
GPU_EMPTY(GparityWilsonImplD);
GPU_EMPTY(GparityWilsonImplDF);
GPU_EMPTY(DomainWallVec5dImplF);
GPU_EMPTY(DomainWallVec5dImplFH);
GPU_EMPTY(DomainWallVec5dImplD);
GPU_EMPTY(DomainWallVec5dImplDF);
GPU_EMPTY(ZDomainWallVec5dImplF);
GPU_EMPTY(ZDomainWallVec5dImplFH);
GPU_EMPTY(ZDomainWallVec5dImplD);
GPU_EMPTY(ZDomainWallVec5dImplDF);

FermOpTemplateInstantiate(WilsonKernels);

NAMESPACE_END(Grid);

