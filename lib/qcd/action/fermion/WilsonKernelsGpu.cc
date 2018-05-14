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
// Gpu implementation; view code at a premium and less unroll
//////////////////////////////////////////////////////////////
  
#define GPU_STENCIL_LEG_PROJ(Dir,spProj)			\
  if (SE->_is_local) {						\
    spProj(chi, in[SE->_offset]);				\
    if (SE->_permute) {						\
      permute(tmp, chi, ptype);					\
      chi = tmp;						\
    }								\
  } else {							\
    chi = buf[SE->_offset];					\
  }								

#define GPU_STENCIL_LEG_RECON(Recon)  Recon(result, Uchi);
// Xp is mu= 0
template <class Impl>
accelerator void WilsonKernels<Impl>::GpuDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U,
						     SiteHalfSpinor *buf, int sF,
						     int sU, const FermionFieldView &in, FermionFieldView &out)
{
  SiteHalfSpinor tmp;
  SiteHalfSpinor chi;
  SiteHalfSpinor Uchi;
  SiteSpinor result;
  StencilEntry *SE;
  int ptype;

  for(int mu=0;mu<2*Nd;mu++) {

    SE = st.GetEntry(ptype, mu, sF);

    switch(mu){
    case Xp:
      GPU_STENCIL_LEG_PROJ(Xp,spProjXp); break;
    case Yp:
      GPU_STENCIL_LEG_PROJ(Yp,spProjYp); break;
    case Zp:
      GPU_STENCIL_LEG_PROJ(Zp,spProjZp); break;
    case Tp:
      GPU_STENCIL_LEG_PROJ(Tp,spProjTp); break;
    case Xm:
      GPU_STENCIL_LEG_PROJ(Xm,spProjXm); break;
    case Ym:
      GPU_STENCIL_LEG_PROJ(Ym,spProjYm); break;
    case Zm:
      GPU_STENCIL_LEG_PROJ(Zm,spProjZm); break;
    case Tm:
    default:
      GPU_STENCIL_LEG_PROJ(Tm,spProjTm); break;
    }

    Impl::multLink(Uchi, U[sU], chi, mu, SE, st);	

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
  }
  vstream(out[sF], result);
};

template <class Impl>
accelerator void WilsonKernels<Impl>::GpuDhopSite(StencilView &st, DoubledGaugeFieldView &U,
						  SiteHalfSpinor *buf, int sF,
						  int sU, const FermionFieldView &in, FermionFieldView &out) 
{
  SiteHalfSpinor tmp;
  SiteHalfSpinor chi;
  SiteHalfSpinor Uchi;
  SiteSpinor result;
  StencilEntry *SE;
  int ptype;
  
  for(int mu=0;mu<2*Nd;mu++) {
 
    SE = st.GetEntry(ptype, mu, sF);

    switch(mu){
    case Xp:
      GPU_STENCIL_LEG_PROJ(Xp,spProjXm); break;
    case Yp:
      GPU_STENCIL_LEG_PROJ(Yp,spProjYm); break;
    case Zp:
      GPU_STENCIL_LEG_PROJ(Zp,spProjZm); break;
    case Tp:
      GPU_STENCIL_LEG_PROJ(Tp,spProjTm); break;
    case Xm:
      GPU_STENCIL_LEG_PROJ(Xm,spProjXp); break;
    case Ym:
      GPU_STENCIL_LEG_PROJ(Ym,spProjYp); break;
    case Zm:
      GPU_STENCIL_LEG_PROJ(Zm,spProjZp); break;
    case Tm:
    default:
      GPU_STENCIL_LEG_PROJ(Tm,spProjTp); break;
    }

    Impl::multLink(Uchi, U[sU], chi, mu, SE, st);	

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
  vstream(out[sF], result);

};

FermOpTemplateInstantiate(WilsonKernels);

NAMESPACE_END(Grid);

