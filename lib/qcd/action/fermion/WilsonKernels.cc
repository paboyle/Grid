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
#include <Grid.h>
namespace Grid {
namespace QCD {

  int WilsonKernelsStatic::Opt   = WilsonKernelsStatic::OptGeneric;
  int WilsonKernelsStatic::Comms = WilsonKernelsStatic::CommsAndCompute;

#ifdef QPX
#include <spi/include/kernel/location.h>
#include <spi/include/l1p/types.h>
#include <hwi/include/bqc/l1p_mmio.h>
#include <hwi/include/bqc/A2_inlines.h>
#endif

void bgq_l1p_optimisation(int mode)
{
#ifdef QPX
#undef L1P_CFG_PF_USR
#define L1P_CFG_PF_USR  (0x3fde8000108ll)   /*  (64 bit reg, 23 bits wide, user/unpriv) */

  uint64_t cfg_pf_usr;
  if ( mode ) { 
    cfg_pf_usr =
        L1P_CFG_PF_USR_ifetch_depth(0)       
      | L1P_CFG_PF_USR_ifetch_max_footprint(1)   
      | L1P_CFG_PF_USR_pf_stream_est_on_dcbt 
      | L1P_CFG_PF_USR_pf_stream_establish_enable
      | L1P_CFG_PF_USR_pf_stream_optimistic
      | L1P_CFG_PF_USR_pf_adaptive_throttle(0xF) ;
    //    if ( sizeof(Float) == sizeof(double) ) {
      cfg_pf_usr |=  L1P_CFG_PF_USR_dfetch_depth(2)| L1P_CFG_PF_USR_dfetch_max_footprint(3)   ;
      //    } else {
      //      cfg_pf_usr |=  L1P_CFG_PF_USR_dfetch_depth(1)| L1P_CFG_PF_USR_dfetch_max_footprint(2)   ;
      //    }
  } else { 
    cfg_pf_usr = L1P_CFG_PF_USR_dfetch_depth(1)
      | L1P_CFG_PF_USR_dfetch_max_footprint(2)   
      | L1P_CFG_PF_USR_ifetch_depth(0)       
      | L1P_CFG_PF_USR_ifetch_max_footprint(1)   
      | L1P_CFG_PF_USR_pf_stream_est_on_dcbt 
      | L1P_CFG_PF_USR_pf_stream_establish_enable
      | L1P_CFG_PF_USR_pf_stream_optimistic
      | L1P_CFG_PF_USR_pf_stream_prefetch_enable;
  }
  *((uint64_t *)L1P_CFG_PF_USR) = cfg_pf_usr;

#endif

}


template <class Impl>
WilsonKernels<Impl>::WilsonKernels(const ImplParams &p) : Base(p){};

////////////////////////////////////////////
// Generic implementation; move to different file?
////////////////////////////////////////////

template <class Impl>
void WilsonKernels<Impl>::GenericDhopSiteDag(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U,
						     SiteHalfSpinor *buf, int sF,
						     int sU, const FermionField &in, FermionField &out,
						     int interior,int exterior) {
  SiteHalfSpinor tmp;
  SiteHalfSpinor chi;
  SiteHalfSpinor *chi_p;
  SiteHalfSpinor Uchi;
  SiteSpinor result;
  StencilEntry *SE;
  int ptype;

  ///////////////////////////
  // Xp
  ///////////////////////////
  SE = st.GetEntry(ptype, Xp, sF);

  if (SE->_is_local) {
    chi_p = &chi;
    if (SE->_permute) {
      spProjXp(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else {
      spProjXp(chi, in._odata[SE->_offset]);
    }
  } else {
    chi_p = &buf[SE->_offset];
  }

  Impl::multLink(Uchi, U._odata[sU], *chi_p, Xp, SE, st);
  spReconXp(result, Uchi);

  ///////////////////////////
  // Yp
  ///////////////////////////
  SE = st.GetEntry(ptype, Yp, sF);

  if (SE->_is_local) {
    chi_p = &chi;
    if (SE->_permute) {
      spProjYp(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else {
      spProjYp(chi, in._odata[SE->_offset]);
    }
  } else {
    chi_p = &buf[SE->_offset];
  }

  Impl::multLink(Uchi, U._odata[sU], *chi_p, Yp, SE, st);
  accumReconYp(result, Uchi);

  ///////////////////////////
  // Zp
  ///////////////////////////
  SE = st.GetEntry(ptype, Zp, sF);

  if (SE->_is_local) {
    chi_p = &chi;
    if (SE->_permute) {
      spProjZp(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else {
      spProjZp(chi, in._odata[SE->_offset]);
    }
  } else {
    chi_p = &buf[SE->_offset];
  }

  Impl::multLink(Uchi, U._odata[sU], *chi_p, Zp, SE, st);
  accumReconZp(result, Uchi);

  ///////////////////////////
  // Tp
  ///////////////////////////
  SE = st.GetEntry(ptype, Tp, sF);

  if (SE->_is_local) {
    chi_p = &chi;
    if (SE->_permute) {
      spProjTp(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else {
      spProjTp(chi, in._odata[SE->_offset]);
    }
  } else {
    chi_p = &buf[SE->_offset];
  }

  Impl::multLink(Uchi, U._odata[sU], *chi_p, Tp, SE, st);
  accumReconTp(result, Uchi);

  ///////////////////////////
  // Xm
  ///////////////////////////
  SE = st.GetEntry(ptype, Xm, sF);

  if (SE->_is_local) {
    chi_p = &chi;
    if (SE->_permute) {
      spProjXm(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else {
      spProjXm(chi, in._odata[SE->_offset]);
    }
  } else {
    chi_p = &buf[SE->_offset];
  }

  Impl::multLink(Uchi, U._odata[sU], *chi_p, Xm, SE, st);
  accumReconXm(result, Uchi);

  ///////////////////////////
  // Ym
  ///////////////////////////
  SE = st.GetEntry(ptype, Ym, sF);

  if (SE->_is_local) {
    chi_p = &chi;
    if (SE->_permute) {
      spProjYm(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else {
      spProjYm(chi, in._odata[SE->_offset]);
    }
  } else {
    chi_p = &buf[SE->_offset];
  }

  Impl::multLink(Uchi, U._odata[sU], *chi_p, Ym, SE, st);
  accumReconYm(result, Uchi);

  ///////////////////////////
  // Zm
  ///////////////////////////
  SE = st.GetEntry(ptype, Zm, sF);

  if (SE->_is_local) {
    chi_p = &chi;
    if (SE->_permute) {
      spProjZm(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else {
      spProjZm(chi, in._odata[SE->_offset]);
    }
  } else {
    chi_p = &buf[SE->_offset];
  }

  Impl::multLink(Uchi, U._odata[sU], *chi_p, Zm, SE, st);
  accumReconZm(result, Uchi);

  ///////////////////////////
  // Tm
  ///////////////////////////
  SE = st.GetEntry(ptype, Tm, sF);

  if (SE->_is_local) {
    chi_p = &chi;
    if (SE->_permute) {
      spProjTm(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else {
      spProjTm(chi, in._odata[SE->_offset]);
    }
  } else {
    chi_p = &buf[SE->_offset];
  }

  Impl::multLink(Uchi, U._odata[sU], *chi_p, Tm, SE, st);
  accumReconTm(result, Uchi);

  vstream(out._odata[sF], result);
};

// Need controls to do interior, exterior, or both
template <class Impl>
void WilsonKernels<Impl>::GenericDhopSite(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U,
						  SiteHalfSpinor *buf, int sF,
						  int sU, const FermionField &in, FermionField &out,int interior,int exterior) {
  SiteHalfSpinor tmp;
  SiteHalfSpinor chi;
  SiteHalfSpinor *chi_p;
  SiteHalfSpinor Uchi;
  SiteSpinor result;
  StencilEntry *SE;
  int ptype;

  ///////////////////////////
  // Xp
  ///////////////////////////
  SE = st.GetEntry(ptype, Xm, sF);

  if (SE->_is_local) {
    chi_p = &chi;
    if (SE->_permute) {
      spProjXp(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else {
      spProjXp(chi, in._odata[SE->_offset]);
    }
  } else {
    chi_p = &buf[SE->_offset];
  }

  Impl::multLink(Uchi, U._odata[sU], *chi_p, Xm, SE, st);
  spReconXp(result, Uchi);

  ///////////////////////////
  // Yp
  ///////////////////////////
  SE = st.GetEntry(ptype, Ym, sF);

  if (SE->_is_local) {
    chi_p = &chi;
    if (SE->_permute) {
      spProjYp(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else {
      spProjYp(chi, in._odata[SE->_offset]);
    }
  } else {
    chi_p = &buf[SE->_offset];
  }

  Impl::multLink(Uchi, U._odata[sU], *chi_p, Ym, SE, st);
  accumReconYp(result, Uchi);

  ///////////////////////////
  // Zp
  ///////////////////////////
  SE = st.GetEntry(ptype, Zm, sF);

  if (SE->_is_local) {
    chi_p = &chi;
    if (SE->_permute) {
      spProjZp(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else {
      spProjZp(chi, in._odata[SE->_offset]);
    }
  } else {
    chi_p = &buf[SE->_offset];
  }

  Impl::multLink(Uchi, U._odata[sU], *chi_p, Zm, SE, st);
  accumReconZp(result, Uchi);

  ///////////////////////////
  // Tp
  ///////////////////////////
  SE = st.GetEntry(ptype, Tm, sF);

  if (SE->_is_local) {
    chi_p = &chi;
    if (SE->_permute) {
      spProjTp(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else {
      spProjTp(chi, in._odata[SE->_offset]);
    }
  } else {
    chi_p = &buf[SE->_offset];
  }

  Impl::multLink(Uchi, U._odata[sU], *chi_p, Tm, SE, st);
  accumReconTp(result, Uchi);

  ///////////////////////////
  // Xm
  ///////////////////////////
  SE = st.GetEntry(ptype, Xp, sF);

  if (SE->_is_local) {
    chi_p = &chi;
    if (SE->_permute) {
      spProjXm(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else {
      spProjXm(chi, in._odata[SE->_offset]);
    }
  } else {
    chi_p = &buf[SE->_offset];
  }

  Impl::multLink(Uchi, U._odata[sU], *chi_p, Xp, SE, st);
  accumReconXm(result, Uchi);

  ///////////////////////////
  // Ym
  ///////////////////////////
  SE = st.GetEntry(ptype, Yp, sF);

  if (SE->_is_local) {
    chi_p = &chi;
    if (SE->_permute) {
      spProjYm(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else {
      spProjYm(chi, in._odata[SE->_offset]);
    }
  } else {
    chi_p = &buf[SE->_offset];
  }

  Impl::multLink(Uchi, U._odata[sU], *chi_p, Yp, SE, st);
  accumReconYm(result, Uchi);

  ///////////////////////////
  // Zm
  ///////////////////////////
  SE = st.GetEntry(ptype, Zp, sF);

  if (SE->_is_local) {
    chi_p = &chi;
    if (SE->_permute) {
      spProjZm(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else {
      spProjZm(chi, in._odata[SE->_offset]);
    }
  } else {
    chi_p = &buf[SE->_offset];
  }

  Impl::multLink(Uchi, U._odata[sU], *chi_p, Zp, SE, st);
  accumReconZm(result, Uchi);

  ///////////////////////////
  // Tm
  ///////////////////////////
  SE = st.GetEntry(ptype, Tp, sF);

  if (SE->_is_local) {
    chi_p = &chi;
    if (SE->_permute) {
      spProjTm(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else {
      spProjTm(chi, in._odata[SE->_offset]);
    }
  } else {
    chi_p = &buf[SE->_offset];
  }

  Impl::multLink(Uchi, U._odata[sU], *chi_p, Tp, SE, st);
  accumReconTm(result, Uchi);

  vstream(out._odata[sF], result);
};

template <class Impl>
void WilsonKernels<Impl>::DhopDir( StencilImpl &st, DoubledGaugeField &U,SiteHalfSpinor *buf, int sF,
					   int sU, const FermionField &in, FermionField &out, int dir, int gamma) {

  SiteHalfSpinor tmp;
  SiteHalfSpinor chi;
  SiteSpinor result;
  SiteHalfSpinor Uchi;
  StencilEntry *SE;
  int ptype;

  SE = st.GetEntry(ptype, dir, sF);

  // Xp
  if (gamma == Xp) {
    if (SE->_is_local && SE->_permute) {
      spProjXp(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else if (SE->_is_local) {
      spProjXp(chi, in._odata[SE->_offset]);
    } else {
      chi = buf[SE->_offset];
    }
    Impl::multLink(Uchi, U._odata[sU], chi, dir, SE, st);
    spReconXp(result, Uchi);
  }

  // Yp
  if (gamma == Yp) {
    if (SE->_is_local && SE->_permute) {
      spProjYp(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else if (SE->_is_local) {
      spProjYp(chi, in._odata[SE->_offset]);
    } else {
      chi = buf[SE->_offset];
    }
    Impl::multLink(Uchi, U._odata[sU], chi, dir, SE, st);
    spReconYp(result, Uchi);
  }

  // Zp
  if (gamma == Zp) {
    if (SE->_is_local && SE->_permute) {
      spProjZp(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else if (SE->_is_local) {
      spProjZp(chi, in._odata[SE->_offset]);
    } else {
      chi = buf[SE->_offset];
    }
    Impl::multLink(Uchi, U._odata[sU], chi, dir, SE, st);
    spReconZp(result, Uchi);
  }

  // Tp
  if (gamma == Tp) {
    if (SE->_is_local && SE->_permute) {
      spProjTp(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else if (SE->_is_local) {
      spProjTp(chi, in._odata[SE->_offset]);
    } else {
      chi = buf[SE->_offset];
    }
    Impl::multLink(Uchi, U._odata[sU], chi, dir, SE, st);
    spReconTp(result, Uchi);
  }

  // Xm
  if (gamma == Xm) {
    if (SE->_is_local && SE->_permute) {
      spProjXm(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else if (SE->_is_local) {
      spProjXm(chi, in._odata[SE->_offset]);
    } else {
      chi = buf[SE->_offset];
    }
    Impl::multLink(Uchi, U._odata[sU], chi, dir, SE, st);
    spReconXm(result, Uchi);
  }

  // Ym
  if (gamma == Ym) {
    if (SE->_is_local && SE->_permute) {
      spProjYm(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else if (SE->_is_local) {
      spProjYm(chi, in._odata[SE->_offset]);
    } else {
      chi = buf[SE->_offset];
    }
    Impl::multLink(Uchi, U._odata[sU], chi, dir, SE, st);
    spReconYm(result, Uchi);
  }

  // Zm
  if (gamma == Zm) {
    if (SE->_is_local && SE->_permute) {
      spProjZm(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else if (SE->_is_local) {
      spProjZm(chi, in._odata[SE->_offset]);
    } else {
      chi = buf[SE->_offset];
    }
    Impl::multLink(Uchi, U._odata[sU], chi, dir, SE, st);
    spReconZm(result, Uchi);
  }

  // Tm
  if (gamma == Tm) {
    if (SE->_is_local && SE->_permute) {
      spProjTm(tmp, in._odata[SE->_offset]);
      permute(chi, tmp, ptype);
    } else if (SE->_is_local) {
      spProjTm(chi, in._odata[SE->_offset]);
    } else {
      chi = buf[SE->_offset];
    }
    Impl::multLink(Uchi, U._odata[sU], chi, dir, SE, st);
    spReconTm(result, Uchi);
  }

  vstream(out._odata[sF], result);
}

FermOpTemplateInstantiate(WilsonKernels);
AdjointFermOpTemplateInstantiate(WilsonKernels);
TwoIndexFermOpTemplateInstantiate(WilsonKernels);

}}

