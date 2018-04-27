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
#include <Grid/qcd/action/fermion/FermionCore.h>

namespace Grid {
namespace QCD {

int WilsonKernelsStatic::Opt   = WilsonKernelsStatic::OptGeneric;
int WilsonKernelsStatic::Comms = WilsonKernelsStatic::CommsAndCompute;

template <class Impl>
WilsonKernels<Impl>::WilsonKernels(const ImplParams &p) : Base(p){};

////////////////////////////////////////////
// Generic implementation; move to different file?
////////////////////////////////////////////
  
#define GENERIC_STENCIL_LEG(Dir,spProj,Recon)			\
  SE = st.GetEntry(ptype, Dir, sF);				\
  if (SE->_is_local) {						\
    chi_p = &chi;						\
    if (SE->_permute) {						\
      spProj(tmp, in._odata[SE->_offset]);			\
      permute(chi, tmp, ptype);					\
    } else {							\
      spProj(chi, in._odata[SE->_offset]);			\
    }								\
  } else {							\
    chi_p = &buf[SE->_offset];					\
  }								\
  Impl::multLink(Uchi, U._odata[sU], *chi_p, Dir, SE, st);	\
  Recon(result, Uchi);
  
#define GENERIC_STENCIL_LEG_INT(Dir,spProj,Recon)		\
  SE = st.GetEntry(ptype, Dir, sF);				\
  if (SE->_is_local) {						\
    chi_p = &chi;						\
    if (SE->_permute) {						\
      spProj(tmp, in._odata[SE->_offset]);			\
      permute(chi, tmp, ptype);					\
    } else {							\
      spProj(chi, in._odata[SE->_offset]);			\
    }								\
  } else if ( st.same_node[Dir] ) {				\
      chi_p = &buf[SE->_offset];				\
  }								\
  if (SE->_is_local || st.same_node[Dir] ) {			\
    Impl::multLink(Uchi, U._odata[sU], *chi_p, Dir, SE, st);	\
    Recon(result, Uchi);					\
  }

#define GENERIC_STENCIL_LEG_EXT(Dir,spProj,Recon)		\
  SE = st.GetEntry(ptype, Dir, sF);				\
  if ((!SE->_is_local) && (!st.same_node[Dir]) ) {		\
    chi_p = &buf[SE->_offset];					\
    Impl::multLink(Uchi, U._odata[sU], *chi_p, Dir, SE, st);	\
    Recon(result, Uchi);					\
    nmu++;							\
  }

#define GENERIC_DHOPDIR_LEG(Dir,spProj,Recon)			\
  if (gamma == Dir) {						\
    if (SE->_is_local && SE->_permute) {			\
      spProj(tmp, in._odata[SE->_offset]);			\
      permute(chi, tmp, ptype);					\
    } else if (SE->_is_local) {					\
      spProj(chi, in._odata[SE->_offset]);			\
    } else {							\
      chi = buf[SE->_offset];					\
    }								\
    Impl::multLink(Uchi, U._odata[sU], chi, dir, SE, st);	\
    Recon(result, Uchi);					\
  }

  ////////////////////////////////////////////////////////////////////
  // All legs kernels ; comms then compute
  ////////////////////////////////////////////////////////////////////
template <class Impl>
void WilsonKernels<Impl>::GenericDhopSiteDag(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U,
					     SiteHalfSpinor *buf, int sF,
					     int sU, const FermionField &in, FermionField &out)
{
  SiteHalfSpinor tmp;
  SiteHalfSpinor chi;
  SiteHalfSpinor *chi_p;
  SiteHalfSpinor Uchi;
  SiteSpinor result;
  StencilEntry *SE;
  int ptype;

  GENERIC_STENCIL_LEG(Xp,spProjXp,spReconXp);
  GENERIC_STENCIL_LEG(Yp,spProjYp,accumReconYp);
  GENERIC_STENCIL_LEG(Zp,spProjZp,accumReconZp);
  GENERIC_STENCIL_LEG(Tp,spProjTp,accumReconTp);
  GENERIC_STENCIL_LEG(Xm,spProjXm,accumReconXm);
  GENERIC_STENCIL_LEG(Ym,spProjYm,accumReconYm);
  GENERIC_STENCIL_LEG(Zm,spProjZm,accumReconZm);
  GENERIC_STENCIL_LEG(Tm,spProjTm,accumReconTm);
  vstream(out._odata[sF], result);
};

template <class Impl>
void WilsonKernels<Impl>::GenericDhopSite(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U,
					  SiteHalfSpinor *buf, int sF,
					  int sU, const FermionField &in, FermionField &out) 
{
  SiteHalfSpinor tmp;
  SiteHalfSpinor chi;
  SiteHalfSpinor *chi_p;
  SiteHalfSpinor Uchi;
  SiteSpinor result;
  StencilEntry *SE;
  int ptype;

  GENERIC_STENCIL_LEG(Xm,spProjXp,spReconXp);
  GENERIC_STENCIL_LEG(Ym,spProjYp,accumReconYp);
  GENERIC_STENCIL_LEG(Zm,spProjZp,accumReconZp);
  GENERIC_STENCIL_LEG(Tm,spProjTp,accumReconTp);
  GENERIC_STENCIL_LEG(Xp,spProjXm,accumReconXm);
  GENERIC_STENCIL_LEG(Yp,spProjYm,accumReconYm);
  GENERIC_STENCIL_LEG(Zp,spProjZm,accumReconZm);
  GENERIC_STENCIL_LEG(Tp,spProjTm,accumReconTm);
  vstream(out._odata[sF], result);
};
  ////////////////////////////////////////////////////////////////////
  // Interior kernels
  ////////////////////////////////////////////////////////////////////
template <class Impl>
void WilsonKernels<Impl>::GenericDhopSiteDagInt(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U,
						SiteHalfSpinor *buf, int sF,
						int sU, const FermionField &in, FermionField &out)
{
  SiteHalfSpinor tmp;
  SiteHalfSpinor chi;
  SiteHalfSpinor *chi_p;
  SiteHalfSpinor Uchi;
  SiteSpinor result;
  StencilEntry *SE;
  int ptype;

  result=zero;
  GENERIC_STENCIL_LEG_INT(Xp,spProjXp,accumReconXp);
  GENERIC_STENCIL_LEG_INT(Yp,spProjYp,accumReconYp);
  GENERIC_STENCIL_LEG_INT(Zp,spProjZp,accumReconZp);
  GENERIC_STENCIL_LEG_INT(Tp,spProjTp,accumReconTp);
  GENERIC_STENCIL_LEG_INT(Xm,spProjXm,accumReconXm);
  GENERIC_STENCIL_LEG_INT(Ym,spProjYm,accumReconYm);
  GENERIC_STENCIL_LEG_INT(Zm,spProjZm,accumReconZm);
  GENERIC_STENCIL_LEG_INT(Tm,spProjTm,accumReconTm);
  vstream(out._odata[sF], result);
};

template <class Impl>
void WilsonKernels<Impl>::GenericDhopSiteInt(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U,
					     SiteHalfSpinor *buf, int sF,
					     int sU, const FermionField &in, FermionField &out) 
{
  SiteHalfSpinor tmp;
  SiteHalfSpinor chi;
  SiteHalfSpinor *chi_p;
  SiteHalfSpinor Uchi;
  SiteSpinor result;
  StencilEntry *SE;
  int ptype;
  result=zero;
  GENERIC_STENCIL_LEG_INT(Xm,spProjXp,accumReconXp);
  GENERIC_STENCIL_LEG_INT(Ym,spProjYp,accumReconYp);
  GENERIC_STENCIL_LEG_INT(Zm,spProjZp,accumReconZp);
  GENERIC_STENCIL_LEG_INT(Tm,spProjTp,accumReconTp);
  GENERIC_STENCIL_LEG_INT(Xp,spProjXm,accumReconXm);
  GENERIC_STENCIL_LEG_INT(Yp,spProjYm,accumReconYm);
  GENERIC_STENCIL_LEG_INT(Zp,spProjZm,accumReconZm);
  GENERIC_STENCIL_LEG_INT(Tp,spProjTm,accumReconTm);
  vstream(out._odata[sF], result);
};
////////////////////////////////////////////////////////////////////
// Exterior kernels
////////////////////////////////////////////////////////////////////
template <class Impl>
void WilsonKernels<Impl>::GenericDhopSiteDagExt(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U,
						SiteHalfSpinor *buf, int sF,
						int sU, const FermionField &in, FermionField &out)
{
  SiteHalfSpinor tmp;
  SiteHalfSpinor chi;
  SiteHalfSpinor *chi_p;
  SiteHalfSpinor Uchi;
  SiteSpinor result;
  StencilEntry *SE;
  int ptype;
  int nmu=0;
  result=zero;
  GENERIC_STENCIL_LEG_EXT(Xp,spProjXp,accumReconXp);
  GENERIC_STENCIL_LEG_EXT(Yp,spProjYp,accumReconYp);
  GENERIC_STENCIL_LEG_EXT(Zp,spProjZp,accumReconZp);
  GENERIC_STENCIL_LEG_EXT(Tp,spProjTp,accumReconTp);
  GENERIC_STENCIL_LEG_EXT(Xm,spProjXm,accumReconXm);
  GENERIC_STENCIL_LEG_EXT(Ym,spProjYm,accumReconYm);
  GENERIC_STENCIL_LEG_EXT(Zm,spProjZm,accumReconZm);
  GENERIC_STENCIL_LEG_EXT(Tm,spProjTm,accumReconTm);
  if ( nmu ) { 
    out._odata[sF] = out._odata[sF] + result; 
  }
};

template <class Impl>
void WilsonKernels<Impl>::GenericDhopSiteExt(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U,
					     SiteHalfSpinor *buf, int sF,
					     int sU, const FermionField &in, FermionField &out) 
{
  SiteHalfSpinor tmp;
  SiteHalfSpinor chi;
  SiteHalfSpinor *chi_p;
  SiteHalfSpinor Uchi;
  SiteSpinor result;
  StencilEntry *SE;
  int ptype;
  int nmu=0;
  result=zero;
  GENERIC_STENCIL_LEG_EXT(Xm,spProjXp,accumReconXp);
  GENERIC_STENCIL_LEG_EXT(Ym,spProjYp,accumReconYp);
  GENERIC_STENCIL_LEG_EXT(Zm,spProjZp,accumReconZp);
  GENERIC_STENCIL_LEG_EXT(Tm,spProjTp,accumReconTp);
  GENERIC_STENCIL_LEG_EXT(Xp,spProjXm,accumReconXm);
  GENERIC_STENCIL_LEG_EXT(Yp,spProjYm,accumReconYm);
  GENERIC_STENCIL_LEG_EXT(Zp,spProjZm,accumReconZm);
  GENERIC_STENCIL_LEG_EXT(Tp,spProjTm,accumReconTm);
  if ( nmu ) { 
    out._odata[sF] = out._odata[sF] + result; 
  }
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
  GENERIC_DHOPDIR_LEG(Xp,spProjXp,spReconXp);
  GENERIC_DHOPDIR_LEG(Yp,spProjYp,spReconYp);
  GENERIC_DHOPDIR_LEG(Zp,spProjZp,spReconZp);
  GENERIC_DHOPDIR_LEG(Tp,spProjTp,spReconTp);
  GENERIC_DHOPDIR_LEG(Xm,spProjXm,spReconXm);
  GENERIC_DHOPDIR_LEG(Ym,spProjYm,spReconYm);
  GENERIC_DHOPDIR_LEG(Zm,spProjZm,spReconZm);
  GENERIC_DHOPDIR_LEG(Tm,spProjTm,spReconTm);
  vstream(out._odata[sF], result);
}

/*******************************************************************************
 * Conserved current utilities for Wilson fermions, for contracting propagators
 * to make a conserved current sink or inserting the conserved current 
 * sequentially. Common to both 4D and 5D.
 ******************************************************************************/
// N.B. Functions below assume a -1/2 factor within U.
#define WilsonCurrentFwd(expr, mu) ((expr - Gamma::gmu[mu]*expr))
#define WilsonCurrentBwd(expr, mu) ((expr + Gamma::gmu[mu]*expr))

/*******************************************************************************
 * Name: ContractConservedCurrentSiteFwd
 * Operation: (1/2) * q2[x] * U(x) * (g[mu] - 1) * q1[x + mu]
 * Notes: - DoubledGaugeField U assumed to contain -1/2 factor.
 *        - Pass in q_in_1 shifted in +ve mu direction.
 ******************************************************************************/
template<class Impl>
void WilsonKernels<Impl>::ContractConservedCurrentSiteFwd(
                                                  const SitePropagator &q_in_1,
                                                  const SitePropagator &q_in_2,
                                                  SitePropagator &q_out,
                                                  DoubledGaugeField &U,
                                                  unsigned int sU,
                                                  unsigned int mu,
                                                  bool switch_sign)
{
    SitePropagator result, tmp;
    Gamma g5(Gamma::Algebra::Gamma5);
    Impl::multLinkProp(tmp, U._odata[sU], q_in_1, mu);
    result = g5 * adj(q_in_2) * g5 * WilsonCurrentFwd(tmp, mu);
    if (switch_sign)
    {
        q_out -= result;
    }
    else
    {
        q_out += result;
    }
}

/*******************************************************************************
 * Name: ContractConservedCurrentSiteBwd
 * Operation: (1/2) * q2[x + mu] * U^dag(x) * (g[mu] + 1) * q1[x]
 * Notes: - DoubledGaugeField U assumed to contain -1/2 factor.
 *        - Pass in q_in_2 shifted in +ve mu direction.
 ******************************************************************************/
template<class Impl>
void WilsonKernels<Impl>::ContractConservedCurrentSiteBwd(
                                                  const SitePropagator &q_in_1,
                                                  const SitePropagator &q_in_2,
                                                  SitePropagator &q_out,
                                                  DoubledGaugeField &U,
                                                  unsigned int sU,
                                                  unsigned int mu,
                                                  bool switch_sign)
{
    SitePropagator result, tmp;
    Gamma g5(Gamma::Algebra::Gamma5);
    Impl::multLinkProp(tmp, U._odata[sU], q_in_1, mu + Nd);
    result = g5 * adj(q_in_2) * g5 * WilsonCurrentBwd(tmp, mu);
    if (switch_sign)
    {
        q_out += result;
    }
    else
    {
        q_out -= result;
    }
}

// G-parity requires more specialised implementation.
#define NO_CURR_SITE(Impl) \
template <> \
void WilsonKernels<Impl>::ContractConservedCurrentSiteFwd( \
                                                  const SitePropagator &q_in_1, \
                                                  const SitePropagator &q_in_2, \
                                                  SitePropagator &q_out,        \
                                                  DoubledGaugeField &U,         \
                                                  unsigned int sU,              \
                                                  unsigned int mu,              \
                                                  bool switch_sign)             \
{ \
    assert(0); \
} \
template <> \
void WilsonKernels<Impl>::ContractConservedCurrentSiteBwd( \
                                                  const SitePropagator &q_in_1, \
                                                  const SitePropagator &q_in_2, \
                                                  SitePropagator &q_out,        \
                                                  DoubledGaugeField &U,         \
                                                  unsigned int mu,              \
                                                  unsigned int sU,              \
                                                  bool switch_sign)             \
{ \
    assert(0); \
}

NO_CURR_SITE(GparityWilsonImplF);
NO_CURR_SITE(GparityWilsonImplD);
NO_CURR_SITE(GparityWilsonImplFH);
NO_CURR_SITE(GparityWilsonImplDF);


/*******************************************************************************
 * Name: SeqConservedCurrentSiteFwd
 * Operation: (1/2) * U(x) * (g[mu] - 1) * q[x + mu]
 * Notes: - DoubledGaugeField U assumed to contain -1/2 factor.
 *        - Pass in q_in shifted in +ve mu direction.
 ******************************************************************************/
template<class Impl>
void WilsonKernels<Impl>::SeqConservedCurrentSiteFwd(const SitePropagator &q_in,
                                                     SitePropagator &q_out,
                                                     DoubledGaugeField &U,
                                                     unsigned int sU,
                                                     unsigned int mu,
                                                     vInteger t_mask,
                                                     bool switch_sign)
{
    SitePropagator result;
    Impl::multLinkProp(result, U._odata[sU], q_in, mu);
    result = WilsonCurrentFwd(result, mu);

    // Zero any unwanted timeslice entries.
    result = predicatedWhere(t_mask, result, 0.*result);

    if (switch_sign)
    {
        q_out -= result;
    }
    else
    {
        q_out += result;
    }
}

/*******************************************************************************
 * Name: SeqConservedCurrentSiteFwd
 * Operation: (1/2) * U^dag(x) * (g[mu] + 1) * q[x - mu]
 * Notes: - DoubledGaugeField U assumed to contain -1/2 factor.
 *        - Pass in q_in shifted in -ve mu direction.
 ******************************************************************************/
template<class Impl>
void WilsonKernels<Impl>::SeqConservedCurrentSiteBwd(const SitePropagator &q_in, 
                                                     SitePropagator &q_out,
                                                     DoubledGaugeField &U,
                                                     unsigned int sU,
                                                     unsigned int mu,
                                                     vInteger t_mask,
                                                     bool switch_sign)
{
    SitePropagator result;
    Impl::multLinkProp(result, U._odata[sU], q_in, mu + Nd);
    result = WilsonCurrentBwd(result, mu);

    // Zero any unwanted timeslice entries.
    result = predicatedWhere(t_mask, result, 0.*result);

    if (switch_sign)
    {
        q_out += result;
    }
    else
    {
        q_out -= result;
    }
}

FermOpTemplateInstantiate(WilsonKernels);
AdjointFermOpTemplateInstantiate(WilsonKernels);
TwoIndexFermOpTemplateInstantiate(WilsonKernels);

}}

