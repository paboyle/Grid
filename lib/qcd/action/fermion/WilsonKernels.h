/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/WilsonKernels.h

Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
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
#ifndef GRID_QCD_DHOP_H
#define GRID_QCD_DHOP_H

namespace Grid {
namespace QCD {

void bgq_l1p_optimisation(int mode);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Helper routines that implement Wilson stencil for a single site.
  // Common to both the WilsonFermion and WilsonFermion5D
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class WilsonKernelsStatic { 
 public:
  enum { OptGeneric, OptHandUnroll, OptInlineAsm };
  // S-direction is INNERMOST and takes no part in the parity.
  static int Opt;  // these are a temporary hack
};
 
template<class Impl> class WilsonKernels : public FermionOperator<Impl> , public WilsonKernelsStatic { 
 public:
   
  INHERIT_IMPL_TYPES(Impl);
  typedef FermionOperator<Impl> Base;
   
public:
   
  template <bool EnableBool = true>
  typename std::enable_if<Impl::Dimension == 3 && Nc == 3 &&EnableBool, void>::type
  DiracOptDhopSite(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
		   int sF, int sU, int Ls, int Ns, const FermionField &in, FermionField &out) 
  {
    bgq_l1p_optimisation(1);
    switch(Opt) {
#if defined(AVX512) || defined (QPX)
    case OptInlineAsm:
       WilsonKernels<Impl>::DiracOptAsmDhopSite(st,lo,U,buf,sF,sU,Ls,Ns,in,out);
       break;
#endif
    case OptHandUnroll:
      for (int site = 0; site < Ns; site++) {
	for (int s = 0; s < Ls; s++) {
	  WilsonKernels<Impl>::DiracOptHandDhopSite(st,lo,U,buf,sF,sU,in,out);
	  sF++;
	}
	sU++;
      }
      break;
    case OptGeneric:
      for (int site = 0; site < Ns; site++) {
	for (int s = 0; s < Ls; s++) {
	  WilsonKernels<Impl>::DiracOptGenericDhopSite(st,lo,U,buf,sF,sU,in,out);
	  sF++;
	}
	sU++;
      }
      break;
    default:
      assert(0);
    }
    bgq_l1p_optimisation(0);
  }
     
  template <bool EnableBool = true>
  typename std::enable_if<(Impl::Dimension != 3 || (Impl::Dimension == 3 && Nc != 3)) && EnableBool, void>::type
  DiracOptDhopSite(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
		   int sF, int sU, int Ls, int Ns, const FermionField &in, FermionField &out) {
    // no kernel choice  
    for (int site = 0; site < Ns; site++) {
      for (int s = 0; s < Ls; s++) {
	WilsonKernels<Impl>::DiracOptGenericDhopSite(st, lo, U, buf, sF, sU, in, out);
	sF++;
      }
      sU++;
    }
  }
     
  template <bool EnableBool = true>
  typename std::enable_if<Impl::Dimension == 3 && Nc == 3 && EnableBool,void>::type
  DiracOptDhopSiteDag(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
		      int sF, int sU, int Ls, int Ns, const FermionField &in, FermionField &out) {

    bgq_l1p_optimisation(1);
    switch(Opt) {
#if defined(AVX512) || defined (QPX)
    case OptInlineAsm:
      WilsonKernels<Impl>::DiracOptAsmDhopSiteDag(st,lo,U,buf,sF,sU,Ls,Ns,in,out);
      break;
#endif
    case OptHandUnroll:
      for (int site = 0; site < Ns; site++) {
	for (int s = 0; s < Ls; s++) {
	  WilsonKernels<Impl>::DiracOptHandDhopSiteDag(st,lo,U,buf,sF,sU,in,out);
	  sF++;
	}
	sU++;
      }
      break;
    case OptGeneric:
      for (int site = 0; site < Ns; site++) {
	for (int s = 0; s < Ls; s++) {
	  WilsonKernels<Impl>::DiracOptGenericDhopSiteDag(st,lo,U,buf,sF,sU,in,out);
	  sF++;
	}
	sU++;
      }
      break;
    default:
      assert(0);
    }
    bgq_l1p_optimisation(0);
  }

  template <bool EnableBool = true>
  typename std::enable_if<(Impl::Dimension != 3 || (Impl::Dimension == 3 && Nc != 3)) && EnableBool,void>::type
  DiracOptDhopSiteDag(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U,SiteHalfSpinor * buf,
		      int sF, int sU, int Ls, int Ns, const FermionField &in, FermionField &out) {

    for (int site = 0; site < Ns; site++) {
      for (int s = 0; s < Ls; s++) {
	WilsonKernels<Impl>::DiracOptGenericDhopSiteDag(st,lo,U,buf,sF,sU,in,out);
	sF++;
      }
      sU++;
    }
  }

  void DiracOptDhopDir(StencilImpl &st, DoubledGaugeField &U,SiteHalfSpinor * buf,
		       int sF, int sU, const FermionField &in, FermionField &out, int dirdisp, int gamma);
      
private:
     // Specialised variants
  void DiracOptGenericDhopSite(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
			       int sF, int sU, const FermionField &in, FermionField &out);
      
  void DiracOptGenericDhopSiteDag(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
				  int sF, int sU, const FermionField &in, FermionField &out);

  void DiracOptAsmDhopSite(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
			   int sF, int sU, int Ls, int Ns, const FermionField &in,FermionField &out);

  void DiracOptAsmDhopSiteDag(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
			      int sF, int sU, int Ls, int Ns, const FermionField &in, FermionField &out);

  void DiracOptHandDhopSite(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
			    int sF, int sU, const FermionField &in, FermionField &out);

  void DiracOptHandDhopSiteDag(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
			       int sF, int sU, const FermionField &in, FermionField &out);
      
public:

  WilsonKernels(const ImplParams &p = ImplParams());

};
    
}}

#endif
