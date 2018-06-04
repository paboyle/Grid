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

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Helper routines that implement Wilson stencil for a single site.
  // Common to both the WilsonFermion and WilsonFermion5D
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class WilsonKernelsStatic { 
 public:
  enum { OptGeneric, OptHandUnroll, OptInlineAsm };
  enum { CommsAndCompute, CommsThenCompute };
  static int Opt;  
  static int Comms;
};
 
template<class Impl> class WilsonKernels : public FermionOperator<Impl> , public WilsonKernelsStatic { 
 public:
   
  INHERIT_IMPL_TYPES(Impl);
  typedef FermionOperator<Impl> Base;
   
public:

  template <bool EnableBool = true>
  typename std::enable_if<Impl::isFundamental==true && Nc == 3 &&EnableBool, void>::type
  DhopSite(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
		   int sF, int sU, int Ls, int Ns, const FermionField &in, FermionField &out,int interior=1,int exterior=1) 
  {
    bgq_l1p_optimisation(1);
    switch(Opt) {
#if defined(AVX512) || defined (QPX)
    case OptInlineAsm:
      if(interior&&exterior) WilsonKernels<Impl>::AsmDhopSite   (st,lo,U,buf,sF,sU,Ls,Ns,in,out);
      else if (interior)     WilsonKernels<Impl>::AsmDhopSiteInt(st,lo,U,buf,sF,sU,Ls,Ns,in,out);
      else if (exterior)     WilsonKernels<Impl>::AsmDhopSiteExt(st,lo,U,buf,sF,sU,Ls,Ns,in,out);
      else assert(0);
      break;
#endif
    case OptHandUnroll:
         for (int site = 0; site < Ns; site++) {
	   for (int s = 0; s < Ls; s++) {
	     if(interior&&exterior) WilsonKernels<Impl>::HandDhopSite(st,lo,U,buf,sF,sU,in,out);
	     else if (interior)     WilsonKernels<Impl>::HandDhopSiteInt(st,lo,U,buf,sF,sU,in,out);
	     else if (exterior)     WilsonKernels<Impl>::HandDhopSiteExt(st,lo,U,buf,sF,sU,in,out);
	     sF++;
	   }
	   sU++;
         }
      break;
    case OptGeneric:
         for (int site = 0; site < Ns; site++) {
	   for (int s = 0; s < Ls; s++) {
	     if(interior&&exterior) WilsonKernels<Impl>::GenericDhopSite(st,lo,U,buf,sF,sU,in,out);
	     else if (interior)     WilsonKernels<Impl>::GenericDhopSiteInt(st,lo,U,buf,sF,sU,in,out);
	     else if (exterior)     WilsonKernels<Impl>::GenericDhopSiteExt(st,lo,U,buf,sF,sU,in,out);
	     else assert(0);
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
  typename std::enable_if<(Impl::isFundamental==false || (Impl::isFundamental==true && Nc != 3)) && EnableBool, void>::type
  DhopSite(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
	   int sF, int sU, int Ls, int Ns, const FermionField &in, FermionField &out,int interior=1,int exterior=1 ) {
    // no kernel choice  
    for (int site = 0; site < Ns; site++) {
      for (int s = 0; s < Ls; s++) {
	if(interior&&exterior) WilsonKernels<Impl>::GenericDhopSite(st,lo,U,buf,sF,sU,in,out);
	else if (interior)     WilsonKernels<Impl>::GenericDhopSiteInt(st,lo,U,buf,sF,sU,in,out);
	else if (exterior)     WilsonKernels<Impl>::GenericDhopSiteExt(st,lo,U,buf,sF,sU,in,out);
	else assert(0);
	sF++;
      }
      sU++;
    }
  }
     
  template <bool EnableBool = true>
  typename std::enable_if<Impl::isFundamental==true && Nc == 3 && EnableBool,void>::type
  DhopSiteDag(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
	      int sF, int sU, int Ls, int Ns, const FermionField &in, FermionField &out,int interior=1,int exterior=1) 
{
    bgq_l1p_optimisation(1);
    switch(Opt) {
#if defined(AVX512) || defined (QPX)
    case OptInlineAsm:
      if(interior&&exterior) WilsonKernels<Impl>::AsmDhopSiteDag   (st,lo,U,buf,sF,sU,Ls,Ns,in,out);
      else if (interior)     WilsonKernels<Impl>::AsmDhopSiteDagInt(st,lo,U,buf,sF,sU,Ls,Ns,in,out);
      else if (exterior)     WilsonKernels<Impl>::AsmDhopSiteDagExt(st,lo,U,buf,sF,sU,Ls,Ns,in,out);
      else assert(0);
      break;
#endif
    case OptHandUnroll:
      for (int site = 0; site < Ns; site++) {
	for (int s = 0; s < Ls; s++) {
	  if(interior&&exterior) WilsonKernels<Impl>::HandDhopSiteDag(st,lo,U,buf,sF,sU,in,out);
	  else if (interior)     WilsonKernels<Impl>::HandDhopSiteDagInt(st,lo,U,buf,sF,sU,in,out);
	  else if (exterior)     WilsonKernels<Impl>::HandDhopSiteDagExt(st,lo,U,buf,sF,sU,in,out);
	  else assert(0);
	  sF++;
	}
	sU++;
      }
      break;
    case OptGeneric:
      for (int site = 0; site < Ns; site++) {
	for (int s = 0; s < Ls; s++) {
	  if(interior&&exterior) WilsonKernels<Impl>::GenericDhopSiteDag(st,lo,U,buf,sF,sU,in,out);
	  else if (interior)     WilsonKernels<Impl>::GenericDhopSiteDagInt(st,lo,U,buf,sF,sU,in,out);
	  else if (exterior)     WilsonKernels<Impl>::GenericDhopSiteDagExt(st,lo,U,buf,sF,sU,in,out);
	  else assert(0);
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
  typename std::enable_if<(Impl::isFundamental==false || (Impl::isFundamental==true && Nc != 3)) && EnableBool,void>::type
  DhopSiteDag(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U,SiteHalfSpinor * buf,
		      int sF, int sU, int Ls, int Ns, const FermionField &in, FermionField &out,int interior=1,int exterior=1) {

    for (int site = 0; site < Ns; site++) {
      for (int s = 0; s < Ls; s++) {
	if(interior&&exterior) WilsonKernels<Impl>::GenericDhopSiteDag(st,lo,U,buf,sF,sU,in,out);
	else if (interior)     WilsonKernels<Impl>::GenericDhopSiteDagInt(st,lo,U,buf,sF,sU,in,out);
	else if (exterior)     WilsonKernels<Impl>::GenericDhopSiteDagExt(st,lo,U,buf,sF,sU,in,out);
	else assert(0);
	sF++;
      }
      sU++;
    }
  }

  void DhopDir(StencilImpl &st, DoubledGaugeField &U,SiteHalfSpinor * buf,
		       int sF, int sU, const FermionField &in, FermionField &out, int dirdisp, int gamma);
      
  //////////////////////////////////////////////////////////////////////////////
  // Utilities for inserting Wilson conserved current.
  //////////////////////////////////////////////////////////////////////////////
  void ContractConservedCurrentSiteFwd(const SitePropagator &q_in_1,
                                       const SitePropagator &q_in_2,
                                       SitePropagator &q_out,
                                       DoubledGaugeField &U,
                                       unsigned int sU,
                                       unsigned int mu,
                                       bool switch_sign = false);
  void ContractConservedCurrentSiteBwd(const SitePropagator &q_in_1,
                                       const SitePropagator &q_in_2,
                                       SitePropagator &q_out,
                                       DoubledGaugeField &U,
                                       unsigned int sU,
                                       unsigned int mu,
                                       bool switch_sign = false);
  void SeqConservedCurrentSiteFwd(const SitePropagator &q_in, 
                                  SitePropagator &q_out,
                                  DoubledGaugeField &U,
                                  unsigned int sU,
                                  unsigned int mu,
                                  vInteger t_mask,
                                  bool switch_sign = false);
  void SeqConservedCurrentSiteBwd(const SitePropagator &q_in,
                                  SitePropagator &q_out,
                                  DoubledGaugeField &U,
                                  unsigned int sU,
                                  unsigned int mu,
                                  vInteger t_mask,
                                  bool switch_sign = false);

private:
     // Specialised variants
  void GenericDhopSite(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
		       int sF, int sU, const FermionField &in, FermionField &out);
      
  void GenericDhopSiteDag(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
			  int sF, int sU, const FermionField &in, FermionField &out);

  void GenericDhopSiteInt(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
			  int sF, int sU, const FermionField &in, FermionField &out);
      
  void GenericDhopSiteDagInt(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
			     int sF, int sU, const FermionField &in, FermionField &out);

  void GenericDhopSiteExt(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
			  int sF, int sU, const FermionField &in, FermionField &out);
      
  void GenericDhopSiteDagExt(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
			     int sF, int sU, const FermionField &in, FermionField &out);


  void AsmDhopSite(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
		   int sF, int sU, int Ls, int Ns, const FermionField &in,FermionField &out);

  void AsmDhopSiteDag(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
		      int sF, int sU, int Ls, int Ns, const FermionField &in, FermionField &out);

  void AsmDhopSiteInt(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
		      int sF, int sU, int Ls, int Ns, const FermionField &in,FermionField &out);

  void AsmDhopSiteDagInt(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
			 int sF, int sU, int Ls, int Ns, const FermionField &in, FermionField &out);

  void AsmDhopSiteExt(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
		      int sF, int sU, int Ls, int Ns, const FermionField &in,FermionField &out);

  void AsmDhopSiteDagExt(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
			 int sF, int sU, int Ls, int Ns, const FermionField &in, FermionField &out);


  void HandDhopSite(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
		    int sF, int sU, const FermionField &in, FermionField &out);

  void HandDhopSiteDag(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
		       int sF, int sU, const FermionField &in, FermionField &out);
      
  void HandDhopSiteInt(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
		       int sF, int sU, const FermionField &in, FermionField &out);
  
  void HandDhopSiteDagInt(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
			  int sF, int sU, const FermionField &in, FermionField &out);
  
  void HandDhopSiteExt(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
		       int sF, int sU, const FermionField &in, FermionField &out);
  
  void HandDhopSiteDagExt(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteHalfSpinor * buf,
			  int sF, int sU, const FermionField &in, FermionField &out);
  
public:

  WilsonKernels(const ImplParams &p = ImplParams());

};
    
}}

#endif
