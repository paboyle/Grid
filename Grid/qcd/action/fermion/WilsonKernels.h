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
#pragma once

NAMESPACE_BEGIN(Grid);

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

  static void DhopKernel(int Opt,StencilImpl &st,  DoubledGaugeField &U, SiteHalfSpinor * buf,
			 int Ls, int Nsite, const FermionField &in, FermionField &out,
			 int interior=1,int exterior=1) ;

  static void DhopDagKernel(int Opt,StencilImpl &st,  DoubledGaugeField &U, SiteHalfSpinor * buf,
			    int Ls, int Nsite, const FermionField &in, FermionField &out,
			    int interior=1,int exterior=1) ;

  static void DhopDirKernel(StencilImpl &st, DoubledGaugeField &U,SiteHalfSpinor * buf,
			    int Ls, int Nsite, const FermionField &in, FermionField &out, int dirdisp, int gamma);

  //////////////////////////////////////////////////////////////////////////////
  // Utilities for inserting Wilson conserved current.
  //////////////////////////////////////////////////////////////////////////////
  static void ContractConservedCurrentSiteFwd(const SitePropagator &q_in_1,
                                       const SitePropagator &q_in_2,
                                       SitePropagator &q_out,
                                       DoubledGaugeFieldView &U,
                                       unsigned int sU,
                                       unsigned int mu,
                                       bool switch_sign = false);

  static void ContractConservedCurrentSiteBwd(const SitePropagator &q_in_1,
                                       const SitePropagator &q_in_2,
                                       SitePropagator &q_out,
                                       DoubledGaugeFieldView &U,
                                       unsigned int sU,
                                       unsigned int mu,
                                       bool switch_sign = false);

  static void SeqConservedCurrentSiteFwd(const SitePropagator &q_in, 
                                  SitePropagator &q_out,
                                  DoubledGaugeFieldView &U,
                                  unsigned int sU,
                                  unsigned int mu,
                                  vPredicate t_mask,
                                  bool switch_sign = false);

  static void SeqConservedCurrentSiteBwd(const SitePropagator &q_in,
                                  SitePropagator &q_out,
                                  DoubledGaugeFieldView &U,
                                  unsigned int sU,
                                  unsigned int mu,
                                  vPredicate t_mask,
                                  bool switch_sign = false);

private:

  static accelerator void DhopDirK(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor * buf,
				   int sF, int sU, const FermionFieldView &in, FermionFieldView &out, int dirdisp, int gamma);
      
  // Specialised variants
  static accelerator void GenericDhopSite(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
					  int sF, int sU, const FermionFieldView &in, FermionFieldView &out);
      
  static accelerator void GenericDhopSiteDag(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
						    int sF, int sU, const FermionFieldView &in, FermionFieldView &out);
  
  static accelerator void GenericDhopSiteInt(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
						    int sF, int sU, const FermionFieldView &in, FermionFieldView &out);
      
  static accelerator void GenericDhopSiteDagInt(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
						int sF, int sU, const FermionFieldView &in, FermionFieldView &out);
  
  static accelerator void GenericDhopSiteExt(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
					     int sF, int sU, const FermionFieldView &in, FermionFieldView &out);
      
  static accelerator void GenericDhopSiteDagExt(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
						       int sF, int sU, const FermionFieldView &in, FermionFieldView &out);

  static void AsmDhopSite(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
			  int sF, int sU, int Ls, int Nsite, const FermionFieldView &in,FermionFieldView &out);
  
  static void AsmDhopSiteDag(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
			     int sF, int sU, int Ls, int Nsite, const FermionFieldView &in, FermionFieldView &out);
  
  static void AsmDhopSiteInt(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
			     int sF, int sU, int Ls, int Nsite, const FermionFieldView &in,FermionFieldView &out);
  
  static void AsmDhopSiteDagInt(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
				int sF, int sU, int Ls, int Nsite, const FermionFieldView &in, FermionFieldView &out);
  
  static void AsmDhopSiteExt(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
			     int sF, int sU, int Ls, int Nsite, const FermionFieldView &in,FermionFieldView &out);
  
  static void AsmDhopSiteDagExt(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
				int sF, int sU, int Ls, int Nsite, const FermionFieldView &in, FermionFieldView &out);

// Keep Hand unrolled temporarily  
  static accelerator void HandDhopSite(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
				       int sF, int sU, const FermionFieldView &in, FermionFieldView &out);
  
  static accelerator void HandDhopSiteDag(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
					  int sF, int sU, const FermionFieldView &in, FermionFieldView &out);
  
  static accelerator void HandDhopSiteInt(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
					  int sF, int sU, const FermionFieldView &in, FermionFieldView &out);
  
  static accelerator void HandDhopSiteDagInt(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
					     int sF, int sU, const FermionFieldView &in, FermionFieldView &out);
  
  static accelerator void HandDhopSiteExt(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
					  int sF, int sU, const FermionFieldView &in, FermionFieldView &out);
  
  static accelerator void HandDhopSiteDagExt(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
					     int sF, int sU, const FermionFieldView &in, FermionFieldView &out);
 public:
 WilsonKernels(const ImplParams &p = ImplParams()) : Base(p){};
};
    
NAMESPACE_END(Grid);


