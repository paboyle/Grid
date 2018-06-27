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

NAMESPACE_BEGIN(Grid);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Helper routines that implement Wilson stencil for a single site.
// Common to both the WilsonFermion and WilsonFermion5D
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class WilsonKernelsStatic { 
public:
  enum { OptGeneric, OptHandUnroll, OptInlineAsm, OptGpu };
  enum { CommsAndCompute, CommsThenCompute };
  static int Opt;  
  static int Comms;
};
 
template<class Impl> class WilsonKernels { 
public:

  INHERIT_IMPL_TYPES(Impl);
   
public:

  static void Dhop(int Opt,StencilImpl &st,  DoubledGaugeField &U, SiteHalfSpinor * buf,
		   int Ls, int Nsite, const FermionField &in, FermionField &out,
		   int interior=1,int exterior=1) 
  {
    auto U_v   = U.View();
    auto in_v  = in.View();
    auto out_v = out.View();
    auto st_v  = st.View();
    if ( (Opt == WilsonKernelsStatic::OptGpu) && interior && exterior ) { 
      const uint64_t nsimd = Simd::Nsimd();
      const uint64_t    NN = Nsite*Ls*Simd::Nsimd();
      accelerator_loopN( sss, NN, {
	  uint64_t cur  = sss;
	  /*	  uint64_t lane = cur % nsimd;  */ cur = cur / nsimd;
	  uint64_t   sF = cur;         cur = cur / Ls;
	  uint64_t   sU = cur;
	  WilsonKernels<Impl>::GpuDhopSite(st_v,U_v,buf,sF,sU,in_v,out_v);
      });
    } else { 
      accelerator_loop( ss, U_v, {
	int sU = ss;
        int sF = Ls * sU;
        DhopSite(Opt,st_v,U_v,st.CommBuf(),sF,sU,Ls,1,in_v,out_v);
      });
    }
  }
  static void DhopDag(int Opt,StencilImpl &st,  DoubledGaugeField &U, SiteHalfSpinor * buf,
		      int Ls, int Nsite, const FermionField &in, FermionField &out,
		      int interior=1,int exterior=1) 
  {
    auto U_v   = U.View();
    auto in_v  = in.View();
    auto out_v = out.View();
    auto st_v  = st.View();

    if ( (Opt == WilsonKernelsStatic::OptGpu) && interior && exterior ) { 
      const uint64_t nsimd = Simd::Nsimd();
      const uint64_t    NN = Nsite*Ls*Simd::Nsimd();
      accelerator_loopN( sss, NN, {
	  uint64_t cur  = sss;
	  /* uint64_t lane = cur % nsimd; */ cur = cur / nsimd;
	  uint64_t   sF = cur;         cur = cur / Ls;
	  uint64_t   sU = cur;
	  WilsonKernels<Impl>::GpuDhopSiteDag(st_v,U_v,buf,sF,sU,in_v,out_v);
      });
    } else { 
      accelerator_loop( ss, U_v, {
	int sU = ss;
        int sF = Ls * sU;
        DhopSiteDag(Opt,st,U_v,st.CommBuf(),sF,sU,Ls,1,in_v,out_v);
      });
    }
  }
   
  template <bool EnableBool = true> static accelerator
  typename std::enable_if<Impl::Dimension == 3 && Nc == 3 &&EnableBool, void>::type
  DhopSite(int Opt,StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
	   int sF, int sU, int Ls, int Nsite, 
	   const FermionFieldView &in, FermionFieldView &out,int interior=1,int exterior=1) 
  {
    //bgq_l1p_optimisation(1);
    switch(Opt) {

#if defined(AVX512) || defined (QPX)
    case WilsonKernelsStatic::OptInlineAsm:
      if(interior&&exterior) WilsonKernels<Impl>::AsmDhopSite   (st,U,buf,sF,sU,Ls,Nsite,in,out);
      else if (interior)     WilsonKernels<Impl>::AsmDhopSiteInt(st,U,buf,sF,sU,Ls,Nsite,in,out);
      else if (exterior)     WilsonKernels<Impl>::AsmDhopSiteExt(st,U,buf,sF,sU,Ls,Nsite,in,out);
      else assert(0);
      break;
#endif
    case WilsonKernelsStatic::OptHandUnroll:
      for (int site = 0; site < Nsite; site++) {
	for (int s = 0; s < Ls; s++) {
	  if(interior&&exterior) WilsonKernels<Impl>::HandDhopSite(st,U,buf,sF,sU,in,out);
	  else if (interior)     WilsonKernels<Impl>::HandDhopSiteInt(st,U,buf,sF,sU,in,out);
	  else if (exterior)     WilsonKernels<Impl>::HandDhopSiteExt(st,U,buf,sF,sU,in,out);
	  sF++;
	}
	sU++;
      }
      break;
    case WilsonKernelsStatic::OptGpu:
    case WilsonKernelsStatic::OptGeneric:
      for (int site = 0; site < Nsite; site++) {
	for (int s = 0; s < Ls; s++) {
	  if(interior&&exterior) WilsonKernels<Impl>::GenericDhopSite(st,U,buf,sF,sU,in,out);
	  else if (interior)     WilsonKernels<Impl>::GenericDhopSiteInt(st,U,buf,sF,sU,in,out);
	  else if (exterior)     WilsonKernels<Impl>::GenericDhopSiteExt(st,U,buf,sF,sU,in,out);
	  else assert(0);
	  sF++;
	}
	sU++;
      }
      break;
    default:
      assert(0);
    }
    //bgq_l1p_optimisation(0);
  }
     
  template <bool EnableBool = true> static accelerator
  typename std::enable_if<(Impl::Dimension != 3 || (Impl::Dimension == 3 && Nc != 3)) && EnableBool, void>::type
  DhopSite(int Opt, StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
	   int sF, int sU, int Ls, int Nsite, const FermionFieldView &in, FermionFieldView &out,int interior=1,int exterior=1 ) {
    // no kernel choice  
    for (int site = 0; site < Nsite; site++) {
      for (int s = 0; s < Ls; s++) {
	if(interior&&exterior) WilsonKernels<Impl>::GenericDhopSite(st,U,buf,sF,sU,in,out);
	else if (interior)     WilsonKernels<Impl>::GenericDhopSiteInt(st,U,buf,sF,sU,in,out);
	else if (exterior)     WilsonKernels<Impl>::GenericDhopSiteExt(st,U,buf,sF,sU,in,out);
	else assert(0);
	sF++;
      }
      sU++;
    }
  }
     
  template <bool EnableBool = true> static accelerator
  typename std::enable_if<Impl::Dimension == 3 && Nc == 3 && EnableBool,void>::type
  DhopSiteDag(int Opt, StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
	      int sF, int sU, int Ls, int Nsite, const FermionFieldView &in, FermionFieldView &out,int interior=1,int exterior=1) 
  {
    //bgq_l1p_optimisation(1);
    switch(Opt) {
#if defined(AVX512) || defined (QPX)
    case WilsonKernelsStatic::OptInlineAsm:
      if(interior&&exterior) WilsonKernels<Impl>::AsmDhopSiteDag   (st,U,buf,sF,sU,Ls,Nsite,in,out);
      else if (interior)     WilsonKernels<Impl>::AsmDhopSiteDagInt(st,U,buf,sF,sU,Ls,Nsite,in,out);
      else if (exterior)     WilsonKernels<Impl>::AsmDhopSiteDagExt(st,U,buf,sF,sU,Ls,Nsite,in,out);
      else assert(0);
      break;
#endif
    case WilsonKernelsStatic::OptHandUnroll:
      for (int site = 0; site < Nsite; site++) {
	for (int s = 0; s < Ls; s++) {
	  if(interior&&exterior) WilsonKernels<Impl>::HandDhopSiteDag(st,U,buf,sF,sU,in,out);
	  else if (interior)     WilsonKernels<Impl>::HandDhopSiteDagInt(st,U,buf,sF,sU,in,out);
	  else if (exterior)     WilsonKernels<Impl>::HandDhopSiteDagExt(st,U,buf,sF,sU,in,out);
	  else assert(0);
	  sF++;
	}
	sU++;
      }
      break;
    case WilsonKernelsStatic::OptGpu:
    case WilsonKernelsStatic::OptGeneric:
      for (int site = 0; site < Nsite; site++) {
	for (int s = 0; s < Ls; s++) {
	  if(interior&&exterior) WilsonKernels<Impl>::GenericDhopSiteDag(st,U,buf,sF,sU,in,out);
	  else if (interior)     WilsonKernels<Impl>::GenericDhopSiteDagInt(st,U,buf,sF,sU,in,out);
	  else if (exterior)     WilsonKernels<Impl>::GenericDhopSiteDagExt(st,U,buf,sF,sU,in,out);
	  else assert(0);
	  sF++;
	}
	sU++;
      }
      break;
    default:
      assert(0);
    }
    //bgq_l1p_optimisation(0);
  }

  template <bool EnableBool = true> static accelerator
  typename std::enable_if<(Impl::Dimension != 3 || (Impl::Dimension == 3 && Nc != 3)) && EnableBool,void>::type
  DhopSiteDag(int Opt,StencilView &st,  DoubledGaugeFieldView &U,SiteHalfSpinor * buf,
	      int sF, int sU, int Ls, int Nsite, const FermionFieldView &in, FermionFieldView &out,int interior=1,int exterior=1) {

    for (int site = 0; site < Nsite; site++) {
      for (int s = 0; s < Ls; s++) {
	if(interior&&exterior) WilsonKernels<Impl>::GenericDhopSiteDag(st,U,buf,sF,sU,in,out);
	else if (interior)     WilsonKernels<Impl>::GenericDhopSiteDagInt(st,U,buf,sF,sU,in,out);
	else if (exterior)     WilsonKernels<Impl>::GenericDhopSiteDagExt(st,U,buf,sF,sU,in,out);
	else assert(0);
	sF++;
      }
      sU++;
    }
  }

  static accelerator void DhopDirK(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor * buf,
				   int sF, int sU, const FermionFieldView &in, FermionFieldView &out, int dirdisp, int gamma);
      
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
                                  vInteger t_mask,
                                  bool switch_sign = false);

  static void SeqConservedCurrentSiteBwd(const SitePropagator &q_in,
                                  SitePropagator &q_out,
                                  DoubledGaugeFieldView &U,
                                  unsigned int sU,
                                  unsigned int mu,
                                  vInteger t_mask,
                                  bool switch_sign = false);

private:
  // Specialised variants
  static accelerator void GpuDhopSite(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
				      int sF,  int sU, const FermionFieldView &in, FermionFieldView &out);
  
  static accelerator void GpuDhopSiteDag(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
					 int sF, int sU, const FermionFieldView &in, FermionFieldView &out);

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
  
  static accelerator void AsmDhopSite(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
				      int sF, int sU, int Ls, int Nsite, const FermionFieldView &in,FermionFieldView &out);
  
  static accelerator void AsmDhopSiteDag(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
					 int sF, int sU, int Ls, int Nsite, const FermionFieldView &in, FermionFieldView &out);
  
  static accelerator void AsmDhopSiteInt(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
					 int sF, int sU, int Ls, int Nsite, const FermionFieldView &in,FermionFieldView &out);
  
  static accelerator void AsmDhopSiteDagInt(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
					    int sF, int sU, int Ls, int Nsite, const FermionFieldView &in, FermionFieldView &out);
  
  static accelerator void AsmDhopSiteExt(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
					 int sF, int sU, int Ls, int Nsite, const FermionFieldView &in,FermionFieldView &out);
  
  static accelerator void AsmDhopSiteDagExt(StencilView &st,  DoubledGaugeFieldView &U, SiteHalfSpinor * buf,
					    int sF, int sU, int Ls, int Nsite, const FermionFieldView &in, FermionFieldView &out);
  

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

};
    
NAMESPACE_END(Grid);

#endif
