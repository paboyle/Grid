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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#ifndef  GRID_QCD_DHOP_H
#define  GRID_QCD_DHOP_H

namespace Grid {

  namespace QCD {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Helper routines that implement Wilson stencil for a single site.
    // Common to both the WilsonFermion and WilsonFermion5D
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class WilsonKernelsStatic { 
    public:
      // S-direction is INNERMOST and takes no part in the parity.
      static int AsmOpt;  // these are a temporary hack
      static int HandOpt; // these are a temporary hack
    };

    template<class Impl> class WilsonKernels : public FermionOperator<Impl> , public WilsonKernelsStatic { 
    public:

     INHERIT_IMPL_TYPES(Impl);
     typedef FermionOperator<Impl> Base;
     
    public:

     void DiracOptDhopSite(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,
			   std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
			   int sF, int sU,int Ls, int Ns, const FermionField &in, FermionField &out);
      
     void DiracOptDhopSiteDag(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,
			      std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
			      int sF,int sU,int Ls, int Ns, const FermionField &in,FermionField &out);

     void DiracOptDhopDir(StencilImpl &st,DoubledGaugeField &U,
			  std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
			  int sF,int sU,const FermionField &in, FermionField &out,int dirdisp,int gamma);

    private:
     // Specialised variants
     void DiracOptGenericDhopSite(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,
			   std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
			   int sF,int sU, const FermionField &in, FermionField &out);
      
     void DiracOptGenericDhopSiteDag(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,
			      std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
			      int sF,int sU,const FermionField &in,FermionField &out);

     void DiracOptAsmDhopSite(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,
			      std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
			      int sF,int sU,int Ls, int Ns, const FermionField &in, FermionField &out);


     void DiracOptHandDhopSite(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,
			      std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
			      int sF,int sU,const FermionField &in, FermionField &out);
     
     void DiracOptHandDhopSiteDag(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,
				 std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
				 int sF,int sU,const FermionField &in, FermionField &out);
    public:

     WilsonKernels(const ImplParams &p= ImplParams());
     
    };

  }
}
#endif
