    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 



    Source file: ./lib/qcd/action/fermion/WilsonKernelsAsm.cc

    Copyright (C) 2015

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

#include <Grid.h>

namespace Grid {
namespace QCD {


  ///////////////////////////////////////////////////////////
  // Default to no assembler implementation
  ///////////////////////////////////////////////////////////
template<class Impl>
void WilsonKernels<Impl >::DiracOptAsmDhopSite(StencilImpl &st,DoubledGaugeField &U,
					       std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
					       int ss,int sU,const FermionField &in, FermionField &out)
{
  assert(0);
}

#if defined(AVX512) 


  ///////////////////////////////////////////////////////////
  // If we are AVX512 specialise the single precision routine
  ///////////////////////////////////////////////////////////

#include <simd/Intel512wilson.h>
#include <simd/Intel512single.h>

static Vector<vComplexF> signs;

int setupSigns(void ){
  Vector<vComplexF> bother(2);
  signs = bother;
  vrsign(signs[0]);
  visign(signs[1]);
  return 1;
}
static int signInit = setupSigns();

#define MAYBEPERM(A,perm) if (perm) { A ; }
#define MULT_2SPIN(ptr,pf) MULT_ADDSUB_2SPIN(ptr)

template<>
void WilsonKernels<WilsonImplF>::DiracOptAsmDhopSite(StencilImpl &st,DoubledGaugeField &U,
					       std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
					       int ss,int sU,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef VMOVIDUP
#undef VMOVRDUP
#undef MAYBEPERM
#undef MULT_2SPIN
#define MAYBEPERM(A,B) 
#define VMOVIDUP(A,B,C)                                  VBCASTIDUPf(A,B,C)
#define VMOVRDUP(A,B,C)                                  VBCASTRDUPf(A,B,C)
#define MULT_2SPIN(ptr,pf) MULT_ADDSUB_2SPIN_LS(ptr,pf)
template<>
void WilsonKernels<DomainWallRedBlack5dImplF>::DiracOptAsmDhopSite(StencilImpl &st,DoubledGaugeField &U,
								   std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
								   int ss,int sU,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#endif

template class WilsonKernels<WilsonImplF>;		
template class WilsonKernels<WilsonImplD>; 
template class WilsonKernels<GparityWilsonImplF>;
template class WilsonKernels<GparityWilsonImplD>;
template class WilsonKernels<DomainWallRedBlack5dImplF>;
template class WilsonKernels<DomainWallRedBlack5dImplD>;
}}

