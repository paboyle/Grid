/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/FermionOperatorImpl.h

Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>

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

template <class S, class Representation = FundamentalRepresentation, class Options=CoeffReal>
class GparityWilsonImpl : public ConjugateGaugeImpl<GaugeImplTypes<S, Representation::Dimension> > {
public:

 static const int Dimension = Representation::Dimension;
 static const bool isFundamental = Representation::isFundamental;
 static const int Nhcs = Options::Nhcs;
 static const bool LsVectorised=false;
 static const bool isGparity=true;

 typedef ConjugateGaugeImpl< GaugeImplTypes<S,Dimension> > Gimpl;
 INHERIT_GIMPL_TYPES(Gimpl);
 
 typedef typename Options::_Coeff_t Coeff_t;
 typedef typename Options::template PrecisionMapper<Simd>::LowerPrecVector SimdL;
      
 template <typename vtype> using iImplSpinor            = iVector<iVector<iVector<vtype, Dimension>, Ns>,   Ngp>;
 template <typename vtype> using iImplPropagator        = iMatrix<iMatrix<iMatrix<vtype, Dimension>, Ns>,   Ngp>;
 template <typename vtype> using iImplHalfSpinor        = iVector<iVector<iVector<vtype, Dimension>, Nhs>,  Ngp>;
 template <typename vtype> using iImplHalfCommSpinor    = iVector<iVector<iVector<vtype, Dimension>, Nhcs>, Ngp>;
 template <typename vtype> using iImplDoubledGaugeField = iVector<iVector<iScalar<iMatrix<vtype, Dimension> >, Nds>, Ngp>;

  typedef iImplSpinor<Simd>            SiteSpinor;
  typedef iImplPropagator<Simd>        SitePropagator;
  typedef iImplHalfSpinor<Simd>        SiteHalfSpinor;
  typedef iImplHalfCommSpinor<SimdL>   SiteHalfCommSpinor;
  typedef iImplDoubledGaugeField<Simd> SiteDoubledGaugeField;

  typedef Lattice<SiteSpinor> FermionField;
  typedef Lattice<SitePropagator> PropagatorField;
  typedef Lattice<SiteDoubledGaugeField> DoubledGaugeField;
 
  typedef GparityWilsonImplParams ImplParams;
  typedef WilsonCompressor<SiteHalfCommSpinor,SiteHalfSpinor, SiteSpinor> Compressor;
  typedef WilsonStencil<SiteSpinor, SiteHalfSpinor, ImplParams> StencilImpl;
  typedef typename StencilImpl::View_type StencilView;
      
  ImplParams Params;

  GparityWilsonImpl(const ImplParams &p = ImplParams()) : Params(p){};

  // provide the multiply by link that is differentiated between Gparity (with
  // flavour index) and non-Gparity
  template<class _Spinor>
  static accelerator_inline void multLink(_Spinor &phi, 
					  const SiteDoubledGaugeField &U,
					  const _Spinor &chi, 
					  int mu) 
  {
    assert(0);
  } 

  template<class _Spinor>
  static accelerator_inline void multLink(_Spinor &phi, 
					  const SiteDoubledGaugeField &U,
					  const _Spinor &chi, 
					  int mu, 
					  StencilEntry *SE,
					  StencilView &St) 
  {
    int direction = St._directions[mu];
    int distance  = St._distances[mu];
    int ptype     = St._permute_type[mu];
    int sl        = St._simd_layout[direction];
    Coordinate icoor;

#ifdef __CUDA_ARCH__
    _Spinor tmp;

    const int Nsimd =SiteDoubledGaugeField::Nsimd();
    int s = SIMTlane(Nsimd);
    St.iCoorFromIindex(icoor,s);

    int mmu = mu % Nd;
    if ( SE->_around_the_world && St.parameters.twists[mmu] ) {
      
      int permute_lane = (sl==1) 
    	|| ((distance== 1)&&(icoor[direction]==1))
	|| ((distance==-1)&&(icoor[direction]==0));

      if ( permute_lane ) { 
	tmp(0) = chi(1);
	tmp(1) = chi(0);
      } else {
	tmp(0) = chi(0);
	tmp(1) = chi(1);
      }

      auto UU0=coalescedRead(U(0)(mu));
      auto UU1=coalescedRead(U(1)(mu));

      mult(&phi(0),&UU0,&tmp(0));
      mult(&phi(1),&UU1,&tmp(1));

    } else {

      auto UU0=coalescedRead(U(0)(mu));
      auto UU1=coalescedRead(U(1)(mu));

      mult(&phi(0),&UU0,&chi(0));
      mult(&phi(1),&UU1,&chi(1));

    }

#else
    typedef _Spinor vobj;
    typedef typename SiteHalfSpinor::scalar_object sobj;
    typedef typename SiteHalfSpinor::vector_type   vector_type;
	
    vobj vtmp;
    sobj stmp;
        
    const int Nsimd =vector_type::Nsimd();
    
    // Fixme X.Y.Z.T hardcode in stencil
    int mmu = mu % Nd;
        
    // assert our assumptions
    assert((distance == 1) || (distance == -1));  // nearest neighbour stencil hard code
    assert((sl == 1) || (sl == 2));

    if ( SE->_around_the_world && St.parameters.twists[mmu] ) {

      if ( sl == 2 ) {
       
	ExtractBuffer<sobj> vals(Nsimd);

	extract(chi,vals);
	for(int s=0;s<Nsimd;s++){

	  St.iCoorFromIindex(icoor,s);
              
	  assert((icoor[direction]==0)||(icoor[direction]==1));
              
	  int permute_lane;
	  if ( distance == 1) {
	    permute_lane = icoor[direction]?1:0;
	  } else {
	    permute_lane = icoor[direction]?0:1;
	  }
              
	  if ( permute_lane ) { 
	    stmp(0) = vals[s](1);
	    stmp(1) = vals[s](0);
	    vals[s] = stmp;
	  }
	}
	merge(vtmp,vals);
            
      } else { 
	vtmp(0) = chi(1);
	vtmp(1) = chi(0);
      }
      mult(&phi(0),&U(0)(mu),&vtmp(0));
      mult(&phi(1),&U(1)(mu),&vtmp(1));
     
    } else { 
      mult(&phi(0),&U(0)(mu),&chi(0));
      mult(&phi(1),&U(1)(mu),&chi(1));
    }
#endif   
  }


  template<class _SpinorField>
  inline void multLinkField(_SpinorField & out,
			    const DoubledGaugeField &Umu,
			    const _SpinorField & phi,
			    int mu)
  {
    assert(0);
  }

  template <class ref>
  static accelerator_inline void loadLinkElement(Simd &reg, ref &memory) 
  {
    reg = memory;
  }

  inline void DoubleStore(GridBase *GaugeGrid,DoubledGaugeField &Uds,const GaugeField &Umu)
  {
    conformable(Uds.Grid(),GaugeGrid);
    conformable(Umu.Grid(),GaugeGrid);
   
    GaugeLinkField Utmp (GaugeGrid);
    GaugeLinkField U    (GaugeGrid);
    GaugeLinkField Uconj(GaugeGrid);
   
    Lattice<iScalar<vInteger> > coor(GaugeGrid);
        
    for(int mu=0;mu<Nd;mu++){
          
      LatticeCoordinate(coor,mu);
          
      U     = PeekIndex<LorentzIndex>(Umu,mu);
      Uconj = conjugate(U);
     
      // This phase could come from a simple bc 1,1,-1,1 ..
      int neglink = GaugeGrid->GlobalDimensions()[mu]-1;
      if ( Params.twists[mu] ) { 
	Uconj = where(coor==neglink,-Uconj,Uconj);
      }
	  
      auto U_v = U.View();
      auto Uds_v = Uds.View();
      auto Uconj_v = Uconj.View();
      auto Utmp_v= Utmp.View();
      thread_foreach(ss,U_v,{
	Uds_v[ss](0)(mu) = U_v[ss]();
	Uds_v[ss](1)(mu) = Uconj_v[ss]();
      });
          
      U     = adj(Cshift(U    ,mu,-1));      // correct except for spanning the boundary
      Uconj = adj(Cshift(Uconj,mu,-1));
 
      Utmp = U;
      if ( Params.twists[mu] ) { 
	Utmp = where(coor==0,Uconj,Utmp);
      }

      thread_foreach(ss,Utmp_v,{
	Uds_v[ss](0)(mu+4) = Utmp_v[ss]();
      });
          
      Utmp = Uconj;
      if ( Params.twists[mu] ) { 
	Utmp = where(coor==0,U,Utmp);
      }
	  
      thread_foreach(ss,Utmp_v,{
        Uds_v[ss](1)(mu+4) = Utmp_v[ss]();
      });
          
    }
  }
      
  inline void InsertForce4D(GaugeField &mat, FermionField &Btilde, FermionField &A, int mu) {

    // DhopDir provides U or Uconj depending on coor/flavour.
    GaugeLinkField link(mat.Grid());
    // use lorentz for flavour as hack.
    auto tmp = TraceIndex<SpinIndex>(outerProduct(Btilde, A));
    auto link_v = link.View();
    auto tmp_v = tmp.View();
    thread_foreach(ss,tmp_v,{
      link_v[ss]() = tmp_v[ss](0, 0) + conjugate(tmp_v[ss](1, 1));
    });
    PokeIndex<LorentzIndex>(mat, link, mu);
    return;
  }
      
 inline void outerProductImpl(PropagatorField &mat, const FermionField &Btilde, const FermionField &A){
   //mat = outerProduct(Btilde, A);
   assert(0);
  }

  inline void TraceSpinImpl(GaugeLinkField &mat, PropagatorField&P) {
    assert(0);
    /*
    auto tmp = TraceIndex<SpinIndex>(P);
    parallel_for(auto ss = tmp.begin(); ss < tmp.end(); ss++) {
      mat[ss]() = tmp[ss](0, 0) + conjugate(tmp[ss](1, 1));
    }
    */
  }

  inline void extractLinkField(std::vector<GaugeLinkField> &mat, DoubledGaugeField &Uds){
    assert(0);
  }
  
  inline void InsertForce5D(GaugeField &mat, FermionField &Btilde, FermionField &Atilde, int mu) {

    int Ls = Btilde.Grid()->_fdimensions[0];
        
    GaugeLinkField tmp(mat.Grid());
    tmp = Zero();
    auto tmp_v = tmp.View();
    auto Atilde_v = Atilde.View();
    auto Btilde_v = Btilde.View();
    thread_for(ss,tmp.Grid()->oSites(),{
      for (int s = 0; s < Ls; s++) {
	int sF = s + Ls * ss;
	auto ttmp = traceIndex<SpinIndex>(outerProduct(Btilde_v[sF], Atilde_v[sF]));
	tmp_v[ss]() = tmp_v[ss]() + ttmp(0, 0) + conjugate(ttmp(1, 1));
      }
    });
    PokeIndex<LorentzIndex>(mat, tmp, mu);
    return;
  }
  
};

typedef GparityWilsonImpl<vComplex , FundamentalRepresentation,CoeffReal> GparityWilsonImplR;  // Real.. whichever prec
typedef GparityWilsonImpl<vComplexF, FundamentalRepresentation,CoeffReal> GparityWilsonImplF;  // Float
typedef GparityWilsonImpl<vComplexD, FundamentalRepresentation,CoeffReal> GparityWilsonImplD;  // Double
 
typedef GparityWilsonImpl<vComplex , FundamentalRepresentation,CoeffRealHalfComms> GparityWilsonImplRL;  // Real.. whichever prec
typedef GparityWilsonImpl<vComplexF, FundamentalRepresentation,CoeffRealHalfComms> GparityWilsonImplFH;  // Float
typedef GparityWilsonImpl<vComplexD, FundamentalRepresentation,CoeffRealHalfComms> GparityWilsonImplDF;  // Double

NAMESPACE_END(Grid);
