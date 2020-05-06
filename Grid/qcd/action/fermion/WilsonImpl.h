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

  
/////////////////////////////////////////////////////////////////////////////
// Single flavour four spinors with colour index
/////////////////////////////////////////////////////////////////////////////
template <class S, class Representation = FundamentalRepresentation,class Options = CoeffReal >
class WilsonImpl : public PeriodicGaugeImpl<GaugeImplTypes<S, Representation::Dimension > > {
public:

  static const int Dimension = Representation::Dimension;
  static const bool isFundamental = Representation::isFundamental;
  static const bool LsVectorised=false;
  static const bool isGparity=false;
  static const int Nhcs = Options::Nhcs;

  typedef PeriodicGaugeImpl<GaugeImplTypes<S, Dimension > > Gimpl;
  INHERIT_GIMPL_TYPES(Gimpl);
      
  //Necessary?
  constexpr bool is_fundamental() const{return Dimension == Nc ? 1 : 0;}
    
  typedef typename Options::_Coeff_t Coeff_t;
  typedef typename Options::template PrecisionMapper<Simd>::LowerPrecVector SimdL;
      
  template <typename vtype> using iImplSpinor            = iScalar<iVector<iVector<vtype, Dimension>, Ns> >;
  template <typename vtype> using iImplPropagator        = iScalar<iMatrix<iMatrix<vtype, Dimension>, Ns> >;
  template <typename vtype> using iImplHalfSpinor        = iScalar<iVector<iVector<vtype, Dimension>, Nhs> >;
  template <typename vtype> using iImplHalfCommSpinor    = iScalar<iVector<iVector<vtype, Dimension>, Nhcs> >;
  template <typename vtype> using iImplDoubledGaugeField = iVector<iScalar<iMatrix<vtype, Dimension> >, Nds>;
    
  typedef iImplSpinor<Simd>            SiteSpinor;
  typedef iImplPropagator<Simd>        SitePropagator;
  typedef iImplHalfSpinor<Simd>        SiteHalfSpinor;
  typedef iImplHalfCommSpinor<SimdL>   SiteHalfCommSpinor;
  typedef iImplDoubledGaugeField<Simd> SiteDoubledGaugeField;
    
  typedef Lattice<SiteSpinor>            FermionField;
  typedef Lattice<SitePropagator>        PropagatorField;
  typedef Lattice<SiteDoubledGaugeField> DoubledGaugeField;
    
  typedef WilsonCompressor<SiteHalfCommSpinor,SiteHalfSpinor, SiteSpinor> Compressor;
  typedef WilsonImplParams ImplParams;
  typedef WilsonStencil<SiteSpinor, SiteHalfSpinor,ImplParams> StencilImpl;
  typedef typename StencilImpl::View_type StencilView;
    
  ImplParams Params;

  WilsonImpl(const ImplParams &p = ImplParams()) : Params(p){
    assert(Params.boundary_phases.size() == Nd);
  };

  template<class _Spinor>
  static accelerator_inline void multLink(_Spinor &phi,
					  const SiteDoubledGaugeField &U,
					  const _Spinor &chi,
					  int mu) 
  {
    auto UU = coalescedRead(U(mu));
    mult(&phi(), &UU, &chi());
  }
  template<class _Spinor>
  static accelerator_inline void multLink(_Spinor &phi,
					  const SiteDoubledGaugeField &U,
					  const _Spinor &chi,
					  int mu,
					  StencilEntry *SE,
					  StencilView &St) 
  {
    multLink(phi,U,chi,mu);
  }

  template<class _SpinorField> 
  inline void multLinkField(_SpinorField & out,
			    const DoubledGaugeField &Umu,
			    const _SpinorField & phi,
			    int mu)
  {
    auto out_v= out.View();
    auto phi_v= phi.View();
    auto Umu_v= Umu.View();
    thread_for(sss,out.Grid()->oSites(),{
	multLink(out_v[sss],Umu_v[sss],phi_v[sss],mu);
    });
  }
					   
  template <class ref>
  static accelerator_inline void loadLinkElement(Simd &reg, ref &memory) 
  {
    reg = memory;
  }
      
  inline void DoubleStore(GridBase *GaugeGrid,
			  DoubledGaugeField &Uds,
			  const GaugeField &Umu) 
  {
    typedef typename Simd::scalar_type scalar_type;

    conformable(Uds.Grid(), GaugeGrid);
    conformable(Umu.Grid(), GaugeGrid);

    GaugeLinkField U(GaugeGrid);
    GaugeLinkField tmp(GaugeGrid);

    Lattice<iScalar<vInteger> > coor(GaugeGrid);
      ////////////////////////////////////////////////////
      // apply any boundary phase or twists
      ////////////////////////////////////////////////////
    for (int mu = 0; mu < Nd; mu++) {

	////////// boundary phase /////////////
      auto pha = Params.boundary_phases[mu];
      scalar_type phase( real(pha),imag(pha) );

	int L   = GaugeGrid->GlobalDimensions()[mu];
        int Lmu = L - 1;

      LatticeCoordinate(coor, mu);

      U = PeekIndex<LorentzIndex>(Umu, mu);

	// apply any twists
	RealD theta = Params.twist_n_2pi_L[mu] * 2*M_PI / L;
	if ( theta != 0.0) { 
	  scalar_type twphase(::cos(theta),::sin(theta));
	  U = twphase*U;
	  std::cout << GridLogMessage << " Twist ["<<mu<<"] "<< Params.twist_n_2pi_L[mu]<< " phase"<<phase <<std::endl;
	}

      tmp = where(coor == Lmu, phase * U, U);
      PokeIndex<LorentzIndex>(Uds, tmp, mu);

      U = adj(Cshift(U, mu, -1));
      U = where(coor == 0, conjugate(phase) * U, U); 
      PokeIndex<LorentzIndex>(Uds, U, mu + 4);
    }
  }

  inline void InsertForce4D(GaugeField &mat, FermionField &Btilde, FermionField &A,int mu){
    GaugeLinkField link(mat.Grid());
    link = TraceIndex<SpinIndex>(outerProduct(Btilde,A)); 
    PokeIndex<LorentzIndex>(mat,link,mu);
  }   
      
    inline void outerProductImpl(PropagatorField &mat, const FermionField &B, const FermionField &A){
      mat = outerProduct(B,A); 
    }  

    inline void TraceSpinImpl(GaugeLinkField &mat, PropagatorField&P) {
      mat = TraceIndex<SpinIndex>(P); 
    }
      
    inline void extractLinkField(std::vector<GaugeLinkField> &mat, DoubledGaugeField &Uds){
      for (int mu = 0; mu < Nd; mu++)
      mat[mu] = PeekIndex<LorentzIndex>(Uds, mu);
    }


  inline void InsertForce5D(GaugeField &mat, FermionField &Btilde, FermionField &Atilde,int mu){
      
    int Ls=Btilde.Grid()->_fdimensions[0];
    GaugeLinkField tmp(mat.Grid());
    tmp = Zero();
    auto tmp_v = tmp.View();
    auto Btilde_v = Btilde.View();
    auto Atilde_v = Atilde.View();
    thread_for(sss,tmp.Grid()->oSites(),{
      int sU=sss;
      for(int s=0;s<Ls;s++){
	int sF = s+Ls*sU;
	tmp_v[sU] = tmp_v[sU]+ traceIndex<SpinIndex>(outerProduct(Btilde_v[sF],Atilde_v[sF])); // ordering here
      }
    });
    PokeIndex<LorentzIndex>(mat,tmp,mu);
      
  }
};


typedef WilsonImpl<vComplex,  FundamentalRepresentation, CoeffReal > WilsonImplR;  // Real.. whichever prec
typedef WilsonImpl<vComplexF, FundamentalRepresentation, CoeffReal > WilsonImplF;  // Float
typedef WilsonImpl<vComplexD, FundamentalRepresentation, CoeffReal > WilsonImplD;  // Double

typedef WilsonImpl<vComplex,  FundamentalRepresentation, CoeffRealHalfComms > WilsonImplRL;  // Real.. whichever prec
typedef WilsonImpl<vComplexF, FundamentalRepresentation, CoeffRealHalfComms > WilsonImplFH;  // Float
typedef WilsonImpl<vComplexD, FundamentalRepresentation, CoeffRealHalfComms > WilsonImplDF;  // Double

typedef WilsonImpl<vComplex,  FundamentalRepresentation, CoeffComplex > ZWilsonImplR; // Real.. whichever prec
typedef WilsonImpl<vComplexF, FundamentalRepresentation, CoeffComplex > ZWilsonImplF; // Float
typedef WilsonImpl<vComplexD, FundamentalRepresentation, CoeffComplex > ZWilsonImplD; // Double

typedef WilsonImpl<vComplex,  FundamentalRepresentation, CoeffComplexHalfComms > ZWilsonImplRL; // Real.. whichever prec
typedef WilsonImpl<vComplexF, FundamentalRepresentation, CoeffComplexHalfComms > ZWilsonImplFH; // Float
typedef WilsonImpl<vComplexD, FundamentalRepresentation, CoeffComplexHalfComms > ZWilsonImplDF; // Double
 
typedef WilsonImpl<vComplex,  AdjointRepresentation, CoeffReal > WilsonAdjImplR;   // Real.. whichever prec
typedef WilsonImpl<vComplexF, AdjointRepresentation, CoeffReal > WilsonAdjImplF;  // Float
typedef WilsonImpl<vComplexD, AdjointRepresentation, CoeffReal > WilsonAdjImplD;  // Double
 
typedef WilsonImpl<vComplex,  TwoIndexSymmetricRepresentation, CoeffReal > WilsonTwoIndexSymmetricImplR;   // Real.. whichever prec
typedef WilsonImpl<vComplexF, TwoIndexSymmetricRepresentation, CoeffReal > WilsonTwoIndexSymmetricImplF;  // Float
typedef WilsonImpl<vComplexD, TwoIndexSymmetricRepresentation, CoeffReal > WilsonTwoIndexSymmetricImplD;  // Double
 
typedef WilsonImpl<vComplex,  TwoIndexAntiSymmetricRepresentation, CoeffReal > WilsonTwoIndexAntiSymmetricImplR;   // Real.. whichever prec
typedef WilsonImpl<vComplexF, TwoIndexAntiSymmetricRepresentation, CoeffReal > WilsonTwoIndexAntiSymmetricImplF;  // Float
typedef WilsonImpl<vComplexD, TwoIndexAntiSymmetricRepresentation, CoeffReal > WilsonTwoIndexAntiSymmetricImplD;  // Double


NAMESPACE_END(Grid);

