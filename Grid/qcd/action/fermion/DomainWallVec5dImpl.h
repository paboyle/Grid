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

template<class S,class Representation = FundamentalRepresentation, class Options=CoeffReal>
class DomainWallVec5dImpl :  public PeriodicGaugeImpl< GaugeImplTypes< S,Representation::Dimension> > { 
public:

  typedef PeriodicGaugeImpl<GaugeImplTypes<S, Representation::Dimension> > Gimpl;
  INHERIT_GIMPL_TYPES(Gimpl);

  static const int Dimension = Representation::Dimension;
  static const bool isFundamental = Representation::isFundamental;
  static const bool LsVectorised=true;
  static const int Nhcs = Options::Nhcs;
      
  typedef typename Options::_Coeff_t Coeff_t;      
  typedef typename Options::template PrecisionMapper<Simd>::LowerPrecVector SimdL;
  
  template <typename vtype> using iImplSpinor            = iScalar<iVector<iVector<vtype, Dimension>, Ns> >;
  template <typename vtype> using iImplPropagator        = iScalar<iMatrix<iMatrix<vtype, Dimension>, Ns> >;
  template <typename vtype> using iImplHalfSpinor        = iScalar<iVector<iVector<vtype, Dimension>, Nhs> >;
  template <typename vtype> using iImplHalfCommSpinor    = iScalar<iVector<iVector<vtype, Dimension>, Nhcs> >;
  template <typename vtype> using iImplDoubledGaugeField = iVector<iScalar<iMatrix<vtype, Dimension> >, Nds>;
  template <typename vtype> using iImplGaugeField        = iVector<iScalar<iMatrix<vtype, Dimension> >, Nd>;
  template <typename vtype> using iImplGaugeLink         = iScalar<iScalar<iMatrix<vtype, Dimension> > >;
  
  typedef iImplSpinor<Simd>            SiteSpinor;
  typedef iImplPropagator<Simd>        SitePropagator;
  typedef iImplHalfSpinor<Simd>        SiteHalfSpinor;
  typedef iImplHalfCommSpinor<SimdL>   SiteHalfCommSpinor;
  typedef Lattice<SiteSpinor>          FermionField;
  typedef Lattice<SitePropagator>      PropagatorField;

  /////////////////////////////////////////////////
  // Make the doubled gauge field a *scalar*
  /////////////////////////////////////////////////
  typedef iImplDoubledGaugeField<typename Simd::scalar_type>  SiteDoubledGaugeField;  // This is a scalar
  typedef iImplGaugeField<typename Simd::scalar_type>         SiteScalarGaugeField;  // scalar
  typedef iImplGaugeLink<typename Simd::scalar_type>          SiteScalarGaugeLink;  // scalar
  typedef Lattice<SiteDoubledGaugeField>                      DoubledGaugeField;
      
  typedef WilsonCompressor<SiteHalfCommSpinor,SiteHalfSpinor, SiteSpinor> Compressor;
  typedef WilsonImplParams ImplParams;
  typedef WilsonStencil<SiteSpinor, SiteHalfSpinor,ImplParams> StencilImpl;
  typedef typename StencilImpl::View_type StencilView;
  
  ImplParams Params;

  DomainWallVec5dImpl(const ImplParams &p = ImplParams()) : Params(p){};
      
  template <class ref>
  static accelerator_inline void loadLinkElement(Simd &reg, ref &memory) 
  {
    vsplat(reg, memory);
  }

  template<class _Spinor>
  static accelerator_inline void multLink(_Spinor &phi, const SiteDoubledGaugeField &U,
					  const _Spinor &chi, int mu, StencilEntry *SE,
					  StencilView &St) 
  {
#ifdef GPU_VEC
    // Gauge link is scalarised
    mult(&phi(), &U(mu), &chi());
#else
    SiteGaugeLink UU;
    for (int i = 0; i < Dimension; i++) {
      for (int j = 0; j < Dimension; j++) {
        vsplat(UU()()(i, j), U(mu)()(i, j));
      }
    }
    mult(&phi(), &UU(), &chi());
#endif
  }

  inline void DoubleStore(GridBase *GaugeGrid, DoubledGaugeField &Uds,const GaugeField &Umu) 
  {
    SiteScalarGaugeField  ScalarUmu;
    SiteDoubledGaugeField ScalarUds;
    
    GaugeLinkField U(Umu.Grid());
    GaugeField  Uadj(Umu.Grid());
    for (int mu = 0; mu < Nd; mu++) {
      U = PeekIndex<LorentzIndex>(Umu, mu);
      U = adj(Cshift(U, mu, -1));
      PokeIndex<LorentzIndex>(Uadj, U, mu);
    }

    autoView(Umu_v,Umu,CpuRead);
    autoView(Uadj_v,Uadj,CpuRead);
    autoView(Uds_v,Uds,CpuWrite);
    thread_for( lidx, GaugeGrid->lSites(), {
      Coordinate lcoor;
      GaugeGrid->LocalIndexToLocalCoor(lidx, lcoor);
      
      peekLocalSite(ScalarUmu, Umu_v, lcoor);
      for (int mu = 0; mu < 4; mu++) ScalarUds(mu) = ScalarUmu(mu);
      
      peekLocalSite(ScalarUmu, Uadj_v, lcoor);
      for (int mu = 0; mu < 4; mu++) ScalarUds(mu + 4) = ScalarUmu(mu);
      
      pokeLocalSite(ScalarUds, Uds_v, lcoor);
    });
  }
      
  inline void InsertForce4D(GaugeField &mat, FermionField &Btilde,FermionField &A, int mu) 
  {
    assert(0);
  }

  inline void outerProductImpl(PropagatorField &mat, const FermionField &Btilde, const FermionField &A){
    assert(0);
  } 

  inline void TraceSpinImpl(GaugeLinkField &mat, PropagatorField&P) {
    assert(0);
  }

  inline void extractLinkField(std::vector<GaugeLinkField> &mat, DoubledGaugeField &Uds){
    assert(0);
  }


  inline void InsertForce5D(GaugeField &mat, FermionField &Btilde, FermionField &Atilde, int mu) {

    assert(0);
    // Following lines to be revised after Peter's addition of half prec
    // missing put lane...
    /*
      typedef decltype(traceIndex<SpinIndex>(outerProduct(Btilde[0], Atilde[0]))) result_type;
      unsigned int LLs = Btilde.Grid()->_rdimensions[0];
      conformable(Atilde.Grid(),Btilde.Grid());
      GridBase* grid = mat.Grid();
      GridBase* Bgrid = Btilde.Grid();
      unsigned int dimU = grid->Nd();
      unsigned int dimF = Bgrid->Nd();
      GaugeLinkField tmp(grid); 
      tmp = Zero();
    
      // FIXME 
      // Current implementation works, thread safe, probably suboptimal
      // Passing through the local coordinate for grid transformation
      // the force grid is in general very different from the Ls vectorized grid

      for (int so = 0; so < grid->oSites(); so++) {
      std::vector<typename result_type::scalar_object> vres(Bgrid->Nsimd());
      std::vector<int> ocoor;  grid->oCoorFromOindex(ocoor,so); 
      for (int si = 0; si < tmp.Grid()->iSites(); si++){
      typename result_type::scalar_object scalar_object; scalar_object = Zero();
      std::vector<int> local_coor;      
      std::vector<int> icoor; grid->iCoorFromIindex(icoor,si);
      grid->InOutCoorToLocalCoor(ocoor, icoor, local_coor);
      for (int s = 0; s < LLs; s++) {
      std::vector<int> slocal_coor(dimF);
      slocal_coor[0] = s;
      for (int s4d = 1; s4d< dimF; s4d++) slocal_coor[s4d] = local_coor[s4d-1];
      int sF = Bgrid->oIndexReduced(slocal_coor);  
      assert(sF < Bgrid->oSites());

      extract(traceIndex<SpinIndex>(outerProduct(Btilde[sF], Atilde[sF])), vres); 
      // sum across the 5d dimension
      for (auto v : vres) scalar_object += v;  
      }
      tmp[so].putlane(scalar_object, si);
      }
      }
      PokeIndex<LorentzIndex>(mat, tmp, mu);
    */
  }
};
typedef DomainWallVec5dImpl<vComplex ,FundamentalRepresentation, CoeffReal> DomainWallVec5dImplR; // Real.. whichever prec
typedef DomainWallVec5dImpl<vComplexF,FundamentalRepresentation, CoeffReal> DomainWallVec5dImplF; // Float
typedef DomainWallVec5dImpl<vComplexD,FundamentalRepresentation, CoeffReal> DomainWallVec5dImplD; // Double
 
typedef DomainWallVec5dImpl<vComplex ,FundamentalRepresentation, CoeffRealHalfComms> DomainWallVec5dImplRL; // Real.. whichever prec
typedef DomainWallVec5dImpl<vComplexF,FundamentalRepresentation, CoeffRealHalfComms> DomainWallVec5dImplFH; // Float
typedef DomainWallVec5dImpl<vComplexD,FundamentalRepresentation, CoeffRealHalfComms> DomainWallVec5dImplDF; // Double
 
typedef DomainWallVec5dImpl<vComplex ,FundamentalRepresentation,CoeffComplex> ZDomainWallVec5dImplR; // Real.. whichever prec
typedef DomainWallVec5dImpl<vComplexF,FundamentalRepresentation,CoeffComplex> ZDomainWallVec5dImplF; // Float
typedef DomainWallVec5dImpl<vComplexD,FundamentalRepresentation,CoeffComplex> ZDomainWallVec5dImplD; // Double
 
typedef DomainWallVec5dImpl<vComplex ,FundamentalRepresentation,CoeffComplexHalfComms> ZDomainWallVec5dImplRL; // Real.. whichever prec
typedef DomainWallVec5dImpl<vComplexF,FundamentalRepresentation,CoeffComplexHalfComms> ZDomainWallVec5dImplFH; // Float
typedef DomainWallVec5dImpl<vComplexD,FundamentalRepresentation,CoeffComplexHalfComms> ZDomainWallVec5dImplDF; // Double

NAMESPACE_END(Grid);
