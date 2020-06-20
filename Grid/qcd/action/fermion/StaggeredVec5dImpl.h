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

template <class S, class Representation = FundamentalRepresentation >
class StaggeredVec5dImpl : public PeriodicGaugeImpl<GaugeImplTypes<S, Representation::Dimension > > {

public:

  static const int Dimension = Representation::Dimension;
    static const bool isFundamental = Representation::isFundamental;
  static const bool LsVectorised=true;
  typedef RealD   Coeff_t ;
  typedef PeriodicGaugeImpl<GaugeImplTypes<S, Dimension > > Gimpl;
      
  //Necessary?
  constexpr bool is_fundamental() const{return Dimension == Nc ? 1 : 0;}


  INHERIT_GIMPL_TYPES(Gimpl);

  template <typename vtype> using iImplSpinor            = iScalar<iScalar<iVector<vtype, Dimension> > >;
  template <typename vtype> using iImplHalfSpinor        = iScalar<iScalar<iVector<vtype, Dimension> > >;
  template <typename vtype> using iImplDoubledGaugeField = iVector<iScalar<iMatrix<vtype, Dimension> >, Nds>;
  template <typename vtype> using iImplGaugeField        = iVector<iScalar<iMatrix<vtype, Dimension> >, Nd>;
  template <typename vtype> using iImplGaugeLink         = iScalar<iScalar<iMatrix<vtype, Dimension> > >;
  template <typename vtype> using iImplPropagator        = iScalar<iScalar<iMatrix<vtype, Dimension> > >;

  // Make the doubled gauge field a *scalar*
  typedef iImplDoubledGaugeField<typename Simd::scalar_type>  SiteDoubledGaugeField;  // This is a scalar
  typedef iImplGaugeField<typename Simd::scalar_type>         SiteScalarGaugeField;  // scalar
  typedef iImplGaugeLink<typename Simd::scalar_type>          SiteScalarGaugeLink;  // scalar
  typedef iImplPropagator<Simd>        SitePropagator;

  typedef Lattice<SiteDoubledGaugeField> DoubledGaugeField;
  typedef Lattice<SitePropagator> PropagatorField;
    
  typedef iImplSpinor<Simd>            SiteSpinor;
  typedef iImplHalfSpinor<Simd>        SiteHalfSpinor;

    
  typedef Lattice<SiteSpinor>            FermionField;
    
  typedef StaggeredImplParams ImplParams;
  typedef SimpleCompressor<SiteSpinor> Compressor;
  typedef CartesianStencil<SiteSpinor, SiteSpinor, ImplParams> StencilImpl;
  typedef typename StencilImpl::View_type StencilView;
    
  ImplParams Params;
    
  StaggeredVec5dImpl(const ImplParams &p = ImplParams()) : Params(p){};

  template <class ref>
  static accelerator_inline void loadLinkElement(Simd &reg, ref &memory) 
  {
    vsplat(reg, memory);
  }

  static accelerator_inline void multLink(SiteHalfSpinor &phi, 
					  const SiteDoubledGaugeField &U,
					  const SiteHalfSpinor &chi, 
					  int mu) 
  {
    SiteGaugeLink UU;
    for (int i = 0; i < Dimension; i++) {
      for (int j = 0; j < Dimension; j++) {
	vsplat(UU()()(i, j), U(mu)()(i, j));
      }
    }
    mult(&phi(), &UU(), &chi());
  }
  static accelerator_inline void multLinkAdd(SiteHalfSpinor &phi, 
					     const SiteDoubledGaugeField &U,
					     const SiteHalfSpinor &chi, 
					     int mu) 
  {
    SiteGaugeLink UU;
    for (int i = 0; i < Dimension; i++) {
      for (int j = 0; j < Dimension; j++) {
	vsplat(UU()()(i, j), U(mu)()(i, j));
      }
    }
    mac(&phi(), &UU(), &chi());
  }
      
  inline void InsertGaugeField(DoubledGaugeField &U_ds,const GaugeLinkField &U,int mu)
  {
    assert(0);
  }
  inline void DoubleStore(GridBase *GaugeGrid,
			  DoubledGaugeField &UUUds, // for Naik term
			  DoubledGaugeField &Uds,
			  const GaugeField &Uthin,
			  const GaugeField &Ufat) 
  {

    GridBase * InputGrid = Uthin.Grid();
    conformable(InputGrid,Ufat.Grid());

    GaugeLinkField U(InputGrid);
    GaugeLinkField UU(InputGrid);
    GaugeLinkField UUU(InputGrid);
    GaugeLinkField Udag(InputGrid);
    GaugeLinkField UUUdag(InputGrid);

    for (int mu = 0; mu < Nd; mu++) {

      // Staggered Phase.
      Lattice<iScalar<vInteger> > coor(InputGrid);
      Lattice<iScalar<vInteger> > x(InputGrid); LatticeCoordinate(x,0);
      Lattice<iScalar<vInteger> > y(InputGrid); LatticeCoordinate(y,1);
      Lattice<iScalar<vInteger> > z(InputGrid); LatticeCoordinate(z,2);
      Lattice<iScalar<vInteger> > t(InputGrid); LatticeCoordinate(t,3);

      Lattice<iScalar<vInteger> > lin_z(InputGrid); lin_z=x+y;
      Lattice<iScalar<vInteger> > lin_t(InputGrid); lin_t=x+y+z;

      ComplexField phases(InputGrid);	phases=1.0;

      if ( mu == 1 ) phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
      if ( mu == 2 ) phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
      if ( mu == 3 ) phases = where( mod(lin_t,2)==(Integer)0, phases,-phases);

      // 1 hop based on fat links
      U      = PeekIndex<LorentzIndex>(Ufat, mu);
      Udag   = adj( Cshift(U, mu, -1));

      U    = U    *phases;
      Udag = Udag *phases;

      InsertGaugeField(Uds,U,mu);
      InsertGaugeField(Uds,Udag,mu+4);

      // 3 hop based on thin links. Crazy huh ?
      U  = PeekIndex<LorentzIndex>(Uthin, mu);
      UU = Gimpl::CovShiftForward(U,mu,U);
      UUU= Gimpl::CovShiftForward(U,mu,UU);
	
      UUUdag = adj( Cshift(UUU, mu, -3));

      UUU    = UUU    *phases;
      UUUdag = UUUdag *phases;

      InsertGaugeField(UUUds,UUU,mu);
      InsertGaugeField(UUUds,UUUdag,mu+4);

    }
  }

  inline void InsertForce4D(GaugeField &mat, FermionField &Btilde, FermionField &A,int mu){
    assert(0);
  }   
      
  inline void InsertForce5D(GaugeField &mat, FermionField &Btilde, FermionField &Atilde,int mu){
    assert (0); 
  }
};
typedef StaggeredVec5dImpl<vComplex,  FundamentalRepresentation > StaggeredVec5dImplR;   // Real.. whichever prec
typedef StaggeredVec5dImpl<vComplexF, FundamentalRepresentation > StaggeredVec5dImplF;  // Float
typedef StaggeredVec5dImpl<vComplexD, FundamentalRepresentation > StaggeredVec5dImplD;  // Double

NAMESPACE_END(Grid);
