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
  typedef const typename StencilImpl::View_type StencilView;
    
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
    const int Nsimd = SiteHalfSpinor::Nsimd();
    autoView( out_v, out, AcceleratorWrite);
    autoView( phi_v, phi, AcceleratorRead);
    autoView( Umu_v, Umu, AcceleratorRead);
    typedef decltype(coalescedRead(out_v[0]))   calcSpinor;
    accelerator_for(sss,out.Grid()->oSites(),Nsimd,{
	calcSpinor tmp;
	multLink(tmp,Umu_v[sss],phi_v(sss),mu);
	coalescedWrite(out_v[sss],tmp);
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
    #ifdef GRID_HIP
    {
      const int Nsimd = SiteSpinor::Nsimd();

      autoView( Btilde_v , Btilde, AcceleratorRead);
      autoView( A_v , A, AcceleratorRead);
      autoView( mat_v , mat, AcceleratorWrite);

      accelerator_for(sss,mat.Grid()->oSites(),Nsimd,{

          for(int ic = 0; ic < Dimension; ic++)
          {
            for(int jc = 0; jc < Dimension; jc++)
            {
              typedef decltype(outerProduct(coalescedRead(Btilde_v[sss]()(0)(ic)),coalescedRead(A_v[sss]()(0)(jc)) ) ) tmpType;
              tmpType tmpentry;
              zeroit(tmpentry);
              for(int spn=0;spn<Ns;spn++){ //sum over spin
                auto bb = coalescedRead(Btilde_v[sss]()(spn)(ic) ); //color vector
                auto aa = coalescedRead(A_v[sss]()(spn)(jc) );
                auto op = outerProduct(bb,aa);
                tmpentry = tmpentry + op;
              }
            coalescedWrite(mat_v[sss](mu)()(ic,jc), tmpentry);
          }
        }
        });
    }
    #else
    GaugeLinkField link(mat.Grid());
    link = TraceIndex<SpinIndex>(outerProduct(Btilde,A)); 
    PokeIndex<LorentzIndex>(mat,link,mu);
    #endif
  }   
      
    inline void outerProductImpl(PropagatorField &mat, const FermionField &B, const FermionField &A){
      #ifdef GRID_HIP
      {
        const int Nsimd = SiteSpinor::Nsimd();
        autoView( B_v , B, AcceleratorRead);
        autoView( A_v , A, AcceleratorRead);
        autoView( mat_v , mat, AcceleratorWrite);

        accelerator_for(sss,mat.Grid()->oSites(),Nsimd,{
        for(int spn=0;spn<Ns;spn++){ //loop over spin
          for(int tpn=0;tpn<Ns;tpn++){ //loop over spin
            for(int ic = 0; ic < Dimension; ic++)
            {
              for(int jc = 0; jc < Dimension; jc++)
              {
                auto bb  = coalescedRead(B_v[sss]()(spn)(ic) ); //color entry
                auto aa  = coalescedRead(A_v[sss]()(tpn)(jc) );
                auto tmp = outerProduct(bb,aa);
                coalescedWrite(mat_v[sss]()(spn,tpn)(ic,jc), tmp);
              }
            }
          }
        }
        });
      }
      #else
      mat = outerProduct(B,A);
      #endif
    }  

    inline void TraceSpinImpl(GaugeLinkField &mat, PropagatorField&P) {
      #ifdef GRID_HIP
      {
        autoView( mat_v , mat, AcceleratorWrite);
        autoView( P_v , P, AcceleratorRead);
        const int Nsimd = SiteSpinor::Nsimd();
        accelerator_for(sss,mat.Grid()->oSites(),Nsimd,{
        for(int ic = 0; ic < Dimension; ic++)
        {
          for(int jc = 0; jc < Dimension; jc++)
          {
              typedef decltype(coalescedRead(P_v[sss]()(0,0)(ic,jc))) tmpType;
              tmpType tmpentry;
              zeroit(tmpentry);
              for(int spn=0;spn<Ns;spn++){ //loop over spin
                auto bb = coalescedRead(P_v[sss]()(spn,spn)(ic,jc) ); //color entry
                tmpentry = tmpentry + bb;
              }
              coalescedWrite(mat_v[sss]()()(ic,jc), tmpentry);
          }
        }
        });
      }
      #else
      mat = TraceIndex<SpinIndex>(P); 
      #endif
    }
      
    inline void extractLinkField(std::vector<GaugeLinkField> &mat, DoubledGaugeField &Uds)
    {
      for (int mu = 0; mu < Nd; mu++)
      {
      #ifndef GRID_HIP
      mat[mu] = PeekIndex<LorentzIndex>(Uds, mu);
      #else
	    {
		   autoView( mat_v, mat[mu], AcceleratorWrite);
		   autoView( U_v, Uds, AcceleratorRead);
		   accelerator_for(ss, mat_v.size(), 1, {
			 for(int ic = 0; ic < Dimension; ic++)
			 {
				for(int jc = 0; jc < Dimension; jc++)
				{
				 auto tmp = coalescedRead(U_v[ss](mu)()(ic,jc));
				 coalescedWrite(mat_v[ss]()()(ic,jc),tmp);
        }
			 }
		   });
	    }
	    #endif
      }
    }

  inline void InsertForce5D(GaugeField &mat, FermionField &Btilde, FermionField &Atilde,int mu)
  {
#undef USE_OLD_INSERT_FORCE    
    int Ls=Btilde.Grid()->_fdimensions[0];
    autoView( mat_v , mat, AcceleratorWrite);
#ifdef USE_OLD_INSERT_FORCE    
    GaugeLinkField tmp(mat.Grid());
    tmp = Zero();
    {
      const int Nsimd = SiteSpinor::Nsimd();
      autoView( tmp_v , tmp, AcceleratorWrite);
      autoView( Btilde_v , Btilde, AcceleratorRead);
      autoView( Atilde_v , Atilde, AcceleratorRead);
      accelerator_for(sss,tmp.Grid()->oSites(),1,{
	  int sU=sss;
	  for(int s=0;s<Ls;s++){
	    int sF = s+Ls*sU;
	    tmp_v[sU] = tmp_v[sU]+ traceIndex<SpinIndex>(outerProduct(Btilde_v[sF],Atilde_v[sF])); // ordering here
	  }
	});
    }
    PokeIndex<LorentzIndex>(mat,tmp,mu);
    #else
    #ifdef GRID_HIP
    {
      const int Nsimd = SiteSpinor::Nsimd();
      autoView( Btilde_v , Btilde, AcceleratorRead);
      autoView( Atilde_v , Atilde, AcceleratorRead);
      accelerator_for(sss,mat.Grid()->oSites(),Nsimd,{
      int sU=sss;

        for(int ic = 0; ic < Dimension; ic++)
        {
          for(int jc = 0; jc < Dimension; jc++)
          {
            typedef decltype(outerProduct(coalescedRead(Btilde_v[sU]()(0)(ic)),coalescedRead(Atilde_v[sU]()(0)(jc)) ) ) tmpType;
            tmpType tmpentry;
            zeroit(tmpentry);
            for(int s=0;s<Ls;s++){
              int sF = s+Ls*sU;
              for(int spn=0;spn<Ns;spn++){ //sum over spin
                auto bb = coalescedRead(Btilde_v[sF]()(spn)(ic) ); //color entry
                auto aa = coalescedRead(Atilde_v[sF]()(spn)(jc) );
                auto op = outerProduct(bb,aa);
                tmpentry = tmpentry + op;
            }
          }
          coalescedWrite(mat_v[sU](mu)()(ic,jc), tmpentry);
        }
      }

      });
      }
#else
    {
      const int Nsimd = SiteSpinor::Nsimd();
      autoView( Btilde_v , Btilde, AcceleratorRead);
      autoView( Atilde_v , Atilde, AcceleratorRead);
      accelerator_for(sss,mat.Grid()->oSites(),Nsimd,{
	  int sU=sss;
  	  typedef decltype(coalescedRead(mat_v[sU](mu)() )) ColorMatrixType;
  	  ColorMatrixType sum;
	  zeroit(sum);  
	  for(int s=0;s<Ls;s++){
	    int sF = s+Ls*sU;
  	    for(int spn=0;spn<Ns;spn++){ //sum over spin
  	      auto bb = coalescedRead(Btilde_v[sF]()(spn) ); //color vector
  	      auto aa = coalescedRead(Atilde_v[sF]()(spn) );
	      auto op = outerProduct(bb,aa);
  	      sum = sum + op;
	    }
	  }
  	  coalescedWrite(mat_v[sU](mu)(), sum);
      });
    }
#endif
#endif    
  }
};


typedef WilsonImpl<vComplex,  FundamentalRepresentation, CoeffReal > WilsonImplR;  // Real.. whichever prec
typedef WilsonImpl<vComplexF, FundamentalRepresentation, CoeffReal > WilsonImplF;  // Float
typedef WilsonImpl<vComplexD, FundamentalRepresentation, CoeffReal > WilsonImplD;  // Double

//typedef WilsonImpl<vComplex,  FundamentalRepresentation, CoeffRealHalfComms > WilsonImplRL;  // Real.. whichever prec
//typedef WilsonImpl<vComplexF, FundamentalRepresentation, CoeffRealHalfComms > WilsonImplFH;  // Float
//typedef WilsonImpl<vComplexD, FundamentalRepresentation, CoeffRealHalfComms > WilsonImplDF;  // Double

typedef WilsonImpl<vComplex,  FundamentalRepresentation, CoeffComplex > ZWilsonImplR; // Real.. whichever prec
typedef WilsonImpl<vComplexF, FundamentalRepresentation, CoeffComplex > ZWilsonImplF; // Float
typedef WilsonImpl<vComplexD, FundamentalRepresentation, CoeffComplex > ZWilsonImplD; // Double

//typedef WilsonImpl<vComplex,  FundamentalRepresentation, CoeffComplexHalfComms > ZWilsonImplRL; // Real.. whichever prec
//typedef WilsonImpl<vComplexF, FundamentalRepresentation, CoeffComplexHalfComms > ZWilsonImplFH; // Float
//typedef WilsonImpl<vComplexD, FundamentalRepresentation, CoeffComplexHalfComms > ZWilsonImplDF; // Double
 
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

