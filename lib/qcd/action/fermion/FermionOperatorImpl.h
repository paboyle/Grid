/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/FermionOperatorImpl.h

Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
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
#ifndef GRID_QCD_FERMION_OPERATOR_IMPL_H
#define GRID_QCD_FERMION_OPERATOR_IMPL_H

namespace Grid {
namespace QCD {

  //////////////////////////////////////////////
  // Template parameter class constructs to package
  // externally control Fermion implementations
  // in orthogonal directions
  //
  // Ultimately need Impl to always define types where XXX is opaque
  //
  //    typedef typename XXX               Simd;
  //    typedef typename XXX     GaugeLinkField;        
  //    typedef typename XXX         GaugeField;
  //    typedef typename XXX      GaugeActField;
  //    typedef typename XXX       FermionField;
  //    typedef typename XXX    PropagatorField;
  //    typedef typename XXX  DoubledGaugeField;
  //    typedef typename XXX         SiteSpinor;
  //    typedef typename XXX     SitePropagator;
  //    typedef typename XXX     SiteHalfSpinor;	
  //    typedef typename XXX         Compressor;	
  //
  // and Methods:
  //    void ImportGauge(GridBase *GaugeGrid,DoubledGaugeField &Uds,const GaugeField &Umu)
  //    void DoubleStore(GridBase *GaugeGrid,DoubledGaugeField &Uds,const GaugeField &Umu)
  //    void multLink(SiteHalfSpinor &phi,const SiteDoubledGaugeField &U,const SiteHalfSpinor &chi,int mu,StencilEntry *SE,StencilImpl &St)
  //    void InsertForce4D(GaugeField &mat,const FermionField &Btilde,const FermionField &A,int mu)
  //    void InsertForce5D(GaugeField &mat,const FermionField &Btilde,const FermionField &A,int mu)
  //
  //
  // To acquire the typedefs from "Base" (either a base class or template param) use:
  //
  // INHERIT_GIMPL_TYPES(Base)
  // INHERIT_FIMPL_TYPES(Base)
  // INHERIT_IMPL_TYPES(Base)
  //
  // The Fermion operators will do the following:
  //
  // struct MyOpParams { 
  //   RealD mass;
  // };
  //
  //
  // template<class Impl>
  // class MyOp : public<Impl> { 
  // public:
  //
  //    INHERIT_ALL_IMPL_TYPES(Impl);
  //
  //    MyOp(MyOpParams Myparm, ImplParams &ImplParam) :  Impl(ImplParam)
  //    {
  //
  //    };
  //    
  //  }
  //////////////////////////////////////////////

  template <class T> struct SamePrecisionMapper {
    typedef T HigherPrecVector ;
    typedef T LowerPrecVector ;
  };
  template <class T> struct LowerPrecisionMapper {  };
  template <> struct LowerPrecisionMapper<vRealF> {
    typedef vRealF HigherPrecVector ;
    typedef vRealH LowerPrecVector ;
  };
  template <> struct LowerPrecisionMapper<vRealD> {
    typedef vRealD HigherPrecVector ;
    typedef vRealF LowerPrecVector ;
  };
  template <> struct LowerPrecisionMapper<vComplexF> {
    typedef vComplexF HigherPrecVector ;
    typedef vComplexH LowerPrecVector ;
  };
  template <> struct LowerPrecisionMapper<vComplexD> {
    typedef vComplexD HigherPrecVector ;
    typedef vComplexF LowerPrecVector ;
  };

  struct CoeffReal {
  public:
    typedef RealD _Coeff_t;
    static const int Nhcs = 2;
    template<class Simd> using PrecisionMapper = SamePrecisionMapper<Simd>;
  };
  struct CoeffRealHalfComms {
  public:
    typedef RealD _Coeff_t;
    static const int Nhcs = 1;
    template<class Simd> using PrecisionMapper = LowerPrecisionMapper<Simd>;
  };
  struct CoeffComplex {
  public:
    typedef ComplexD _Coeff_t;
    static const int Nhcs = 2;
    template<class Simd> using PrecisionMapper = SamePrecisionMapper<Simd>;
  };
  struct CoeffComplexHalfComms {
  public:
    typedef ComplexD _Coeff_t;
    static const int Nhcs = 1;
    template<class Simd> using PrecisionMapper = LowerPrecisionMapper<Simd>;
  };

  ////////////////////////////////////////////////////////////////////////
  // Implementation dependent fermion types
  ////////////////////////////////////////////////////////////////////////
  
#define INHERIT_FIMPL_TYPES(Impl)\
  typedef typename Impl::FermionField           FermionField;		\
  typedef typename Impl::PropagatorField     PropagatorField;		\
  typedef typename Impl::DoubledGaugeField DoubledGaugeField;		\
  typedef typename Impl::SiteSpinor               SiteSpinor;		\
  typedef typename Impl::SitePropagator       SitePropagator;		\
  typedef typename Impl::SiteHalfSpinor       SiteHalfSpinor;		\
  typedef typename Impl::Compressor               Compressor;		\
  typedef typename Impl::StencilImpl             StencilImpl;		\
  typedef typename Impl::ImplParams               ImplParams;	        \
  typedef typename Impl::Coeff_t                     Coeff_t;           \
  
#define INHERIT_IMPL_TYPES(Base) \
  INHERIT_GIMPL_TYPES(Base)      \
  INHERIT_FIMPL_TYPES(Base)
  
  /////////////////////////////////////////////////////////////////////////////
  // Single flavour four spinors with colour index
  /////////////////////////////////////////////////////////////////////////////
  template <class S, class Representation = FundamentalRepresentation,class Options = CoeffReal >
  class WilsonImpl : public PeriodicGaugeImpl<GaugeImplTypes<S, Representation::Dimension > > {
    public:

    static const int Dimension = Representation::Dimension;
    static const bool LsVectorised=false;
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
    typedef WilsonStencil<SiteSpinor, SiteHalfSpinor> StencilImpl;
    
    ImplParams Params;
    
    WilsonImpl(const ImplParams &p = ImplParams()) : Params(p){
      assert(Params.boundary_phases.size() == Nd);
    };
      
    bool overlapCommsCompute(void) { return Params.overlapCommsCompute; };
      
    inline void multLink(SiteHalfSpinor &phi,
                         const SiteDoubledGaugeField &U,
                         const SiteHalfSpinor &chi,
                         int mu,
                         StencilEntry *SE,
                         StencilImpl &St) {
      mult(&phi(), &U(mu), &chi());
    }
      
    template <class ref>
    inline void loadLinkElement(Simd &reg, ref &memory) {
      reg = memory;
    }
      
    inline void DoubleStore(GridBase *GaugeGrid,
                            DoubledGaugeField &Uds,
                            const GaugeField &Umu) 
    {
      typedef typename Simd::scalar_type scalar_type;

      conformable(Uds._grid, GaugeGrid);
      conformable(Umu._grid, GaugeGrid);

      GaugeLinkField U(GaugeGrid);
      GaugeLinkField tmp(GaugeGrid);

      Lattice<iScalar<vInteger> > coor(GaugeGrid);
      for (int mu = 0; mu < Nd; mu++) {

	      auto pha = Params.boundary_phases[mu];
	      scalar_type phase( real(pha),imag(pha) );

        int Lmu = GaugeGrid->GlobalDimensions()[mu] - 1;

        LatticeCoordinate(coor, mu);

        U = PeekIndex<LorentzIndex>(Umu, mu);
        tmp = where(coor == Lmu, phase * U, U);
        PokeIndex<LorentzIndex>(Uds, tmp, mu);

        U = adj(Cshift(U, mu, -1));
        U = where(coor == 0, conjugate(phase) * U, U); 
        PokeIndex<LorentzIndex>(Uds, U, mu + 4);
      }
    }

    inline void InsertForce4D(GaugeField &mat, FermionField &Btilde, FermionField &A,int mu){
      GaugeLinkField link(mat._grid);
      link = TraceIndex<SpinIndex>(outerProduct(Btilde,A)); 
      PokeIndex<LorentzIndex>(mat,link,mu);
    }   
      
    inline void InsertForce5D(GaugeField &mat, FermionField &Btilde, FermionField &Atilde,int mu){
      
      int Ls=Btilde._grid->_fdimensions[0];
      GaugeLinkField tmp(mat._grid);
      tmp = zero;
      
      parallel_for(int sss=0;sss<tmp._grid->oSites();sss++){
	    int sU=sss;
	    for(int s=0;s<Ls;s++){
	        int sF = s+Ls*sU;
	        tmp[sU] = tmp[sU]+ traceIndex<SpinIndex>(outerProduct(Btilde[sF],Atilde[sF])); // ordering here
	    }
      }
      PokeIndex<LorentzIndex>(mat,tmp,mu);
      
    }
  };

  ////////////////////////////////////////////////////////////////////////////////////
  // Single flavour four spinors with colour index, 5d redblack
  ////////////////////////////////////////////////////////////////////////////////////
template<class S,int Nrepresentation=Nc, class Options=CoeffReal>
class DomainWallVec5dImpl :  public PeriodicGaugeImpl< GaugeImplTypes< S,Nrepresentation> > { 
  public:

  typedef PeriodicGaugeImpl<GaugeImplTypes<S, Nrepresentation> > Gimpl;
  INHERIT_GIMPL_TYPES(Gimpl);

  static const int Dimension = Nrepresentation;
  static const bool LsVectorised=true;
  static const int Nhcs = Options::Nhcs;
      
  typedef typename Options::_Coeff_t Coeff_t;      
  typedef typename Options::template PrecisionMapper<Simd>::LowerPrecVector SimdL;
  
  template <typename vtype> using iImplSpinor            = iScalar<iVector<iVector<vtype, Nrepresentation>, Ns> >;
  template <typename vtype> using iImplPropagator        = iScalar<iMatrix<iMatrix<vtype, Nrepresentation>, Ns> >;
  template <typename vtype> using iImplHalfSpinor        = iScalar<iVector<iVector<vtype, Nrepresentation>, Nhs> >;
  template <typename vtype> using iImplHalfCommSpinor    = iScalar<iVector<iVector<vtype, Nrepresentation>, Nhcs> >;
  template <typename vtype> using iImplDoubledGaugeField = iVector<iScalar<iMatrix<vtype, Nrepresentation> >, Nds>;
  template <typename vtype> using iImplGaugeField        = iVector<iScalar<iMatrix<vtype, Nrepresentation> >, Nd>;
  template <typename vtype> using iImplGaugeLink         = iScalar<iScalar<iMatrix<vtype, Nrepresentation> > >;
  
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
  typedef WilsonStencil<SiteSpinor, SiteHalfSpinor> StencilImpl;
  
  ImplParams Params;
  
  DomainWallVec5dImpl(const ImplParams &p = ImplParams()) : Params(p){};
      
  bool overlapCommsCompute(void) { return false; };
      
  template <class ref>
  inline void loadLinkElement(Simd &reg, ref &memory) {
    vsplat(reg, memory);
  }

  inline void multLink(SiteHalfSpinor &phi, const SiteDoubledGaugeField &U,
                       const SiteHalfSpinor &chi, int mu, StencilEntry *SE,
                       StencilImpl &St) {
    SiteGaugeLink UU;
    for (int i = 0; i < Nrepresentation; i++) {
      for (int j = 0; j < Nrepresentation; j++) {
        vsplat(UU()()(i, j), U(mu)()(i, j));
      }
    }
    mult(&phi(), &UU(), &chi());
  }
      
  inline void DoubleStore(GridBase *GaugeGrid, DoubledGaugeField &Uds,const GaugeField &Umu) 
  {
    SiteScalarGaugeField  ScalarUmu;
    SiteDoubledGaugeField ScalarUds;
    
    GaugeLinkField U(Umu._grid);
    GaugeField  Uadj(Umu._grid);
    for (int mu = 0; mu < Nd; mu++) {
      U = PeekIndex<LorentzIndex>(Umu, mu);
      U = adj(Cshift(U, mu, -1));
      PokeIndex<LorentzIndex>(Uadj, U, mu);
    }
    
    for (int lidx = 0; lidx < GaugeGrid->lSites(); lidx++) {
      std::vector<int> lcoor;
      GaugeGrid->LocalIndexToLocalCoor(lidx, lcoor);
      
      peekLocalSite(ScalarUmu, Umu, lcoor);
      for (int mu = 0; mu < 4; mu++) ScalarUds(mu) = ScalarUmu(mu);
      
      peekLocalSite(ScalarUmu, Uadj, lcoor);
      for (int mu = 0; mu < 4; mu++) ScalarUds(mu + 4) = ScalarUmu(mu);
      
      pokeLocalSite(ScalarUds, Uds, lcoor);
    }
  }
      
  inline void InsertForce4D(GaugeField &mat, FermionField &Btilde,FermionField &A, int mu) 
  {
    assert(0);
  }

  inline void InsertForce5D(GaugeField &mat, FermionField &Btilde, FermionField &Atilde, int mu) {

    assert(0);
    // Following lines to be revised after Peter's addition of half prec
    // missing put lane...
    /*
    typedef decltype(traceIndex<SpinIndex>(outerProduct(Btilde[0], Atilde[0]))) result_type;
    unsigned int LLs = Btilde._grid->_rdimensions[0];
    conformable(Atilde._grid,Btilde._grid);
    GridBase* grid = mat._grid;
    GridBase* Bgrid = Btilde._grid;
    unsigned int dimU = grid->Nd();
    unsigned int dimF = Bgrid->Nd();
    GaugeLinkField tmp(grid); 
    tmp = zero;
    
    // FIXME 
    // Current implementation works, thread safe, probably suboptimal
    // Passing through the local coordinate for grid transformation
    // the force grid is in general very different from the Ls vectorized grid

    PARALLEL_FOR_LOOP
    for (int so = 0; so < grid->oSites(); so++) {
      std::vector<typename result_type::scalar_object> vres(Bgrid->Nsimd());
      std::vector<int> ocoor;  grid->oCoorFromOindex(ocoor,so); 
      for (int si = 0; si < tmp._grid->iSites(); si++){
        typename result_type::scalar_object scalar_object; scalar_object = zero;
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
        tmp._odata[so].putlane(scalar_object, si);
      }
    }
    PokeIndex<LorentzIndex>(mat, tmp, mu);
    */
  }
};
    
    ////////////////////////////////////////////////////////////////////////////////////////
    // Flavour doubled spinors; is Gparity the only? what about C*?
    ////////////////////////////////////////////////////////////////////////////////////////
template <class S, int Nrepresentation, class Options=CoeffReal>
class GparityWilsonImpl : public ConjugateGaugeImpl<GaugeImplTypes<S, Nrepresentation> > {
 public:

 static const int Dimension = Nrepresentation;
 static const int Nhcs = Options::Nhcs;
 static const bool LsVectorised=false;

 typedef ConjugateGaugeImpl< GaugeImplTypes<S,Nrepresentation> > Gimpl;
 INHERIT_GIMPL_TYPES(Gimpl);

 typedef typename Options::_Coeff_t Coeff_t;
 typedef typename Options::template PrecisionMapper<Simd>::LowerPrecVector SimdL;
      
 template <typename vtype> using iImplSpinor            = iVector<iVector<iVector<vtype, Nrepresentation>, Ns>,   Ngp>;
 template <typename vtype> using iImplPropagator        = iVector<iMatrix<iMatrix<vtype, Nrepresentation>, Ns>,   Ngp>;
 template <typename vtype> using iImplHalfSpinor        = iVector<iVector<iVector<vtype, Nrepresentation>, Nhs>,  Ngp>;
 template <typename vtype> using iImplHalfCommSpinor    = iVector<iVector<iVector<vtype, Nrepresentation>, Nhcs>, Ngp>;
 template <typename vtype> using iImplDoubledGaugeField = iVector<iVector<iScalar<iMatrix<vtype, Nrepresentation> >, Nds>, Ngp>;

 typedef iImplSpinor<Simd>            SiteSpinor;
 typedef iImplPropagator<Simd>        SitePropagator;
 typedef iImplHalfSpinor<Simd>        SiteHalfSpinor;
 typedef iImplHalfCommSpinor<SimdL>   SiteHalfCommSpinor;
 typedef iImplDoubledGaugeField<Simd> SiteDoubledGaugeField;

 typedef Lattice<SiteSpinor> FermionField;
 typedef Lattice<SitePropagator> PropagatorField;
 typedef Lattice<SiteDoubledGaugeField> DoubledGaugeField;
 
 typedef WilsonCompressor<SiteHalfCommSpinor,SiteHalfSpinor, SiteSpinor> Compressor;
 typedef WilsonStencil<SiteSpinor, SiteHalfSpinor> StencilImpl;
 
 typedef GparityWilsonImplParams ImplParams;
      
 ImplParams Params;

 GparityWilsonImpl(const ImplParams &p = ImplParams()) : Params(p){};

 bool overlapCommsCompute(void) { return Params.overlapCommsCompute; };

 // provide the multiply by link that is differentiated between Gparity (with
 // flavour index) and non-Gparity
 inline void multLink(SiteHalfSpinor &phi, const SiteDoubledGaugeField &U,
                      const SiteHalfSpinor &chi, int mu, StencilEntry *SE,
                      StencilImpl &St) {

   typedef SiteHalfSpinor vobj;
   typedef typename SiteHalfSpinor::scalar_object sobj;
	
   vobj vtmp;
   sobj stmp;
        
   GridBase *grid = St._grid;
        
   const int Nsimd = grid->Nsimd();
        
   int direction = St._directions[mu];
   int distance = St._distances[mu];
   int ptype = St._permute_type[mu];
   int sl = St._grid->_simd_layout[direction];
   
   // Fixme X.Y.Z.T hardcode in stencil
   int mmu = mu % Nd;
        
   // assert our assumptions
   assert((distance == 1) || (distance == -1));  // nearest neighbour stencil hard code
   assert((sl == 1) || (sl == 2));
   
   std::vector<int> icoor;
        
   if ( SE->_around_the_world && Params.twists[mmu] ) {

     if ( sl == 2 ) {
       
       std::vector<sobj> vals(Nsimd);

       extract(chi,vals);
       for(int s=0;s<Nsimd;s++){

         grid->iCoorFromIindex(icoor,s);
              
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
   
 }


 template <class ref>
 inline void loadLinkElement(Simd &reg, ref &memory) {
   reg = memory;
 }

 inline void DoubleStore(GridBase *GaugeGrid,DoubledGaugeField &Uds,const GaugeField &Umu)
 {
   conformable(Uds._grid,GaugeGrid);
   conformable(Umu._grid,GaugeGrid);
   
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
	  
     parallel_for(auto ss=U.begin();ss<U.end();ss++){
       Uds[ss](0)(mu) = U[ss]();
       Uds[ss](1)(mu) = Uconj[ss]();
     }
          
     U     = adj(Cshift(U    ,mu,-1));      // correct except for spanning the boundary
     Uconj = adj(Cshift(Uconj,mu,-1));
 
     Utmp = U;
     if ( Params.twists[mu] ) { 
       Utmp = where(coor==0,Uconj,Utmp);
     }

	  
     parallel_for(auto ss=U.begin();ss<U.end();ss++){
       Uds[ss](0)(mu+4) = Utmp[ss]();
     }
          
     Utmp = Uconj;
     if ( Params.twists[mu] ) { 
       Utmp = where(coor==0,U,Utmp);
     }
	  
     parallel_for(auto ss=U.begin();ss<U.end();ss++){
       Uds[ss](1)(mu+4) = Utmp[ss]();
     }
          
   }
 }
      
 inline void InsertForce4D(GaugeField &mat, FermionField &Btilde, FermionField &A, int mu) {

   // DhopDir provides U or Uconj depending on coor/flavour.
   GaugeLinkField link(mat._grid);
   // use lorentz for flavour as hack.
   auto tmp = TraceIndex<SpinIndex>(outerProduct(Btilde, A));
   parallel_for(auto ss = tmp.begin(); ss < tmp.end(); ss++) {
     link[ss]() = tmp[ss](0, 0) + conjugate(tmp[ss](1, 1));
   }
   PokeIndex<LorentzIndex>(mat, link, mu);
   return;
 }
      
 inline void InsertForce5D(GaugeField &mat, FermionField &Btilde, FermionField &Atilde, int mu) {

   int Ls = Btilde._grid->_fdimensions[0];
        
   GaugeLinkField tmp(mat._grid);
   tmp = zero;
   parallel_for(int ss = 0; ss < tmp._grid->oSites(); ss++) {
     for (int s = 0; s < Ls; s++) {
       int sF = s + Ls * ss;
       auto ttmp = traceIndex<SpinIndex>(outerProduct(Btilde[sF], Atilde[sF]));
       tmp[ss]() = tmp[ss]() + ttmp(0, 0) + conjugate(ttmp(1, 1));
     }
   }
   PokeIndex<LorentzIndex>(mat, tmp, mu);
   return;
 }

};

/////////////////////////////////////////////////////////////////////////////
// Single flavour one component spinors with colour index
/////////////////////////////////////////////////////////////////////////////
template <class S, class Representation = FundamentalRepresentation >
class StaggeredImpl : public PeriodicGaugeImpl<GaugeImplTypes<S, Representation::Dimension > > {

    public:

    typedef RealD  _Coeff_t ;
    static const int Dimension = Representation::Dimension;
    static const bool LsVectorised=false;
    typedef PeriodicGaugeImpl<GaugeImplTypes<S, Dimension > > Gimpl;
      
    //Necessary?
    constexpr bool is_fundamental() const{return Dimension == Nc ? 1 : 0;}
    
    typedef _Coeff_t Coeff_t;

    INHERIT_GIMPL_TYPES(Gimpl);
      
    template <typename vtype> using iImplSpinor            = iScalar<iScalar<iVector<vtype, Dimension> > >;
    template <typename vtype> using iImplHalfSpinor        = iScalar<iScalar<iVector<vtype, Dimension> > >;
    template <typename vtype> using iImplDoubledGaugeField = iVector<iScalar<iMatrix<vtype, Dimension> >, Nds>;
    template <typename vtype> using iImplPropagator        = iScalar<iScalar<iMatrix<vtype, Dimension> > >;
    
    typedef iImplSpinor<Simd>            SiteSpinor;
    typedef iImplHalfSpinor<Simd>        SiteHalfSpinor;
    typedef iImplDoubledGaugeField<Simd> SiteDoubledGaugeField;
    typedef iImplPropagator<Simd>        SitePropagator;
    
    typedef Lattice<SiteSpinor>            FermionField;
    typedef Lattice<SiteDoubledGaugeField> DoubledGaugeField;
    typedef Lattice<SitePropagator> PropagatorField;
    
    typedef SimpleCompressor<SiteSpinor> Compressor;
    typedef StaggeredImplParams ImplParams;
    typedef CartesianStencil<SiteSpinor, SiteSpinor> StencilImpl;
    
    ImplParams Params;
    
    StaggeredImpl(const ImplParams &p = ImplParams()) : Params(p){};
      
    inline void multLink(SiteSpinor &phi,
			 const SiteDoubledGaugeField &U,
			 const SiteSpinor &chi,
			 int mu){
      mult(&phi(), &U(mu), &chi());
    }
    inline void multLinkAdd(SiteSpinor &phi,
			    const SiteDoubledGaugeField &U,
			    const SiteSpinor &chi,
			    int mu){
      mac(&phi(), &U(mu), &chi());
    }
      
    template <class ref>
    inline void loadLinkElement(Simd &reg, ref &memory) {
      reg = memory;
    }
      
    inline void DoubleStore(GridBase *GaugeGrid,
			    DoubledGaugeField &UUUds, // for Naik term
			    DoubledGaugeField &Uds,
			    const GaugeField &Uthin,
			    const GaugeField &Ufat) {
      conformable(Uds._grid, GaugeGrid);
      conformable(Uthin._grid, GaugeGrid);
      conformable(Ufat._grid, GaugeGrid);
      GaugeLinkField U(GaugeGrid);
      GaugeLinkField UU(GaugeGrid);
      GaugeLinkField UUU(GaugeGrid);
      GaugeLinkField Udag(GaugeGrid);
      GaugeLinkField UUUdag(GaugeGrid);
      for (int mu = 0; mu < Nd; mu++) {

	// Staggered Phase.
	Lattice<iScalar<vInteger> > coor(GaugeGrid);
	Lattice<iScalar<vInteger> > x(GaugeGrid); LatticeCoordinate(x,0);
	Lattice<iScalar<vInteger> > y(GaugeGrid); LatticeCoordinate(y,1);
	Lattice<iScalar<vInteger> > z(GaugeGrid); LatticeCoordinate(z,2);
	Lattice<iScalar<vInteger> > t(GaugeGrid); LatticeCoordinate(t,3);

	Lattice<iScalar<vInteger> > lin_z(GaugeGrid); lin_z=x+y;
	Lattice<iScalar<vInteger> > lin_t(GaugeGrid); lin_t=x+y+z;

	ComplexField phases(GaugeGrid);	phases=1.0;

	if ( mu == 1 ) phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
	if ( mu == 2 ) phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
	if ( mu == 3 ) phases = where( mod(lin_t,2)==(Integer)0, phases,-phases);

	// 1 hop based on fat links
	U      = PeekIndex<LorentzIndex>(Ufat, mu);
	Udag   = adj( Cshift(U, mu, -1));

	U    = U    *phases;
	Udag = Udag *phases;

	PokeIndex<LorentzIndex>(Uds, U, mu);
	PokeIndex<LorentzIndex>(Uds, Udag, mu + 4);

	// 3 hop based on thin links. Crazy huh ?
	U  = PeekIndex<LorentzIndex>(Uthin, mu);
	UU = Gimpl::CovShiftForward(U,mu,U);
	UUU= Gimpl::CovShiftForward(U,mu,UU);
	
	UUUdag = adj( Cshift(UUU, mu, -3));

	UUU    = UUU    *phases;
	UUUdag = UUUdag *phases;

	PokeIndex<LorentzIndex>(UUUds, UUU, mu);
	PokeIndex<LorentzIndex>(UUUds, UUUdag, mu+4);

      }
    }

    inline void InsertForce4D(GaugeField &mat, FermionField &Btilde, FermionField &A,int mu){
      GaugeLinkField link(mat._grid);
      link = TraceIndex<SpinIndex>(outerProduct(Btilde,A)); 
      PokeIndex<LorentzIndex>(mat,link,mu);
    }   
      
    inline void InsertForce5D(GaugeField &mat, FermionField &Btilde, FermionField &Atilde,int mu){
      assert (0); 
      // Must never hit
    }
  };

  /////////////////////////////////////////////////////////////////////////////
  // Single flavour one component spinors with colour index. 5d vec
  /////////////////////////////////////////////////////////////////////////////
  template <class S, class Representation = FundamentalRepresentation >
  class StaggeredVec5dImpl : public PeriodicGaugeImpl<GaugeImplTypes<S, Representation::Dimension > > {

    public:

    static const int Dimension = Representation::Dimension;
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
    
    typedef SimpleCompressor<SiteSpinor> Compressor;
    typedef StaggeredImplParams ImplParams;
    typedef CartesianStencil<SiteSpinor, SiteSpinor> StencilImpl;
    
    ImplParams Params;
    
    StaggeredVec5dImpl(const ImplParams &p = ImplParams()) : Params(p){};

    template <class ref>
    inline void loadLinkElement(Simd &reg, ref &memory) {
      vsplat(reg, memory);
    }

    inline void multLink(SiteHalfSpinor &phi, const SiteDoubledGaugeField &U,
			 const SiteHalfSpinor &chi, int mu) {
      SiteGaugeLink UU;
      for (int i = 0; i < Dimension; i++) {
	for (int j = 0; j < Dimension; j++) {
	  vsplat(UU()()(i, j), U(mu)()(i, j));
	}
      }
      mult(&phi(), &UU(), &chi());
    }
    inline void multLinkAdd(SiteHalfSpinor &phi, const SiteDoubledGaugeField &U,
			    const SiteHalfSpinor &chi, int mu) {
      SiteGaugeLink UU;
      for (int i = 0; i < Dimension; i++) {
	for (int j = 0; j < Dimension; j++) {
	  vsplat(UU()()(i, j), U(mu)()(i, j));
	}
      }
      mac(&phi(), &UU(), &chi());
    }
      
    inline void DoubleStore(GridBase *GaugeGrid,
			    DoubledGaugeField &UUUds, // for Naik term
			    DoubledGaugeField &Uds,
			    const GaugeField &Uthin,
			    const GaugeField &Ufat) 
    {

      GridBase * InputGrid = Uthin._grid;
      conformable(InputGrid,Ufat._grid);

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


	for (int lidx = 0; lidx < GaugeGrid->lSites(); lidx++) {
	  SiteScalarGaugeLink   ScalarU;
	  SiteDoubledGaugeField ScalarUds;
	  
	  std::vector<int> lcoor;
	  GaugeGrid->LocalIndexToLocalCoor(lidx, lcoor);
	  peekLocalSite(ScalarUds, Uds, lcoor);

	  peekLocalSite(ScalarU, U, lcoor);
	  ScalarUds(mu) = ScalarU();

	  peekLocalSite(ScalarU, Udag, lcoor);
	  ScalarUds(mu + 4) = ScalarU();

	  pokeLocalSite(ScalarUds, Uds, lcoor);
	}

	// 3 hop based on thin links. Crazy huh ?
	U  = PeekIndex<LorentzIndex>(Uthin, mu);
	UU = Gimpl::CovShiftForward(U,mu,U);
	UUU= Gimpl::CovShiftForward(U,mu,UU);
	
	UUUdag = adj( Cshift(UUU, mu, -3));

	UUU    = UUU    *phases;
	UUUdag = UUUdag *phases;

	for (int lidx = 0; lidx < GaugeGrid->lSites(); lidx++) {

	  SiteScalarGaugeLink  ScalarU;
	  SiteDoubledGaugeField ScalarUds;
	  
	  std::vector<int> lcoor;
	  GaugeGrid->LocalIndexToLocalCoor(lidx, lcoor);
      
	  peekLocalSite(ScalarUds, UUUds, lcoor);

	  peekLocalSite(ScalarU, UUU, lcoor);
	  ScalarUds(mu) = ScalarU();

	  peekLocalSite(ScalarU, UUUdag, lcoor);
	  ScalarUds(mu + 4) = ScalarU();
	  
	  pokeLocalSite(ScalarUds, UUUds, lcoor);
	}

      }
    }

    inline void InsertForce4D(GaugeField &mat, FermionField &Btilde, FermionField &A,int mu){
      assert(0);
    }   
      
    inline void InsertForce5D(GaugeField &mat, FermionField &Btilde, FermionField &Atilde,int mu){
      assert (0); 
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
 
typedef DomainWallVec5dImpl<vComplex ,Nc, CoeffReal> DomainWallVec5dImplR; // Real.. whichever prec
typedef DomainWallVec5dImpl<vComplexF,Nc, CoeffReal> DomainWallVec5dImplF; // Float
typedef DomainWallVec5dImpl<vComplexD,Nc, CoeffReal> DomainWallVec5dImplD; // Double
 
typedef DomainWallVec5dImpl<vComplex ,Nc, CoeffRealHalfComms> DomainWallVec5dImplRL; // Real.. whichever prec
typedef DomainWallVec5dImpl<vComplexF,Nc, CoeffRealHalfComms> DomainWallVec5dImplFH; // Float
typedef DomainWallVec5dImpl<vComplexD,Nc, CoeffRealHalfComms> DomainWallVec5dImplDF; // Double
 
typedef DomainWallVec5dImpl<vComplex ,Nc,CoeffComplex> ZDomainWallVec5dImplR; // Real.. whichever prec
typedef DomainWallVec5dImpl<vComplexF,Nc,CoeffComplex> ZDomainWallVec5dImplF; // Float
typedef DomainWallVec5dImpl<vComplexD,Nc,CoeffComplex> ZDomainWallVec5dImplD; // Double
 
typedef DomainWallVec5dImpl<vComplex ,Nc,CoeffComplexHalfComms> ZDomainWallVec5dImplRL; // Real.. whichever prec
typedef DomainWallVec5dImpl<vComplexF,Nc,CoeffComplexHalfComms> ZDomainWallVec5dImplFH; // Float
typedef DomainWallVec5dImpl<vComplexD,Nc,CoeffComplexHalfComms> ZDomainWallVec5dImplDF; // Double
 
typedef GparityWilsonImpl<vComplex , Nc,CoeffReal> GparityWilsonImplR;  // Real.. whichever prec
typedef GparityWilsonImpl<vComplexF, Nc,CoeffReal> GparityWilsonImplF;  // Float
typedef GparityWilsonImpl<vComplexD, Nc,CoeffReal> GparityWilsonImplD;  // Double
 
typedef GparityWilsonImpl<vComplex , Nc,CoeffRealHalfComms> GparityWilsonImplRL;  // Real.. whichever prec
typedef GparityWilsonImpl<vComplexF, Nc,CoeffRealHalfComms> GparityWilsonImplFH;  // Float
typedef GparityWilsonImpl<vComplexD, Nc,CoeffRealHalfComms> GparityWilsonImplDF;  // Double

typedef StaggeredImpl<vComplex,  FundamentalRepresentation > StaggeredImplR;   // Real.. whichever prec
typedef StaggeredImpl<vComplexF, FundamentalRepresentation > StaggeredImplF;  // Float
typedef StaggeredImpl<vComplexD, FundamentalRepresentation > StaggeredImplD;  // Double

typedef StaggeredVec5dImpl<vComplex,  FundamentalRepresentation > StaggeredVec5dImplR;   // Real.. whichever prec
typedef StaggeredVec5dImpl<vComplexF, FundamentalRepresentation > StaggeredVec5dImplF;  // Float
typedef StaggeredVec5dImpl<vComplexD, FundamentalRepresentation > StaggeredVec5dImplD;  // Double

}}

#endif
