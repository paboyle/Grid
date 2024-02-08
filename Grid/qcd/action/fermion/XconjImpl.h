#pragma once

NAMESPACE_BEGIN(Grid);

  
/////////////////////////////////////////////////////////////////////////////
// Single flavour four spinors with colour index
/////////////////////////////////////////////////////////////////////////////
template <class S, class Representation = FundamentalRepresentation,class Options = CoeffReal >
class XconjugateWilsonImpl : public ConjugateGaugeImpl<GaugeImplTypes<S, Representation::Dimension> > {
public:

  static const int Dimension = Representation::Dimension;
  static const bool isFundamental = Representation::isFundamental;
  static const bool LsVectorised=false;
  static const bool isGparity=false;
  static const int Nhcs = Options::Nhcs;

  typedef ConjugateGaugeImpl< GaugeImplTypes<S,Dimension> > Gimpl;
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
  typedef XconjWilsonImplParams ImplParams;
  typedef WilsonStencil<SiteSpinor, SiteHalfSpinor,ImplParams> StencilImpl;
  typedef const typename StencilImpl::View_type StencilView;
    
  ImplParams Params;

  XconjugateWilsonImpl(const ImplParams &p = ImplParams()) : Params(p){
  };

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
    Gamma X = Gamma(Gamma::Algebra::MinusSigmaXZ);
    typedef decltype(coalescedRead( *( (SiteSpinor*)NULL ) )) _FullSpinor;

    int direction = St._directions[mu];
    int distance  = St._distances[mu];
    int ptype     = St._permute_type[mu];
    int sl        = St._simd_layout[direction];
    Coordinate icoor;

    const int Nsimd =SiteDoubledGaugeField::Nsimd();
    int s = acceleratorSIMTlane(Nsimd);
    St.iCoorFromIindex(icoor,s);

    int mmu = mu % Nd;
    int dag = mu / Nd;

    auto UU=coalescedRead(U(mu));
    
    //Decide whether we do a G-parity flavor twist
    //Note: this assumes (but does not check) that sl==1 || sl==2 i.e. max 2 SIMD lanes in G-parity dir
    //It also assumes (but does not check) that abs(distance) == 1
    int permute_lane = (sl==1) 
    || ((distance== 1)&&(icoor[direction]==1))
    || ((distance==-1)&&(icoor[direction]==0));

    permute_lane = permute_lane && SE->_around_the_world && St.parameters.twists[mmu] && mmu < Nd-1; //only if we are going around the world in a spatial direction

    //Implement the X-conjugate BC
    /*
      Note 
      C g^mu C = g^mu,T
      g^mu C = -C g^mu,T

      We need
      (1+g^mu) Cg^5 psi^*
      = g^5(1-g^mu) C psi^*
      = g^5 (C - g^mu C) psi^*
      = g^5 (C + C g^mu,T ) psi^*
      = g^5 C (1 + g^mu,T ) psi^*
      = [ g^5 C (1 + g^mu,dag ) psi ]^*
      = [ g^5 C (1 + g^mu ) psi ]^*

      We also need to return a half spinor but we can apply the projection again once we have
      (1+g^mu) Cg^5 psi^*

      Note:
      Applying Proj after Recon picks up an extra factor of 2
      e.g.  (from TwoSpinor.h)
      spProjXp =  [0]+i[3]
                  [1]+i[2]
		  
      spReconXp =    [0]
                     [1]
                   -i[1]
                   -i[0]

       spProjXp spReconXp =    [0] + [0]
                               [1] + [1]
    */
    if(permute_lane){
      _FullSpinor fs; 
      WilsonProjector::Recon(fs,chi,mmu,dag);
      fs = distance == 1 ? -0.5*(X*conjugate(fs)) : 0.5*(X*conjugate(fs)); //extra 0.5 due to Recon/Proj cycle
      fs = fs * St.parameters.boundary_phase;
      _Spinor XchiStar;
      WilsonProjector::Proj(XchiStar,fs,mmu,dag);
      mult(&phi(),&UU,&XchiStar());
    }else{
      mult(&phi(),&UU,&chi());
    }
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
    for (int mu = 0; mu <= Nd-1; mu++) { //mu=Nd-1 is either periodic or antiperiodic
      int L   = GaugeGrid->GlobalDimensions()[mu];
      int Lmu = L - 1;

      LatticeCoordinate(coor, mu);

      U = PeekIndex<LorentzIndex>(Umu, mu);
      
      if(mu == Nd-1 && Params.twists[mu] ){ //antiperiodic BCs
	tmp = where(coor == Lmu, -U, U);
	PokeIndex<LorentzIndex>(Uds, tmp, mu);

	U = adj(Cshift(U, mu, -1));
	U = where(coor == 0, -U, U); 
	PokeIndex<LorentzIndex>(Uds, U, mu + 4);
      }else if(mu < Nd-1 && Params.twists[mu] ){ //charge conj BCs
	PokeIndex<LorentzIndex>(Uds, U, mu);
	
	U = adj(Cshift(U, mu, -1));
	tmp = conjugate(U);
	U = where(coor == 0, tmp, U);

	PokeIndex<LorentzIndex>(Uds, U, mu + 4);
      }else{                                     //periodic BCs
	PokeIndex<LorentzIndex>(Uds, U, mu);

	U = adj(Cshift(U, mu, -1));
	PokeIndex<LorentzIndex>(Uds, U, mu + 4);
      }
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
      
    inline void extractLinkField(std::vector<GaugeLinkField> &mat, DoubledGaugeField &Uds)
    {
      for (int mu = 0; mu < Nd; mu++)
      mat[mu] = PeekIndex<LorentzIndex>(Uds, mu);
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
  }
};

typedef XconjugateWilsonImpl<vComplex,  FundamentalRepresentation, CoeffReal > XconjugateWilsonImplR;  // Real.. whichever prec
typedef XconjugateWilsonImpl<vComplexF, FundamentalRepresentation, CoeffReal > XconjugateWilsonImplF;  // Float
typedef XconjugateWilsonImpl<vComplexD, FundamentalRepresentation, CoeffReal > XconjugateWilsonImplD;  // Double

NAMESPACE_END(Grid);
