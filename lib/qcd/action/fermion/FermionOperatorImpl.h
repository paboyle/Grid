#ifndef  GRID_QCD_FERMION_OPERATOR_IMPL_H
#define  GRID_QCD_FERMION_OPERATOR_IMPL_H

namespace Grid {

  namespace QCD {

    // Variable precision "S" and variable Nc
    template<class S,int Nrepresentation=Nc>
    class WilsonImpl { 
    public:

      typedef S Simd;

      template<typename vtype> using iImplSpinor             = iScalar<iVector<iVector<vtype, Nrepresentation>, Ns> >;
      template<typename vtype> using iImplHalfSpinor         = iScalar<iVector<iVector<vtype, Nrepresentation>, Nhs> >;
      template<typename vtype> using iImplGaugeLink          = iScalar<iScalar<iMatrix<vtype, Nrepresentation> > >;
      template<typename vtype> using iImplGaugeField         = iVector<iScalar<iMatrix<vtype, Nrepresentation> >, Nd  >;
      template<typename vtype> using iImplDoubledGaugeField  = iVector<iScalar<iMatrix<vtype, Nrepresentation> >, Nds >;
    
      typedef iImplSpinor    <Simd>           SiteSpinor;
      typedef iImplHalfSpinor<Simd>           SiteHalfSpinor;
      typedef iImplGaugeLink <Simd>           SiteGaugeLink;
      typedef iImplGaugeField<Simd>           SiteGaugeField;
      typedef iImplDoubledGaugeField<Simd>    SiteDoubledGaugeField;

      typedef Lattice<SiteSpinor>                 FermionField;
      typedef Lattice<SiteGaugeLink>            GaugeLinkField; // bit ugly naming; polarised gauge field, lorentz... all ugly
      typedef Lattice<SiteGaugeField>               GaugeField;
      typedef Lattice<SiteDoubledGaugeField> DoubledGaugeField;

      typedef WilsonCompressor<SiteHalfSpinor,SiteSpinor> Compressor;

      static inline void multLink(SiteHalfSpinor &phi,const SiteDoubledGaugeField &U,const SiteHalfSpinor &chi,int mu,StencilEntry *SE,CartesianStencil &St){
        mult(&phi(),&U(mu),&chi());
      }
      inline void DoubleStore(GridBase *GaugeGrid,DoubledGaugeField &Uds,const GaugeField &Umu)
      {
        conformable(Uds._grid,GaugeGrid);
        conformable(Umu._grid,GaugeGrid);
        GaugeLinkField U(GaugeGrid);
        for(int mu=0;mu<Nd;mu++){
  	  U = PeekIndex<LorentzIndex>(Umu,mu);
	  PokeIndex<LorentzIndex>(Uds,U,mu);
	  U = adj(Cshift(U,mu,-1));
	  PokeIndex<LorentzIndex>(Uds,U,mu+4);
	}
      }

      inline void InsertForce(GaugeField &mat,const FermionField &Btilde,const FermionField &A,int mu){
	GaugeLinkField link(mat._grid);
	link = TraceIndex<SpinIndex>(outerProduct(Btilde,A)); 
	PokeIndex<LorentzIndex>(mat,link,mu);
      }   

    };

    template<class S,int Nrepresentation=Nc>
    class GparityWilsonImpl { 
    public:

      typedef S Simd;

      template<typename vtype> using iImplSpinor             = iVector<iVector<iVector<vtype, Nrepresentation>, Ns>, Ngp >;
      template<typename vtype> using iImplHalfSpinor         = iVector<iVector<iVector<vtype, Nrepresentation>, Nhs>, Ngp >;
      template<typename vtype> using iImplGaugeField         = iVector<iScalar<iMatrix<vtype, Nrepresentation> >, Nd  >;

      template<typename vtype> using iImplGaugeLink          = iScalar<iScalar<iMatrix<vtype, Nrepresentation> > >;
      template<typename vtype> using iImplDoubledGaugeField  = iVector<iVector<iScalar<iMatrix<vtype, Nrepresentation> >, Nds >, Ngp >;
    
      typedef iImplSpinor    <Simd>           SiteSpinor;
      typedef iImplHalfSpinor<Simd>           SiteHalfSpinor;
      typedef iImplGaugeLink <Simd>           SiteGaugeLink;
      typedef iImplGaugeField<Simd>           SiteGaugeField;
      typedef iImplDoubledGaugeField<Simd>    SiteDoubledGaugeField;

      typedef Lattice<SiteSpinor>                 FermionField;
      typedef Lattice<SiteGaugeLink>            GaugeLinkField; // bit ugly naming; polarised gauge field, lorentz... all ugly
      typedef Lattice<SiteGaugeField>               GaugeField;
      typedef Lattice<SiteDoubledGaugeField> DoubledGaugeField;

      //      typedef GparityWilsonCompressor<SiteHalfSpinor,SiteSpinor> Compressor;
      typedef WilsonCompressor<SiteHalfSpinor,SiteSpinor> Compressor;

      // provide the multiply by link that is differentiated between Gparity (with flavour index) and 
      // non-Gparity
      static inline void multLink(SiteHalfSpinor &phi,const SiteDoubledGaugeField &U,const SiteHalfSpinor &chi,int mu,StencilEntry *SE,CartesianStencil &St){

	typedef SiteHalfSpinor vobj;
	typedef typename SiteHalfSpinor::scalar_object sobj;

	vobj vtmp;
	sobj stmp;
	std::vector<int> gpbc({0,0,0,1,0,0,0,1});
	
	GridBase *grid = St._grid;
      
	const int Nsimd = grid->Nsimd();
	
	int direction    = St._directions[mu];
	int distance     = St._distances[mu];
	int ptype     = St._permute_type[mu]; 
	int sl        = St._grid->_simd_layout[direction];
	
	// assert our assumptions
	assert((distance==1)||(distance==-1)); // nearest neighbour stencil hard code
	assert((sl==1)||(sl==2));
	
	std::vector<int> icoor;
      
	if ( SE->_around_the_world && gpbc[mu] ) {

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

    static inline void InsertForce(GaugeField &mat,const FermionField &Btilde,const FermionField &A,int mu){
      // Fixme
      return;
    }
    static inline void DoubleStore(GridBase *GaugeGrid,DoubledGaugeField &Uds,const GaugeField &Umu)
    {

      conformable(Uds._grid,GaugeGrid);
      conformable(Umu._grid,GaugeGrid);

      GaugeLinkField Utmp(GaugeGrid);
      GaugeLinkField U(GaugeGrid);
      GaugeLinkField Uconj(GaugeGrid);

      Lattice<iScalar<vInteger> > coor(GaugeGrid);
      
      std::vector<int> gpdirs({0,0,0,1});
      
      for(int mu=0;mu<Nd;mu++){
	
	LatticeCoordinate(coor,mu);
	
	U     = PeekIndex<LorentzIndex>(Umu,mu);
	Uconj = conjugate(U);

	int neglink = GaugeGrid->GlobalDimensions()[mu]-1;
	
	if ( gpdirs[mu] ) { 
	  Uconj = where(coor==neglink,-Uconj,Uconj);
	}

PARALLEL_FOR_LOOP
        for(auto ss=U.begin();ss<U.end();ss++){
	  Uds[ss](0)(mu) = U[ss]();
	  Uds[ss](1)(mu) = Uconj[ss]();
	}
 
        U     = adj(Cshift(U    ,mu,-1));      // correct except for spanning the boundary
        Uconj = adj(Cshift(Uconj,mu,-1));

	Utmp = U;
	if ( gpdirs[mu] ) { 
	  Utmp = where(coor==0,Uconj,Utmp);
	}

PARALLEL_FOR_LOOP
        for(auto ss=U.begin();ss<U.end();ss++){
	  Uds[ss](0)(mu+4) = Utmp[ss]();
	}

        Utmp = Uconj;
	if ( gpdirs[mu] ) { 
	  Utmp = where(coor==0,U,Utmp);
	}

PARALLEL_FOR_LOOP
	for(auto ss=U.begin();ss<U.end();ss++){
	  Uds[ss](1)(mu+4) = Utmp[ss]();
	}

      }
    }

    };

    // Could QCD, SU2, QED etc...
    // Fund, adj, etc...
    // Composition with smeared link, bc's etc.. probably need multiple inheritance
    typedef WilsonImpl<vComplex ,Nc> WilsonImplR; // Real.. whichever prec
    typedef WilsonImpl<vComplexF,Nc> WilsonImplF; // Float
    typedef WilsonImpl<vComplexD,Nc> WilsonImplD; // Double

    typedef GparityWilsonImpl<vComplex ,Nc> GparityWilsonImplR; // Real.. whichever prec
    typedef GparityWilsonImpl<vComplexF,Nc> GparityWilsonImplF; // Float
    typedef GparityWilsonImpl<vComplexD,Nc> GparityWilsonImplD; // Double


  }
}
#endif
