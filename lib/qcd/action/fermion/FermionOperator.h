#ifndef  GRID_QCD_FERMION_OPERATOR_H
#define  GRID_QCD_FERMION_OPERATOR_H

namespace Grid {

  namespace QCD {

    ////////////////////////////////////////////////////////////////
    // Hardwire to four spinors, allow to select 
    // between gauge representation rank, and gparity/flavour index,
    // and single/double precision.
    ////////////////////////////////////////////////////////////////

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

      // provide the multiply by link that is differentiated between Gparity (with flavour index) and non-Gparity
      static inline void multLink(SiteHalfSpinor &phi,const SiteDoubledGaugeField &U,const SiteHalfSpinor &chi,int mu,StencilEntry *SE){
        mult(&phi(),&U(mu),&chi());
      }
      static inline void DoubleStore(GridBase *GaugeGrid,DoubledGaugeField &Uds,const GaugeField &Umu)
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
      static inline void InsertForce(GaugeField &mat,const FermionField &Btilde,const FermionField &A,int mu){
	GaugeLinkField link(mat._grid);
	link = TraceIndex<SpinIndex>(outerProduct(Btilde,A)); 
	PokeIndex<LorentzIndex>(mat,link,mu);
      }   

    };

    typedef WilsonImpl<vComplex,Nc>  WilsonImplR; // Real.. whichever prec
    typedef WilsonImpl<vComplexF,Nc> WilsonImplF; // Float
    typedef WilsonImpl<vComplexD,Nc> WilsonImplD; // Double

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
    static inline void multLink(SiteHalfSpinor &phi,const SiteDoubledGaugeField &U,const SiteHalfSpinor &chi,int mu,StencilEntry *SE){
      // FIXME; need to be more careful. If this is a simd direction we are still stuffed
      if ( SE->_around_the_world && ((mu==Xp)||(mu==Xm)) ) {
	mult(&phi(0),&U(0)(mu),&chi(1));
	mult(&phi(1),&U(1)(mu),&chi(0));
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
	
	std::vector<int> gpdirs({1,0,0,0});

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

    typedef GparityWilsonImpl<vComplex,Nc>  GparityWilsonImplR; // Real.. whichever prec
    typedef GparityWilsonImpl<vComplexF,Nc> GparityWilsonImplF; // Float
    typedef GparityWilsonImpl<vComplexD,Nc> GparityWilsonImplD; // Double


    
    //////////////////////////////////////////////////////////////////////////////
    // Four component fermions
    // Should type template the vector and gauge types
    // Think about multiple representations
    //////////////////////////////////////////////////////////////////////////////
    template<class Impl>
    class FermionOperator : public CheckerBoardedSparseMatrixBase<typename Impl::FermionField>
    {
    public:
#include <qcd/action/fermion/FermionImplTypedefs.h>
    public:

      GridBase * Grid(void)   { return FermionGrid(); };   // this is all the linalg routines need to know
      GridBase * RedBlackGrid(void) { return FermionRedBlackGrid(); };

      virtual GridBase *FermionGrid(void)         =0;
      virtual GridBase *FermionRedBlackGrid(void) =0;
      virtual GridBase *GaugeGrid(void)           =0;
      virtual GridBase *GaugeRedBlackGrid(void)   =0;

      // override multiply
      virtual RealD  M    (const FermionField &in, FermionField &out)=0;
      virtual RealD  Mdag (const FermionField &in, FermionField &out)=0;

      // half checkerboard operaions
      virtual void   Meooe       (const FermionField &in, FermionField &out)=0;
      virtual void   MeooeDag    (const FermionField &in, FermionField &out)=0;
      virtual void   Mooee       (const FermionField &in, FermionField &out)=0;
      virtual void   MooeeDag    (const FermionField &in, FermionField &out)=0;
      virtual void   MooeeInv    (const FermionField &in, FermionField &out)=0;
      virtual void   MooeeInvDag (const FermionField &in, FermionField &out)=0;

      // non-hermitian hopping term; half cb or both
      virtual void Dhop  (const FermionField &in, FermionField &out,int dag)=0;
      virtual void DhopOE(const FermionField &in, FermionField &out,int dag)=0;
      virtual void DhopEO(const FermionField &in, FermionField &out,int dag)=0;
      virtual void DhopDir(const FermionField &in, FermionField &out,int dir,int disp)=0; // implemented by WilsonFermion and WilsonFermion5D

      // force terms; five routines; default to Dhop on diagonal
      virtual void MDeriv  (GaugeField &mat,const FermionField &U,const FermionField &V,int dag){DhopDeriv(mat,U,V,dag);};
      virtual void MoeDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag){DhopDerivOE(mat,U,V,dag);};
      virtual void MeoDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag){DhopDerivEO(mat,U,V,dag);};
      virtual void MooDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag){mat=zero;};
      virtual void MeeDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag){mat=zero;};

      virtual void DhopDeriv  (GaugeField &mat,const FermionField &U,const FermionField &V,int dag)=0;
      virtual void DhopDerivEO(GaugeField &mat,const FermionField &U,const FermionField &V,int dag)=0;
      virtual void DhopDerivOE(GaugeField &mat,const FermionField &U,const FermionField &V,int dag)=0;


      virtual void  Mdiag  (const FermionField &in, FermionField &out) { Mooee(in,out);};   // Same as Mooee applied to both CB's
      virtual void  Mdir   (const FermionField &in, FermionField &out,int dir,int disp)=0;   // case by case Wilson, Clover, Cayley, ContFrac, PartFrac

      ///////////////////////////////////////////////
      // Updates gauge field during HMC
      ///////////////////////////////////////////////
      virtual void ImportGauge(const GaugeField & _U)=0;

    };

  }
}

#endif
