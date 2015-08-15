#ifndef  GRID_QCD_WILSON_FERMION_5D_H
#define  GRID_QCD_WILSON_FERMION_5D_H

namespace Grid {

  namespace QCD {

    ////////////////////////////////////////////////////////////////////////////////
    // This is the 4d red black case appropriate to support
    //
    // parity = (x+y+z+t)|2;
    // generalised five dim fermions like mobius, zolotarev etc..	
    //
    // i.e. even even contains fifth dim hopping term.
    //
    // [DIFFERS from original CPS red black implementation parity = (x+y+z+t+s)|2 ]
    ////////////////////////////////////////////////////////////////////////////////

    class WilsonFermion5DStatic { 
    public:
      // S-direction is INNERMOST and takes no part in the parity.
      static int HandOptDslash; // these are a temporary hack
      static const std::vector<int> directions;
      static const std::vector<int> displacements;
      const int npoint = 8;
    };

    template<class Impl>
    class WilsonFermion5D : public WilsonKernels<Impl>, public WilsonFermion5DStatic
    {
     INHERIT_IMPL_TYPES(Impl);
     typedef WilsonKernels<Impl> Kernels;

    public:
      ///////////////////////////////////////////////////////////////
      // Implement the abstract base
      ///////////////////////////////////////////////////////////////
      GridBase *GaugeGrid(void)              { return _FourDimGrid ;}
      GridBase *GaugeRedBlackGrid(void)      { return _FourDimRedBlackGrid ;}
      GridBase *FermionGrid(void)            { return _FiveDimGrid;}
      GridBase *FermionRedBlackGrid(void)    { return _FiveDimRedBlackGrid;}

      // full checkerboard operations; leave unimplemented as abstract for now
      virtual RealD  M    (const FermionField &in, FermionField &out){assert(0); return 0.0;};
      virtual RealD  Mdag (const FermionField &in, FermionField &out){assert(0); return 0.0;};

      // half checkerboard operations; leave unimplemented as abstract for now
      virtual void   Meooe       (const FermionField &in, FermionField &out){assert(0);};
      virtual void   Mooee       (const FermionField &in, FermionField &out){assert(0);};
      virtual void   MooeeInv    (const FermionField &in, FermionField &out){assert(0);};

      virtual void   MeooeDag    (const FermionField &in, FermionField &out){assert(0);};
      virtual void   MooeeDag    (const FermionField &in, FermionField &out){assert(0);};
      virtual void   MooeeInvDag (const FermionField &in, FermionField &out){assert(0);};

      // These can be overridden by fancy 5d chiral action
      virtual void DhopDeriv  (GaugeField &mat,const FermionField &U,const FermionField &V,int dag);
      virtual void DhopDerivEO(GaugeField &mat,const FermionField &U,const FermionField &V,int dag);
      virtual void DhopDerivOE(GaugeField &mat,const FermionField &U,const FermionField &V,int dag);

      // Implement hopping term non-hermitian hopping term; half cb or both
      // Implement s-diagonal DW
      void DW    (const FermionField &in, FermionField &out,int dag);
      void Dhop  (const FermionField &in, FermionField &out,int dag);
      void DhopOE(const FermionField &in, FermionField &out,int dag);
      void DhopEO(const FermionField &in, FermionField &out,int dag);

      // add a DhopComm
      // -- suboptimal interface will presently trigger multiple comms.
      void DhopDir(const FermionField &in, FermionField &out,int dir,int disp);

      ///////////////////////////////////////////////////////////////
      // New methods added 
      ///////////////////////////////////////////////////////////////
      void DerivInternal(CartesianStencil & st,
			 DoubledGaugeField & U,
			 GaugeField &mat,
			 const FermionField &A,
			 const FermionField &B,
			 int dag);

      void DhopInternal(CartesianStencil & st,
			LebesgueOrder &lo,
			DoubledGaugeField &U,
			const FermionField &in, 
			FermionField &out,
			int dag);

      // Constructors
      WilsonFermion5D(GaugeField &_Umu,
		      GridCartesian         &FiveDimGrid,
		      GridRedBlackCartesian &FiveDimRedBlackGrid,
		      GridCartesian         &FourDimGrid,
		      GridRedBlackCartesian &FourDimRedBlackGrid,
		      double _M5,const ImplParams &p= ImplParams());

      // DoubleStore
      void ImportGauge(const GaugeField &_Umu);

      ///////////////////////////////////////////////////////////////
      // Data members require to support the functionality
      ///////////////////////////////////////////////////////////////
    protected:

      // Add these to the support from Wilson
      GridBase *_FourDimGrid;
      GridBase *_FourDimRedBlackGrid;
      GridBase *_FiveDimGrid;
      GridBase *_FiveDimRedBlackGrid;

      double                        M5;
      int Ls;

      //Defines the stencils for even and odd
      CartesianStencil Stencil; 
      CartesianStencil StencilEven; 
      CartesianStencil StencilOdd; 

      // Copy of the gauge field , with even and odd subsets
      DoubledGaugeField Umu;
      DoubledGaugeField UmuEven;
      DoubledGaugeField UmuOdd;

      LebesgueOrder Lebesgue;
      LebesgueOrder LebesgueEvenOdd;

      // Comms buffer
      std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  comm_buf;
      
    };
  }
}

#endif
