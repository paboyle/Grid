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
    ////////////////////////////
    //ContFrac:
    //  Ls always odd. Rational poly deg is either Ls or Ls-1
    //PartFrac 
    //  Ls always odd. Rational poly deg is either Ls or Ls-1
    //
    //Cayley: Ls always even, Rational poly deg is Ls
    // 
    // Just set nrational as Ls. Forget about Ls-1 cases.
    //
    // Require odd Ls for cont and part frac
    ////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    class WilsonFermion5D : public FermionOperator<LatticeFermion,LatticeGaugeField>
    {
    public:
      ///////////////////////////////////////////////////////////////
      // Implement the abstract base
      ///////////////////////////////////////////////////////////////
      GridBase *GaugeGrid(void)              { return _FourDimGrid ;}
      GridBase *GaugeRedBlackGrid(void)      { return _FourDimRedBlackGrid ;}
      GridBase *FermionGrid(void)            { return _FiveDimGrid;}
      GridBase *FermionRedBlackGrid(void)    { return _FiveDimRedBlackGrid;}

      // full checkerboard operations; leave unimplemented as abstract for now
      //virtual RealD  M    (const LatticeFermion &in, LatticeFermion &out)=0;
      //virtual RealD  Mdag (const LatticeFermion &in, LatticeFermion &out)=0;

      // half checkerboard operations; leave unimplemented as abstract for now
      //      virtual void   Meooe       (const LatticeFermion &in, LatticeFermion &out)=0;
      //      virtual void   MeooeDag    (const LatticeFermion &in, LatticeFermion &out)=0;
      //      virtual void   Mooee       (const LatticeFermion &in, LatticeFermion &out)=0;
      //      virtual void   MooeeDag    (const LatticeFermion &in, LatticeFermion &out)=0;
      //      virtual void   MooeeInv    (const LatticeFermion &in, LatticeFermion &out)=0;
      //      virtual void   MooeeInvDag (const LatticeFermion &in, LatticeFermion &out)=0;

      // Implement hopping term non-hermitian hopping term; half cb or both
      // Implement s-diagonal DW
      void DW    (const LatticeFermion &in, LatticeFermion &out,int dag);
      void Dhop  (const LatticeFermion &in, LatticeFermion &out,int dag);
      void DhopOE(const LatticeFermion &in, LatticeFermion &out,int dag);
      void DhopEO(const LatticeFermion &in, LatticeFermion &out,int dag);

      // add a DhopComm
      // -- suboptimal interface will presently trigger multiple comms.
      void DhopDir(const LatticeFermion &in, LatticeFermion &out,int dir,int disp);

      ///////////////////////////////////////////////////////////////
      // New methods added 
      ///////////////////////////////////////////////////////////////
      void DhopInternal(CartesianStencil & st,
			LebesgueOrder &lo,
			LatticeDoubledGaugeField &U,
			const LatticeFermion &in, 
			LatticeFermion &out,
			int dag);

      // Constructors
      WilsonFermion5D(LatticeGaugeField &_Umu,
			  GridCartesian         &FiveDimGrid,
			  GridRedBlackCartesian &FiveDimRedBlackGrid,
			  GridCartesian         &FourDimGrid,
			  GridRedBlackCartesian &FourDimRedBlackGrid,
			  double _M5);

      // DoubleStore
      void DoubleStore(LatticeDoubledGaugeField &Uds,const LatticeGaugeField &Umu);

      ///////////////////////////////////////////////////////////////
      // Data members require to support the functionality
      ///////////////////////////////////////////////////////////////
      static int HandOptDslash; // these are a temporary hack

    protected:

      // Add these to the support from Wilson
      GridBase *_FourDimGrid;
      GridBase *_FourDimRedBlackGrid;
      GridBase *_FiveDimGrid;
      GridBase *_FiveDimRedBlackGrid;

      static const int npoint=8;
      static const std::vector<int> directions   ;
      static const std::vector<int> displacements;

      double                        M5;
      int Ls;

      //Defines the stencils for even and odd
      CartesianStencil Stencil; 
      CartesianStencil StencilEven; 
      CartesianStencil StencilOdd; 

      // Copy of the gauge field , with even and odd subsets
      LatticeDoubledGaugeField Umu;
      LatticeDoubledGaugeField UmuEven;
      LatticeDoubledGaugeField UmuOdd;

      LebesgueOrder Lebesgue;
      LebesgueOrder LebesgueEvenOdd;

      // Comms buffer
      std::vector<vHalfSpinColourVector,alignedAllocator<vHalfSpinColourVector> >  comm_buf;
      
    };
  }
}

#endif
