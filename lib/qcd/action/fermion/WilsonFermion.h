#ifndef  GRID_QCD_WILSON_FERMION_H
#define  GRID_QCD_WILSON_FERMION_H

namespace Grid {

  namespace QCD {

    class WilsonFermion : public FermionOperator<LatticeFermion,LatticeGaugeField>
    {
    public:

      ///////////////////////////////////////////////////////////////
      // Implement the abstract base
      ///////////////////////////////////////////////////////////////
      GridBase *GaugeGrid(void)              { return _grid ;}
      GridBase *GaugeRedBlackGrid(void)      { return _cbgrid ;}
      GridBase *FermionGrid(void)            { return _grid;}
      GridBase *FermionRedBlackGrid(void)    { return _cbgrid;}

      // override multiply
      virtual RealD  M    (const LatticeFermion &in, LatticeFermion &out);
      virtual RealD  Mdag (const LatticeFermion &in, LatticeFermion &out);

      // half checkerboard operaions
      void   Meooe       (const LatticeFermion &in, LatticeFermion &out);
      void   MeooeDag    (const LatticeFermion &in, LatticeFermion &out);
      virtual void   Mooee       (const LatticeFermion &in, LatticeFermion &out); // remain virtual so we 
      virtual void   MooeeDag    (const LatticeFermion &in, LatticeFermion &out); // can derive Clover
      virtual void   MooeeInv    (const LatticeFermion &in, LatticeFermion &out); // from Wilson bas
      virtual void   MooeeInvDag (const LatticeFermion &in, LatticeFermion &out);

      // non-hermitian hopping term; half cb or both
      void Dhop  (const LatticeFermion &in, LatticeFermion &out,int dag);
      void DhopOE(const LatticeFermion &in, LatticeFermion &out,int dag);
      void DhopEO(const LatticeFermion &in, LatticeFermion &out,int dag);

      ///////////////////////////////////////////////////////////////
      // Extra methods added by derived
      ///////////////////////////////////////////////////////////////
      void DhopInternal(CartesianStencil & st,
			LatticeDoubledGaugeField &U,
			const LatticeFermion &in, 
			LatticeFermion &out,
			int dag);

      // Constructor
      WilsonFermion(LatticeGaugeField &_Umu,GridCartesian &Fgrid,GridRedBlackCartesian &Hgrid,RealD _mass);

      // DoubleStore
      void DoubleStore(LatticeDoubledGaugeField &Uds,const LatticeGaugeField &Umu);

      ///////////////////////////////////////////////////////////////
      // Data members require to support the functionality
      ///////////////////////////////////////////////////////////////
      static int HandOptDslash; // these are a temporary hack
      static int MortonOrder;

    protected:

      RealD                        mass;

      GridBase                     *    _grid; 
      GridBase                     *  _cbgrid;

      static const int npoint=8;
      static const std::vector<int> directions   ;
      static const std::vector<int> displacements;

      //Defines the stencils for even and odd
      CartesianStencil Stencil; 
      CartesianStencil StencilEven; 
      CartesianStencil StencilOdd; 

      // Copy of the gauge field , with even and odd subsets
      LatticeDoubledGaugeField Umu;
      LatticeDoubledGaugeField UmuEven;
      LatticeDoubledGaugeField UmuOdd;

      // Comms buffer
      std::vector<vHalfSpinColourVector,alignedAllocator<vHalfSpinColourVector> >  comm_buf;

      
    };

  }
}
#endif
