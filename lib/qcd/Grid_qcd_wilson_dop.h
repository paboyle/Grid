#ifndef  GRID_QCD_WILSON_DOP_H
#define  GRID_QCD_WILSON_DOP_H


namespace Grid {

  namespace QCD {

    class WilsonMatrix : public CheckerBoardedSparseMatrixBase<LatticeFermion>
    {
      //NB r=1;
    public:
      double                        mass;
      GridBase                     *grid;

      // Copy of the gauge field 
      LatticeDoubledGaugeField             Umu;

      //Defines the stencil
      CartesianStencil              Stencil; 
      static const int npoint=9;
      static const std::vector<int> directions   ;
      static const std::vector<int> displacements;
      static const int Xp,Xm,Yp,Ym,Zp,Zm,Tp,Tm;

      // Comms buffer
      std::vector<vHalfSpinColourVector,alignedAllocator<vHalfSpinColourVector> >  comm_buf;

      // Constructor
      WilsonMatrix(LatticeGaugeField &Umu,double mass);

      // DoubleStore
      void DoubleStore(LatticeDoubledGaugeField &Uds,const LatticeGaugeField &Umu);

      // override multiply
      virtual RealD  M    (const LatticeFermion &in, LatticeFermion &out);
      virtual RealD  Mdag (const LatticeFermion &in, LatticeFermion &out);

      // half checkerboard operaions
      virtual void   Meooe       (const LatticeFermion &in, LatticeFermion &out);
      virtual void   MeooeDag    (const LatticeFermion &in, LatticeFermion &out);
      virtual void   Mooee       (const LatticeFermion &in, LatticeFermion &out);
      virtual void   MooeeDag    (const LatticeFermion &in, LatticeFermion &out);
      virtual void   MooeeInv    (const LatticeFermion &in, LatticeFermion &out);
      virtual void   MooeeInvDag (const LatticeFermion &in, LatticeFermion &out);

      // non-hermitian hopping term; half cb or both
      void Dhop(const LatticeFermion &in, LatticeFermion &out,int dag);
      void DhopSite   (int ss,const LatticeFermion &in, LatticeFermion &out);
      void DhopSiteDag(int ss,const LatticeFermion &in, LatticeFermion &out);

      typedef iScalar<iMatrix<vComplex, Nc> > matrix;

      
    };


  }
}
#endif
