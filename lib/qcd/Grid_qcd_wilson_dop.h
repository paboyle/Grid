#ifndef  GRID_QCD_WILSON_DOP_H
#define  GRID_QCD_WILSON_DOP_H

#include <Grid.h>

#include <algorithms/LinearOperator.h>

namespace Grid {

  namespace QCD {

    class WilsonMatrix : public SparseMatrixBase<LatticeFermion>
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
      virtual void M    (const LatticeFermion &in, LatticeFermion &out);
      virtual void Mdag (const LatticeFermion &in, LatticeFermion &out);
      virtual void MdagM(const LatticeFermion &in, LatticeFermion &out);

      // half checkerboard operaions
      void Mpc      (const LatticeFermion &in, LatticeFermion &out);
      void MpcDag   (const LatticeFermion &in, LatticeFermion &out);
      void MpcDagMpc(const LatticeFermion &in, LatticeFermion &out);

      // non-hermitian hopping term; half cb or both
      void Dhop(const LatticeFermion &in, LatticeFermion &out);

      // m+4r -1/2 Dhop; both cb's
      void Dw(const LatticeFermion &in, LatticeFermion &out);

      typedef iScalar<iMatrix<vComplex, Nc> > matrix;

      
    };


  }
}
#endif
