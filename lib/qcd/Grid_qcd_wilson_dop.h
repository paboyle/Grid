#ifnfdef GRID_QCD_WILSON_DOP_H
#define  GRID_QCD_WILSON_DOP_H

#include <Grid.h>

namespace Grid {

  namespace QCD {


    template<class vtype> class LinearOperatorBase {
    public:
      void multiply(const Lattice<vtype> &in, Lattice<vtype> &out){ assert(0);}
    };

    class WilsonMatrix : public LinearOperatorBase<LatticeFermion>
    {
      //NB r=1;
    public:
      double                        mass;
      GridBase                     *grid;

      // Copy of the gauge field 
      LatticeGaugeField             Umu;

      //Defines the stencil
      CartesianStencil              Stencil; 
      static const int npoint=9;
      static const std::vector<int> directions   ;
      static const std::vector<int> displacements;

      static const int Xp,Xm,Yp,Ym,Zp,Zm,Tp,Tm;

      // Comms buffer
      std::vector<vSpinColourVector,alignedAllocator<vSpinColourVector> >  comm_buf;

      // Constructor
      WilsonMatrix(LatticeGaugeField &Umu,int mass);

      // override multiply
      void multiply(const LatticeFermion &in, LatticeFermion &out);

      // non-hermitian hopping term; half cb or both
      void Dhop(const LatticeFermion &in, LatticeFermion &out);

      // m+4r -1/2 Dhop; both cb's
      void Dw(const LatticeFermion &in, LatticeFermion &out);

      // half checkerboard operaions
      void MpcDag   (const LatticeFermion &in, LatticeFermion &out);
      void Mpc      (const LatticeFermion &in, LatticeFermion &out);
      void MpcDagMpc(const LatticeFermion &in, LatticeFermion &out);

      // full checkerboard hermitian
      void MDagM    (const LatticeFermion &in, LatticeFermion &out);

      
    };

  }
}
#endif
