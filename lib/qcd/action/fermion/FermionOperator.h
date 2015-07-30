#ifndef  GRID_QCD_FERMION_OPERATOR_H
#define  GRID_QCD_FERMION_OPERATOR_H

namespace Grid {

  namespace QCD {

    //////////////////////////////////////////////////////////////////////////////
    // Four component fermions
    // Should type template the vector and gauge types
    // Think about multiple representations
    //////////////////////////////////////////////////////////////////////////////
    template<class FermionField,class GaugeField>
    class FermionOperator : public CheckerBoardedSparseMatrixBase<FermionField>
    {
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
      virtual void MDeriv  (LatticeGaugeField &mat,const FermionField &U,const FermionField &V,int dag){DhopDeriv(mat,U,V,dag);};
      virtual void MoeDeriv(LatticeGaugeField &mat,const FermionField &U,const FermionField &V,int dag){DhopDerivOE(mat,U,V,dag);};
      virtual void MeoDeriv(LatticeGaugeField &mat,const FermionField &U,const FermionField &V,int dag){DhopDerivEO(mat,U,V,dag);};
      virtual void MooDeriv(LatticeGaugeField &mat,const FermionField &U,const FermionField &V,int dag){mat=zero;};
      virtual void MeeDeriv(LatticeGaugeField &mat,const FermionField &U,const FermionField &V,int dag){mat=zero;};

      virtual void DhopDeriv  (LatticeGaugeField &mat,const FermionField &U,const FermionField &V,int dag)=0;
      virtual void DhopDerivEO(LatticeGaugeField &mat,const FermionField &U,const FermionField &V,int dag)=0;
      virtual void DhopDerivOE(LatticeGaugeField &mat,const FermionField &U,const FermionField &V,int dag)=0;


      virtual void  Mdiag  (const FermionField &in, FermionField &out) { Mooee(in,out);};   // Same as Mooee applied to both CB's
      virtual void  Mdir   (const FermionField &in, FermionField &out,int dir,int disp)=0;   // case by case Wilson, Clover, Cayley, ContFrac, PartFrac

      ///////////////////////////////////////////////
      // Updates gauge field during HMC
      ///////////////////////////////////////////////
      virtual void ImportGauge(const GaugeField & _U);

    };

  }
}
#endif
