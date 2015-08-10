#ifndef  GRID_QCD_CAYLEY_FERMION_H
#define  GRID_QCD_CAYLEY_FERMION_H

namespace Grid {

  namespace QCD {

    template<class Impl>
    class CayleyFermion5D : public WilsonFermion5D<Impl>
    {
    public:
#include <qcd/action/fermion/FermionImplTypedefs.h>
    public:

      // override multiply
      virtual RealD  M    (const FermionField &in, FermionField &out);
      virtual RealD  Mdag (const FermionField &in, FermionField &out);

      // half checkerboard operations
      virtual void   Meooe       (const FermionField &in, FermionField &out);
      virtual void   MeooeDag    (const FermionField &in, FermionField &out);
      virtual void   Mooee       (const FermionField &in, FermionField &out);
      virtual void   MooeeDag    (const FermionField &in, FermionField &out);
      virtual void   MooeeInv    (const FermionField &in, FermionField &out);
      virtual void   MooeeInvDag (const FermionField &in, FermionField &out);
      virtual void   Instantiatable(void)=0;

      // force terms; five routines; default to Dhop on diagonal
      virtual void MDeriv  (GaugeField &mat,const FermionField &U,const FermionField &V,int dag);
      virtual void MoeDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag);
      virtual void MeoDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag);

      // Efficient support for multigrid coarsening
      virtual void  Mdir (const FermionField &in, FermionField &out,int dir,int disp);
      
      void   Meooe5D       (const FermionField &in, FermionField &out);
      void   MeooeDag5D    (const FermionField &in, FermionField &out);

      //    protected:
      RealD mass;

      // Cayley form Moebius (tanh and zolotarev)
      std::vector<RealD> omega; 
      std::vector<RealD> bs;    // S dependent coeffs
      std::vector<RealD> cs;    
      std::vector<RealD> as;    
      // For preconditioning Cayley form
      std::vector<RealD> bee;    
      std::vector<RealD> cee;    
      std::vector<RealD> aee;    
      std::vector<RealD> beo;    
      std::vector<RealD> ceo;    
      std::vector<RealD> aeo;    
      // LDU factorisation of the eeoo matrix
      std::vector<RealD> lee;    
      std::vector<RealD> leem;    
      std::vector<RealD> uee;    
      std::vector<RealD> ueem;    
      std::vector<RealD> dee;    

      // Constructors
      CayleyFermion5D(GaugeField &_Umu,
		      GridCartesian         &FiveDimGrid,
		      GridRedBlackCartesian &FiveDimRedBlackGrid,
		      GridCartesian         &FourDimGrid,
		      GridRedBlackCartesian &FourDimRedBlackGrid,
		      RealD _mass,RealD _M5);

    protected:
      void SetCoefficientsZolotarev(RealD zolohi,Approx::zolotarev_data *zdata,RealD b,RealD c);
      void SetCoefficientsTanh(Approx::zolotarev_data *zdata,RealD b,RealD c);
    };

  }
}

#endif
