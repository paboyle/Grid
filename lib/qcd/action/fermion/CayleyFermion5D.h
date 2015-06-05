#ifndef  GRID_QCD_CAYLEY_FERMION_H
#define  GRID_QCD_CAYLEY_FERMION_H

namespace Grid {

  namespace QCD {

    class CayleyFermion5D : public WilsonFermion5D
    {
    public:

      // override multiply
      virtual RealD  M    (const LatticeFermion &in, LatticeFermion &out);
      virtual RealD  Mdag (const LatticeFermion &in, LatticeFermion &out);

      // half checkerboard operations
      virtual void   Meooe       (const LatticeFermion &in, LatticeFermion &out);
      virtual void   MeooeDag    (const LatticeFermion &in, LatticeFermion &out);
      virtual void   Mooee       (const LatticeFermion &in, LatticeFermion &out);
      virtual void   MooeeDag    (const LatticeFermion &in, LatticeFermion &out);
      virtual void   MooeeInv    (const LatticeFermion &in, LatticeFermion &out);
      virtual void   MooeeInvDag (const LatticeFermion &in, LatticeFermion &out);
      virtual void   Instantiatable(void)=0;
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
      CayleyFermion5D(LatticeGaugeField &_Umu,
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
