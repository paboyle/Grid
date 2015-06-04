#ifndef  GRID_QCD_PARTIAL_FRACTION_H
#define  GRID_QCD_PARTIAL_FRACTION_H

namespace Grid {

  namespace QCD {

    class PartialFractionFermion5D : public WilsonFermion5D
    {
    public:

      const int part_frac_chroma_convention=1;

      void   Meooe_internal(const LatticeFermion &in, LatticeFermion &out,int dag);
      void   Mooee_internal(const LatticeFermion &in, LatticeFermion &out,int dag);
      void   MooeeInv_internal(const LatticeFermion &in, LatticeFermion &out,int dag);
      void   M_internal(const LatticeFermion &in, LatticeFermion &out,int dag);

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

      virtual void   Instantiatable(void) =0; // ensure no make-eee

      // Constructors
      PartialFractionFermion5D(LatticeGaugeField &_Umu,
				    GridCartesian         &FiveDimGrid,
				    GridRedBlackCartesian &FiveDimRedBlackGrid,
				    GridCartesian         &FourDimGrid,
				    GridRedBlackCartesian &FourDimRedBlackGrid,
				    RealD _mass,RealD M5);

    protected:

      virtual void SetCoefficientsTanh(Approx::zolotarev_data *zdata,RealD scale);
      virtual void SetCoefficientsZolotarev(RealD zolo_hi,Approx::zolotarev_data *zdata);

      // Part frac
      RealD mass;
      RealD dw_diag;
      RealD R;
      RealD amax;
      RealD scale;
      std::vector<double> p; 
      std::vector<double> q;

    };


  }
}

#endif
