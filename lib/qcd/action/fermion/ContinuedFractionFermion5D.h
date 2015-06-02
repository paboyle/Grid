#ifndef  GRID_QCD_CONTINUED_FRACTION_H
#define  GRID_QCD_CONTINUED_FRACTION_H

namespace Grid {

  namespace QCD {

    class ContinuedFractionFermion5D : public WilsonFermion5D
    {
    public:

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

    private:

      Approx::zolotarev_data *zdata;

      // Cont frac
      RealD mass;
      RealD R;
      RealD scale;
      std::vector<double> Beta;
      std::vector<double> cc;;
      std::vector<double> cc_d;;
      std::vector<double> sqrt_cc;
      std::vector<double> See;
      std::vector<double> Aee;

      // Constructors
      ContinuedFractionFermion5D(LatticeGaugeField &_Umu,
				 GridCartesian         &FiveDimGrid,
				 GridRedBlackCartesian &FiveDimRedBlackGrid,
				 GridCartesian         &FourDimGrid,
				 GridRedBlackCartesian &FourDimRedBlackGrid,
				 RealD _mass,RealD M5);

    };


  }
}

#endif
