#ifndef  GRID_QCD_PARTIAL_FRACTION_H
#define  GRID_QCD_PARTIAL_FRACTION_H

namespace Grid {

  namespace QCD {

    class PartialFractionFermion5D : public WilsonFermion5D
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

      virtual void PartialFractionCoefficients(void);

      zolotarev_data *zdata;

      // Part frac
      double R=(1+this->mass)/(1-this->mass);
      std::vector<double> p; 
      std::vector<double> q;

      // Constructors
      PartialFractionFermion5D(LatticeGaugeField &_Umu,
				    GridCartesian         &FiveDimGrid,
				    GridRedBlackCartesian &FiveDimRedBlackGrid,
				    GridCartesian         &FourDimGrid,
				    GridRedBlackCartesian &FourDimRedBlackGrid,
				    RealD _mass,RealD M5);

    };


  }
}

#endif
