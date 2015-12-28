#ifndef QCD_PLAQ_PLUS_RECTANGLE_ACTION_H
#define QCD_PLAQ_PLUS_RECTANGLE_ACTION_H

namespace Grid{
  namespace QCD{
    
    ////////////////////////////////////////////////////////////////////////
    // PlaqPlusRectangleActoin
    ////////////////////////////////////////////////////////////////////////
    template<class GaugeField>
    class PlaqPlusRectangleAction : public Action<GaugeField> {
    public:

      typedef LorentzScalar<GaugeField> GaugeLinkField;

    private:
      RealD c_plaq;
      RealD c_rect;

    public:
    PlaqPlusRectangleAction(RealD b,RealD c): c_plaq(b),c_rect(c){};
      
      virtual void refresh(const GaugeField &U, GridParallelRNG& pRNG) {}; // noop as no pseudoferms
      
      virtual RealD S(const GaugeField &U) {
	RealD vol = U._grid->gSites();

	RealD plaq = WilsonLoops<GaugeField>::avgPlaquette(U);
	RealD rect = WilsonLoops<GaugeField>::avgRectangle(U);

	RealD action=c_plaq*(1.0 -plaq)*(Nd*(Nd-1.0))*vol*0.5
	            +c_rect*(1.0 -rect)*(Nd*(Nd-1.0))*vol;

	return action;
      };

      virtual void deriv(const GaugeField &U,GaugeField & dSdU) {
	//extend Ta to include Lorentz indexes
	RealD factor_p = c_plaq/RealD(Nc)*0.5;
	RealD factor_r =   c_rect/RealD(Nc)*0.5;

	GaugeLinkField Umu(U._grid);
	GaugeLinkField dSdU_mu(U._grid);
	GaugeLinkField staple(U._grid);
	for (int mu=0; mu < Nd; mu++){

	  Umu = PeekIndex<LorentzIndex>(U,mu);

	  // Staple in direction mu
	  WilsonLoops<GaugeField>::Staple(staple,U,mu);
	  dSdU_mu = Ta(Umu*staple)*factor_p;
	  
	  WilsonLoops<GaugeField>::RectStaple(staple,U,mu);
	  dSdU_mu = dSdU_mu + Ta(Umu*staple)*factor_r;
	  
	  PokeIndex<LorentzIndex>(dSdU, dSdU_mu, mu);

	}
      };
    };

    // Convenience for common physically defined cases.
    //
    // RBC c1 parameterisation is not really RBC but don't have good
    // reference and we are happy to change name if prior use of this plaq coeff
    // parameterisation is made known to us. 
    template<class GaugeField>
    class RBCGaugeAction : public PlaqPlusRectangleAction<GaugeField> {
    public:
      RBCGaugeAction(RealD beta,RealD c1) : PlaqPlusRectangleAction<GaugeField>(beta*(1.0-8.0*c1), beta*c1) {
      };
    };

    template<class GaugeField>
    class IwasakiGaugeAction : public RBCGaugeAction<GaugeField> {
    public:
      IwasakiGaugeAction(RealD beta) : RBCGaugeAction<GaugeField>(beta,-0.331) {
      };
    };

    template<class GaugeField>
    class SymanzikGaugeAction : public RBCGaugeAction<GaugeField> {
    public:
      SymanzikGaugeAction(RealD beta) : RBCGaugeAction<GaugeField>(beta,-1.0/12.0) {
      };
    };

    template<class GaugeField>
    class DBW2GaugeAction : public RBCGaugeAction<GaugeField> {
    public:
      DBW2GaugeAction(RealD beta) : RBCGaugeAction<GaugeField>(beta,-1.4067) {
      };
    };

  }
}

#endif
