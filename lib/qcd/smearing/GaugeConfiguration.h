/*!
  @file GaugeConfiguration.h

  @brief Declares the GaugeConfiguration class
*/
#ifndef GAUGE_CONFIG_
#define GAUGE_CONFIG_

namespace Grid {
  
  namespace QCD {
    
    /*!
      @brief Smeared configuration container
      
      It will behave like a configuration from the point of view of
      the HMC update and integrators.
      An "advanced configuration" object that can provide not only the 
      data to store the gauge configuration but also operations to manipulate
      it like smearing.
      
      It stores a list of smeared configurations.
    */
    template <class Gimpl>
    class SmearedConfiguration {
    public:
      INHERIT_GIMPL_TYPES(Gimpl)
      private:
      const unsigned int smearingLevels;
      Smear_Stout StoutSmearing;
      std::vector<GaugeField> SmearedSet;
      
      // Member functions
      void fill_smearedSet();
      GaugeField AnalyticSmearedForce(const GaugeField&, 
				      const GaugeField&) const;
      const GaugeField& get_smeared_conf(int) const;
      
      void set_iLambda(GaugeField& iLambda, 
		       GaugeField& e_iQ,
		       const GaugeField& iQ, 
		       const GaugeField& Sigmap,
		       const GaugeField& U)const;
      
      /* Check these types (do I need to pass iQ1,2 ? )
      void set_uw(RealD& u, RealD& w,
		  const SUNmat& iQ1, const SUNmat& iQ2)const ;
      void set_fj(ComplexD& f0, ComplexD& f1,
		  CompledD& f2, const RealD& u,
		  const RealD& w)const;
      */
      
      RealD func_xi0(RealD w)const;
      RealD func_xi1(RealD w)const; 
      
    public:
      GaugeField* ThinLinks;      /*!< @brief Pointer to the thin 
				    links configuration */
      
      /*! @brief Standard constructor */
      SmearedConfiguration(GridCartesian * UGrid,
			 unsigned int Nsmear, 
			 Smear_Stout& Stout):
	smearingLevels(Nsmear),
	StoutSmearing(Stout),
	ThinLinks(new GaugeField){
	for (unsigned int i=0; i< smearingLevels; ++i)
	  SmearedSet.push_back(*(new GaugeField(UGrid)));
      }
      
      /*! For just thin links */
      SmearedConfiguration(GridCartesian * UGrid):
	smearingLevels(0),
	StoutSmearing(),
	SmearedSet(0),
	ThinLinks(new GaugeField(UGrid)){}
      
      void set_GaugeField(){ fill_smearedSet(); }
      void smeared_force(GaugeField&) const;
      GaugeField* get_SmearedU() const{ 
	return const_cast<GaugeField*>(&(SmearedSet[smearingLevels-1]));
      }

      GaugeField* get_U(bool smeared=false) const { 
	// get the config, thin links by default
	if (smeared){
	  if (smearingLevels) return get_SmearedU();
	  else                return ThinLinks;
	}
	else return ThinLinks;
      }
      
    };
    
    
  }

}





#endif
