/*
  @file stoutSmear.hpp
  @brief Declares Stout smearing class
*/
#ifndef STOUT_SMEAR_
#define STOUT_SMEAR_

namespace Grid {
  namespace QCD {

    /*!  @brief Stout smearing of link variable. */
    template <class Gimpl> 
    class Smear_Stout: public Smear<Gimpl> {
    private:
      const std::vector<double> d_rho;
      const Smear < Gimpl > * SmearBase;
      
      LatticeReal func_xi0(LatticeReal w) const{
	// Define a function to do the check
	//if( w < 1e-4 ) std::cout << GridLogWarning << "[Smear_stout] w too small: "<< w <<"\n";
	return  sin(w)/w;
      }
      
    public:
      INHERIT_GIMPL_TYPES(Gimpl)
      
      Smear_Stout(Smear < Gimpl >* base):SmearBase(base){}
      
      /*! Default constructor */
      Smear_Stout():SmearBase(new Smear_APE < Gimpl > ()){}
      
      ~Smear_Stout(){}
      
      void smear(GaugeField& u_smr,const GaugeField& U) const{
	long double timing;
	
	GaugeField u_tmp1, q_mu;
	
	std::cout<< GridLogDebug << "Stout smearing started\n";
	
	//Smear the configurations
	SmearBase->smear(u_tmp1, U);
	
	q_mu = Ta(u_tmp1*adj(u_tmp1)); // q_mu = Ta(Omega_mu)
	
	exponentiate_iQ(u_tmp1, q_mu);

	u_smr = u_tmp1*U;
	
	std::cout<< GridLogDebug << "Stout smearing completed\n";
      }
      void derivative(GaugeField& SigmaTerm,
		      const GaugeField& iLambda,
		      const GaugeField& Gauge) const{
	SmearBase->derivative(SigmaTerm, iLambda, Gauge);
      }
      
      
      void BaseSmear(GaugeField& C,
		     const GaugeField& U) const{
	SmearBase->smear(C, U);
      }
      
      void exponentiate_iQ(GaugeField& e_iQ,
			   const GaugeField& iQ) const{
	// Put this outside 
	// only valid for SU(3) matrices

	GridBase *grid = iQ._grid;
	Real one_over_three = 1.0/3.0;
	Real one_over_two = 1.0/2.0;

	GaugeField unity;
	GaugeLinkField Umu(iQ._grid);
	Umu=1.0;
	for(int mu=0;mu<Nd;mu++){
	  pokeLorentz(unity,Umu,mu);
	}

	
	GaugeField iQ2, iQ3;
	LatticeReal c0(grid), c1(grid), c0max(grid), u_val(grid), tmp(grid);
	LatticeReal w(grid), theta(grid), xi0(grid), u2(grid), w2(grid), cosw(grid);
	LatticeComplex fden(grid);
	LatticeComplex f0(grid), f1(grid), f2(grid), h0(grid), h1(grid), h2(grid);
	LatticeComplex e2iu(grid), emiu(grid), ixi0(grid), qt(grid);

	iQ2 = iQ * iQ;
	iQ3 = iQ * iQ2;

	c0    = - imag(trace(iQ3)) * one_over_three;
	c1    = - real(trace(iQ2)) * one_over_two;
	tmp   = c1 * one_over_three;
    	c0max = 2.0 * pow(tmp, 1.5);

	theta = acos(c0/c0max);

	u_val = sqrt(tmp) * cos( theta * one_over_three);
	w     = sqrt(c1) * sin ( theta * one_over_three);
	xi0   = func_xi0(w);
	u2    = u_val * u_val;
	w2    = w * w;
	cosw  = cos(w);

	ixi0  = timesI(toComplex(xi0));
	emiu  = toComplex(cos(u_val)) - timesI(toComplex(u_val));
	e2iu  = toComplex(cos(2.0*u_val)) + timesI(toComplex(2.0*u_val));
	
	h0    = e2iu * toComplex(u2 - w2) + emiu *( toComplex(8.0*u2*cosw) +
						    toComplex(2.0*u_val*(3.0*u2 + w2))*ixi0);
	
	h1    = toComplex(2.0*u_val) * e2iu - emiu*( toComplex(2.0*u_val*cosw) -
						     toComplex(3.0*u2-w2)*ixi0);
	
	h2    = e2iu - emiu * (toComplex(cosw) + toComplex(3.0*u_val)*ixi0);

	
	tmp   = 9.0*u2 - w2;
	fden  = toComplex(pow(tmp, -1.0));
	f0    = h0 * fden;
	f1    = h1 * fden;
	f2    = h2 * fden;

	
	e_iQ = f0*unity + f1 * timesMinusI(iQ) - f2 * iQ2;
	
	
	
      };
      
    };

  }
}

#endif  
