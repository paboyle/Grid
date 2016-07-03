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

  		public:
  			INHERIT_GIMPL_TYPES(Gimpl)

  			Smear_Stout(Smear < Gimpl >* base):SmearBase(base){
  				static_assert(Nc==3, "Stout smearing currently implemented only for Nc==3");
  			}

      /*! Default constructor */
  			Smear_Stout(double rho = 1.0):SmearBase(new Smear_APE < Gimpl > (rho)){
  				static_assert(Nc==3, "Stout smearing currently implemented only for Nc==3");
  			}

  			~Smear_Stout(){}

  			void smear(GaugeField& u_smr,const GaugeField& U) const{
  				GaugeField C(U._grid);
  				GaugeLinkField tmp(U._grid), iq_mu(U._grid), Umu(U._grid);

  				std::cout<< GridLogDebug << "Stout smearing started\n";

	//Smear the configurations
  				SmearBase->smear(C, U);
  				for (int mu = 0; mu<Nd; mu++)
  				{
  					tmp = peekLorentz(C,mu);
  					Umu = peekLorentz(U,mu);
		  			iq_mu = Ta(tmp * adj(Umu)); // iq_mu = Ta(Omega_mu) to match the signs with the paper
		  			exponentiate_iQ(tmp, iq_mu);  
		  			GaugeLinkField check = adj(tmp) * tmp - 1.0;
					pokeLorentz(u_smr, tmp*Umu, mu);// u_smr = exp(iQ_mu)*U_mu
		  		}

		  		std::cout<< GridLogDebug << "Stout smearing completed\n";


		  	};


		  	void derivative(GaugeField& SigmaTerm,
		  		const GaugeField& iLambda,
		  		const GaugeField& Gauge) const{
		  		SmearBase->derivative(SigmaTerm, iLambda, Gauge);
		  	};


		  	void BaseSmear(GaugeField& C,
		  		const GaugeField& U) const{
		  		SmearBase->smear(C, U);
		  	};

		  	void exponentiate_iQ(GaugeLinkField& e_iQ,
		  		const GaugeLinkField& iQ) const{
		// Put this outside 
		// only valid for SU(3) matrices

		// only one Lorentz direction at a time 

		// notice that it actually computes
		// exp ( input matrix )  
		// the i sign is coming from outside
		// input matrix is anti-hermitian NOT hermitian

		  		GridBase *grid = iQ._grid;
		  		GaugeLinkField unity(grid);
		  		unity=1.0;

		  		GaugeLinkField iQ2(grid), iQ3(grid);
		  		LatticeComplex u(grid), w(grid);
		  		LatticeComplex f0(grid), f1(grid), f2(grid);

		  		iQ2 = iQ * iQ;
		  		iQ3 = iQ * iQ2;

		  		set_uw_complex(u, w, iQ2, iQ3);
		  		set_fj_complex(f0, f1, f2, u, w);
		  		

		  		e_iQ = f0*unity + timesMinusI(f1) * iQ - f2 * iQ2;
		  	};


		  	void set_uw(LatticeReal& u, LatticeReal& w,
		  		GaugeLinkField& iQ2, GaugeLinkField& iQ3) const{
		  		Real one_over_three = 1.0/3.0;
		  		Real one_over_two = 1.0/2.0;

		  		GridBase *grid = u._grid;
		  		LatticeReal c0(grid), c1(grid), tmp(grid), c0max(grid), theta(grid);

				// sign in c0 from the conventions on the Ta
				//	c0    = - toReal(imag(trace(iQ3))) * one_over_three;
				c0    = - toReal(real(timesMinusI(trace(iQ3)))) * one_over_three; //slow and temporary, FIX the bug in imag
				c1    = - toReal(real(trace(iQ2))) * one_over_two;
				tmp   = c1 * one_over_three;
				c0max = 2.0 * pow(tmp, 1.5);

				theta = acos(c0/c0max);
				u = sqrt(tmp) * cos( theta * one_over_three);
				w = sqrt(c1)  * sin( theta * one_over_three);

			}

			void set_uw_complex(LatticeComplex& u, LatticeComplex& w,
		  		GaugeLinkField& iQ2, GaugeLinkField& iQ3) const{
		  		Complex one_over_three = 1.0/3.0;
		  		Complex one_over_two = 1.0/2.0;

		  		GridBase *grid = u._grid;
		  		LatticeComplex c0(grid), c1(grid), tmp(grid), c0max(grid), theta(grid);

				// sign in c0 from the conventions on the Ta
				c0    = - real(timesMinusI(trace(iQ3))) * one_over_three; //temporary hack 
				c1    = - real(trace(iQ2)) * one_over_two;
				tmp   = c1 * one_over_three;
				c0max = 2.0 * pow(tmp, 1.5);

				theta = acos(c0/c0max);
				u = sqrt(tmp) * cos( theta * one_over_three);
				w = sqrt(c1)  * sin( theta * one_over_three);

			}


			void set_fj(LatticeComplex& f0, LatticeComplex& f1, LatticeComplex& f2,
				const LatticeReal& u, const LatticeReal& w) const{

				GridBase *grid = u._grid;
				LatticeReal xi0(grid), u2(grid), w2(grid), cosw(grid), tmp(grid);
				LatticeComplex fden(grid);
				LatticeComplex h0(grid), h1(grid), h2(grid);
				LatticeComplex e2iu(grid), emiu(grid), ixi0(grid), qt(grid);

				xi0   = func_xi0(w);
				u2    = u * u;
				w2    = w * w;
				cosw  = cos(w);

				ixi0  = timesI(toComplex(xi0));
				emiu  = toComplex(cos(u)) - timesI(toComplex(sin(u)));
				e2iu  = toComplex(cos(2.0*u)) + timesI(toComplex(sin(2.0*u)));

				h0    = e2iu * toComplex(u2 - w2) + emiu *( toComplex(8.0*u2*cosw) +
					toComplex(2.0*u*(3.0*u2 + w2))*ixi0);

				h1    = toComplex(2.0*u) * e2iu - emiu*( toComplex(2.0*u*cosw) -
					toComplex(3.0*u2-w2)*ixi0);

				h2    = e2iu - emiu * (toComplex(cosw) + toComplex(3.0*u)*ixi0);

				tmp   = 9.0*u2 - w2;
				fden  = toComplex(pow(tmp, -1.0));
				f0    = h0 * fden;
				f1    = h1 * fden;
				f2    = h2 * fden;	


			}

			void set_fj_complex(LatticeComplex& f0, LatticeComplex& f1, LatticeComplex& f2,
				const LatticeComplex& u, const LatticeComplex& w) const{

				GridBase *grid = u._grid;
				LatticeComplex xi0(grid), u2(grid), w2(grid), cosw(grid), tmp(grid);
				LatticeComplex fden(grid);
				LatticeComplex h0(grid), h1(grid), h2(grid);
				LatticeComplex e2iu(grid), emiu(grid), ixi0(grid), qt(grid);

				xi0   = sin(w)/w;//func_xi0(w);
				u2    = u * u;
				w2    = w * w;
				cosw  = cos(w);

				ixi0  = timesI(xi0);
				emiu  = cos(u) - timesI(sin(u));
				e2iu  = cos(2.0*u) + timesI(sin(2.0*u));

				h0    = e2iu * (u2 - w2) + emiu *( (8.0*u2*cosw) +
					(2.0*u*(3.0*u2 + w2)*ixi0));

				h1    = (2.0*u) * e2iu - emiu*( (2.0*u*cosw) -
					(3.0*u2-w2)*ixi0);

				h2    = e2iu - emiu * (cosw + (3.0*u)*ixi0);

				tmp   = 9.0*u2 - w2;
				fden  = pow(tmp, -1.0);
				f0    = h0 * fden;
				f1    = h1 * fden;
				f2    = h2 * fden;	


			}




			LatticeReal func_xi0(const LatticeReal& w) const{
	// Define a function to do the check
	//if( w < 1e-4 ) std::cout << GridLogWarning<< "[Smear_stout] w too small: "<< w <<"\n";
				return  sin(w)/w;
			}

			LatticeReal func_xi1(const LatticeReal& w) const{
	// Define a function to do the check
	//if( w < 1e-4 ) std::cout << GridLogWarning << "[Smear_stout] w too small: "<< w <<"\n";
				return  cos(w)/(w*w) - sin(w)/(w*w*w);
			}

		};

	}
}

#endif  
