/*!
  @brief Declaration of Smear_APE class for APE smearing
*/

#ifndef APE_SMEAR_
#define APE_SMEAR_

namespace Grid {
  namespace QCD {


    /*!  @brief APE type smearing of link variables. */
    template <class Gimpl> 
    class Smear_APE: public Smear<Gimpl>{
    private:
      const std::vector<double> rho;/*!< Array of weights */

      //This member must be private - we do not want to control from outside 
      std::vector<double> set_rho(const double common_rho)const {
	std::vector<double> res;
    
	for(int mn=0; mn<Nd*Nd; ++mn) res.push_back(common_rho);
	for(int mu=0; mu<Nd; ++mu) res[mu + mu*Nd] = 0.0;
	return res;
      }

    public:
      // Defines the gauge field types
      INHERIT_GIMPL_TYPES(Gimpl)

  
      // Constructors and destructors
      Smear_APE(const std::vector<double>& rho_):rho(rho_){}
      Smear_APE(double rho_val):rho(set_rho(rho_val)){}
      Smear_APE():rho(set_rho(1.0)){}
      ~Smear_APE(){}


      void smear(GaugeField& u_smr, const GaugeField& U)const{
	GridBase *grid = U._grid;
	double d_rho;
	GaugeLinkField Cup(grid), tmp_stpl(grid);
	WilsonLoops<Gimpl> WL;
	u_smr = zero; // probably unecessary

	for(int mu=0; mu<Nd; ++mu){
	  Cup = zero;
	  for(int nu=0; nu<Nd; ++nu){
	    d_rho = rho[mu + Nd * nu];
	    WL.Staple(tmp_stpl, U, mu, nu);  //nb staple conventions of IroIro and Grid differ by a dag
	    Cup += tmp_stpl*d_rho;
	  }
	  pokeLorentz(u_smr, adj(Cup), mu); // u_smr[mu] = Cup^dag
	}
      }

      void derivative(GaugeField& SigmaTerm,
		      const GaugeField& iLambda,
		      const GaugeField& U)const{
    
	// Reference 
	// Morningstar, Peardon, Phys.Rev.D69,054501(2004)
	// Equation 75

    	GridBase *grid = U._grid;
	int vol = U._grid->gSites();
	
	WilsonLoops<Gimpl> WL;
	GaugeLinkField staple(grid), u_tmp(grid), iLambda_mu(grid), iLambda_nu(grid);
	GaugeLinkField U_mu(grid), U_nu(grid), sh_field(grid), temp_Sigma(grid);
	Real rho_munu, rho_numu;
    
	for(int mu = 0; mu < Nd; ++mu){
	  U_mu       = PeekIndex<LorentzIndex>(      U, mu);
	  iLambda_mu = PeekIndex<LorentzIndex>(iLambda, mu);
      
	  for(int nu = 0; nu < Nd; ++nu){
	    if(nu==mu) continue;
	    U_nu       = PeekIndex<LorentzIndex>(      U, nu);
	    iLambda_nu = PeekIndex<LorentzIndex>(iLambda, nu);
	
	    rho_munu = rho[mu + Nd * nu];
	    rho_numu = rho[nu + Nd * mu];
	
	    WL.StapleUpper(staple, U, mu, nu);

	    temp_Sigma = -rho_numu*staple*iLambda_nu;
	    //-r_numu*U_nu(x+mu)*Udag_mu(x+nu)*Udag_nu(x)*Lambda_nu(x)
	    AddGaugeLink(SigmaTerm, temp_Sigma, mu);

	    sh_field = Cshift(iLambda_nu, mu, 1);// general also for Gparity?
	
	    temp_Sigma = rho_numu*sh_field*staple;
	    //r_numu*Lambda_nu(mu)*U_nu(x+mu)*Udag_mu(x+nu)*Udag_nu(x)
	    AddGaugeLink(SigmaTerm, temp_Sigma, mu);

	    sh_field = Cshift(iLambda_mu, nu, 1);

	    temp_Sigma = -rho_munu*staple*U_nu*sh_field*adj(U_nu);
	    //-r_munu*U_nu(x+mu)*Udag_mu(x+nu)*Lambda_mu(x+nu)*Udag_nu(x)
	    AddGaugeLink(SigmaTerm, temp_Sigma, mu);

	    staple = zero;
	    sh_field = Cshift(U_nu, mu, 1);

	    temp_Sigma = -rho_munu*adj(sh_field)*adj(U_mu)*iLambda_mu*U_nu;
	    temp_Sigma += rho_numu*adj(sh_field)*adj(U_mu)*iLambda_nu*U_nu;

	    u_tmp = adj(U_nu)*iLambda_nu;
	    sh_field = Cshift(u_tmp, mu, 1);
	    temp_Sigma += -rho_numu*sh_field*adj(U_mu)*U_nu;
	    sh_field = Cshift(temp_Sigma, nu, -1);
	    AddGaugeLink(SigmaTerm, sh_field, mu);
	    
	  }
	}
      }
    };



  }// namespace QCD
}//namespace Grid
#endif  
