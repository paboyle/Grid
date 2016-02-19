/*!
  @brief Declaration of Smear_APE class for APE smearing
*/

#ifndef APE_SMEAR_
#define APE_SMEAR_

/*!  @brief APE type smearing of link variables. */

template <class Gimpl> 
class Smear_APE: public Smear<Gimpl>{
private:
  const std::vector<double> rho;/*!< Array of weights */

  //This member must be private - we do not want to control from outside 
  std::vector<double> set_rho(const double)const {
    std::vector<double> res;
    
    for(int mn=0; mn<Nd*Nd; ++mn) res.push_back(common_rho);
    for(int mu=0; mu<Nd; ++mu) res[mu + mu*Nd] = 0.0;
    return res;
  }

public:
  INHERIT_GIMPL_TYPES(Gimpl)

  Smear_APE(const std::vector<double>& rho_):rho(rho_){}
  Smear_APE(double rho_val):rho(set_rho(rho_val)){}
  Smear_APE():rho(set_rho(1.0)){}
  ~Smear_APE(){}

  void smear(GaugeField& u_smr, const GaugeField& U)const{
    double d_rho;
    GaugeLinkField Cup, tmp_stpl;
    WilsonLoops<Gimpl> WL;
    u_smr = zero;

    for(int mu=0; mu<Nd; ++mu){
      Cup = zero;
      for(int nu=0; nu<Nd; ++nu){
	d_rho = rho[mu + Nd * nu];
	WL.Staple(tmp_stpl, U, mu, nu);
	Cup += tmp_stpl*d_rho;
      }
      pokeLorentz(u_smr, Cup, mu);
    }
  }

  void derivative(GaugeField& SigmaTerm,
		  const GaugeField& iLambda,
		  const GaugeField& U)const{

    /*


    WilsonLoops<Gimpl> WL;
    GaugeLinkField staple, u_tmp, iLambda_mu, iLambda_nu;
    GaugeLinkField U_mu, U_nu;
    GaugeLinkField sh_field ;
    GaugeLinkField temp_Sigma;

    SU<N>::Matrix temp_mat, temp_mat2;
    Real rho_munu, rho_numu;
    
    // to be completed 
    int Nvol = CommonPrms::instance()->Nvol();
    
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
	
	temp_Sigma = adj(staple)*iLambda_nu;
	temp_Sigma *= - rho_numu;
	//-r_numu*U_nu(x+mu)*Udag_mu(x+nu)*Udag_nu(x)*Lambda_nu(x)
	SigmaTerm .................

	for (int site = 0; site < Nvol; ++site){
	  temp_mat = mat_dag(staple,site) * mat(iLambda_nu,site);
	  temp_mat *= - rho_numu;
	  AddMat(SigmaTerm, temp_mat, site, mu);
	}
	sh_field = shiftField(iLambda_nu, mu, Forward());
	
	for (int site = 0; site < Nvol; ++site){
	  temp_mat = mat(sh_field,site) * mat_dag(staple,site);
	  temp_mat *= rho_numu;
	  AddMat(SigmaTerm, temp_mat, site, mu);
	}//r_numu*Lambda_nu(mu)*U_nu(x+mu)*Udag_mu(x+nu)*Udag_nu(x)
	
	sh_field = shiftField(iLambda_mu, nu, Forward());
	
	for (int site = 0; site < Nvol; ++site){
	  temp_mat = mat(U_nu,site) * mat(sh_field,site) * mat_dag(U_nu,site);
	  temp_mat = mat_dag(staple,site) * temp_mat;
	  temp_mat *= - rho_munu;
	  AddMat(SigmaTerm, temp_mat, site, mu);
	}//-r_munu*U_nu(x+mu)*Udag_mu(x+nu)*Lambda_mu(x+nu)*Udag_nu(x)
	
	staple = 0.0;
	sh_field = shiftField(U_nu, mu, Forward());
	
	for (int site = 0; site < Nvol; ++site){
	  temp_mat2 = mat_dag(sh_field,site) * mat_dag(U_mu,site);
	  temp_mat = temp_mat2  * mat(iLambda_mu,site) * mat(U_nu,site);
	  temp_mat *= - rho_munu;
	  AddMat(staple, temp_mat, site);
	  temp_mat = temp_mat2 * mat(iLambda_nu,site) * mat(U_nu,site);
	  temp_mat *= rho_numu;
	  AddMat(staple, temp_mat, site);
	} 
	
	for (int site = 0; site < Nvol; ++site){
	  temp_mat = mat_dag(U_nu,site) * mat(iLambda_nu,site);
	  SetMat(u_tmp, temp_mat, site);
	}
	
	sh_field = shiftField(u_tmp, mu, Forward()); 
	
	for (int site = 0; site < Nvol; ++site){
	  temp_mat = mat(sh_field,site) * mat_dag(U_mu,site) * mat(U_nu,site);
	  temp_mat *= - rho_numu;
	  AddMat(staple, temp_mat, site);
	}   
	
	sh_field = shiftField(staple, nu, Backward());
	
	AddSlice(SigmaTerm, sh_field, mu);
      }
    }

    */
  }




};

#endif  
