
/*!
  @file GaugeConfiguration.h
  @brief Declares the GaugeConfiguration class
*/
#pragma once

NAMESPACE_BEGIN(Grid);


template<class T> void Dump(const Lattice<T> & lat,
			    std::string s,
			    Coordinate site = Coordinate({0,0,0,0}))
{
  typename T::scalar_object tmp;
  if (lat.Checkerboard() != lat.Grid()->CheckerBoard(site))
    site = Coordinate({0,0,0,1});
  peekSite(tmp,lat,site);
  std::cout << " Dump "<<s<<" "<<tmp<<std::endl;
}


/*!
  @brief Smeared configuration masked container
  Modified for a multi-subset smearing (aka Luscher Flowed HMC)
*/
template <class Gimpl>
class SmearedConfigurationMasked : public SmearedConfiguration<Gimpl>
{
public:
  INHERIT_GIMPL_TYPES(Gimpl);

private:
  // These live in base class
  //  const unsigned int smearingLevels;
  //  Smear_Stout<Gimpl> *StoutSmearing;
  //  std::vector<GaugeField> SmearedSet;

  GridCartesian         * UGrid; // keep a copy of the grid 
  GridRedBlackCartesian * UrbGrid; // keep a copy of the redblack grid for life of object
  PaddedCell              Ghost;

  std::vector<LatticeLorentzComplex> masks;
  std::vector<int> cbs;
  std::vector<GeneralLocalStencil> gStencils;
  std::vector<GeneralLocalStencil> gStencils_smear;
  
  typedef typename SU3Adjoint::AMatrix AdjMatrix;
  typedef typename SU3Adjoint::LatticeAdjMatrix  AdjMatrixField;
  typedef typename SU3Adjoint::LatticeAdjVector  AdjVectorField;
  typedef typename SU3::vAlgebraMatrix vAlgebraMatrix;
  
  
  void BaseSmearDerivative(GaugeField& SigmaTerm,
			   const GaugeField& iLambda,
			   const GaugeField& U,
			   int mmu, RealD rho)
  {GRID_TRACE("BaseSmearDerivative");
    // Reference
    // Morningstar, Peardon, Phys.Rev.D69,054501(2004)
    // Equation 75
    // Computing Sigma_mu, derivative of S[fat links] with respect to the thin links
    // Output SigmaTerm

    GridBase *grid = U.Grid();

    WilsonLoops<Gimpl> WL;
    GaugeLinkField staple(grid), u_tmp(grid);
    GaugeLinkField iLambda_mu(grid), iLambda_nu(grid);
    GaugeLinkField U_mu(grid), U_nu(grid);
    GaugeLinkField sh_field(grid), temp_Sigma(grid);
    Real rho_munu, rho_numu;

    RealD t = 0;
    t-=usecond();
    
    rho_munu = rho;
    rho_numu = rho;
    for(int mu = 0; mu < Nd; ++mu){
      U_mu       = peekLorentz(      U, mu);
      iLambda_mu = peekLorentz(iLambda, mu);

      for(int nu = 0; nu < Nd; ++nu){
	if(nu==mu) continue;

	U_nu       = peekLorentz(      U, nu);

	// Nd(nd-1) = 12 staples normally.
	// We must compute 6 of these
	// in FTHMC case
	if ( (mu==mmu)||(nu==mmu) )
	  WL.StapleUpper(staple, U, mu, nu);
	
	if(nu==mmu) {
	  iLambda_nu = peekLorentz(iLambda, nu);

	  temp_Sigma = -rho_numu*staple*iLambda_nu;  //ok
	  //-r_numu*U_nu(x+mu)*Udag_mu(x+nu)*Udag_nu(x)*Lambda_nu(x)
	  Gimpl::AddLink(SigmaTerm, temp_Sigma, mu);

	  sh_field = Cshift(iLambda_nu, mu, 1);// general also for Gparity?

	  temp_Sigma = rho_numu*sh_field*staple; //ok
	  //r_numu*Lambda_nu(mu)*U_nu(x+mu)*Udag_mu(x+nu)*Udag_nu(x)
	  Gimpl::AddLink(SigmaTerm, temp_Sigma, mu);
	}

	if ( mu == mmu ) { 
	  sh_field = Cshift(iLambda_mu, nu, 1);

	  temp_Sigma = -rho_munu*staple*U_nu*sh_field*adj(U_nu); //ok
	  //-r_munu*U_nu(x+mu)*Udag_mu(x+nu)*Lambda_mu(x+nu)*Udag_nu(x)
	  Gimpl::AddLink(SigmaTerm, temp_Sigma, mu);
	}

	//	staple = Zero();
	sh_field = Cshift(U_nu, mu, 1);

	temp_Sigma = Zero();

	if ( mu == mmu )
	  temp_Sigma = -rho_munu*adj(sh_field)*adj(U_mu)*iLambda_mu*U_nu;

	if ( nu == mmu ) {
	  temp_Sigma += rho_numu*adj(sh_field)*adj(U_mu)*iLambda_nu*U_nu;

	  u_tmp = adj(U_nu)*iLambda_nu;
	  sh_field = Cshift(u_tmp, mu, 1);
	  temp_Sigma += -rho_numu*sh_field*adj(U_mu)*U_nu;
	}
	
	sh_field = Cshift(temp_Sigma, nu, -1);
	Gimpl::AddLink(SigmaTerm, sh_field, mu);

      }
    }

    t+=usecond();
    std::cout << GridLogMessage << " BaseSmearDerivative " << t/1e3 << " ms " << std::endl;
    
  }
  //for debugging
  void BaseSmear_cb(GaugeLinkField& Cup, const GaugeField& U,int mu,RealD rho) {
    GRID_TRACE("BaseSmear_cb");
    GridBase *grid = U.Grid();
    GridBase *hgrid = Cup.Grid();
    GaugeLinkField tmp_stpl(grid);
    GaugeLinkField tmp_stpl_eo(hgrid);
    WilsonLoops<Gimpl> WL;
    int cb = Cup.Checkerboard();
    RealD t = 0;

    t-=usecond();
    Cup = Zero();
    for(int nu=0; nu<Nd; ++nu){
      if (nu != mu) {
        // get the staple in direction mu, nu
        WL.Staple(tmp_stpl, U, mu, nu);  //nb staple conventions of IroIro and Grid differ by a dagger
        pickCheckerboard(cb,tmp_stpl_eo,tmp_stpl); // ideally, compute tmp_stpl only on the current checkerboard
        Cup += adj(tmp_stpl_eo*rho);
      }
    }
    t+=usecond();
    std::cout << GridLogMessage << " BaseSmear " << t/1e3 << " ms " << std::endl;
  }
  void BaseSmear(GaugeLinkField& Cup, const GaugeField& U,int mu,RealD rho) {
    GRID_TRACE("BaseSmear");
    GridBase *grid = U.Grid();
    GaugeLinkField tmp_stpl(grid);
    WilsonLoops<Gimpl> WL;
    RealD t = 0;

    t-=usecond();
    Cup = Zero();
    for(int nu=0; nu<Nd; ++nu){
      if (nu != mu) {
	// get the staple in direction mu, nu
	WL.Staple(tmp_stpl, U, mu, nu);  //nb staple conventions of IroIro and Grid differ by a dagger
	Cup += adj(tmp_stpl*rho);
      }
    }
    t+=usecond();
    std::cout << GridLogMessage << " BaseSmear " << t/1e3 << " ms " << std::endl;
  }
  
  // Assume: gU is extended gauge field with its boundary of size 1
  void BaseSmear_ghost(GaugeLinkField& Cup, const GaugeField& gU,int mu,RealD rho) {
    GRID_TRACE("BaseSmear_ghost");
    GridBase *ggrid = gU.Grid();
    GaugeLinkField tmp_stpl(ggrid);
    int cb = Cup.Checkerboard();
    RealD t = 0;

    t-=usecond();
    autoView( tmp_stpl_v , tmp_stpl, AcceleratorWrite);
    autoView( gU_v , gU, AcceleratorRead);
    autoView( gStencil_v, gStencils_smear[mu], AcceleratorRead);
    accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
	typedef decltype(coalescedRead(tmp_stpl_v[0])) LinkMat;

	LinkMat tmp = Zero();
	for(int nu=0; nu<Nd; ++nu){
	  int inc = 6*(nu - (mu<=nu));
	  if (nu != mu) {
	    GeneralStencilEntry const* e = gStencil_v.GetEntry(0+inc,ss);
	    auto U_nu_x = coalescedReadGeneralPermute(gU_v[e->_offset], e->_permute, Nd)(nu)();
	    e = gStencil_v.GetEntry(1+inc,ss);
	    auto U_mu_xpnu = coalescedReadGeneralPermute(gU_v[e->_offset], e->_permute, Nd)(mu)();
	    e = gStencil_v.GetEntry(2+inc,ss);
	    auto Udag_nu_xpmu = adj(coalescedReadGeneralPermute(gU_v[e->_offset], e->_permute, Nd))(nu)();
      
	    tmp()() = tmp()() + U_nu_x * U_mu_xpnu * Udag_nu_xpmu;

	    e = gStencil_v.GetEntry(3+inc,ss);
	    auto Udag_nu_xmnu = adj(coalescedReadGeneralPermute(gU_v[e->_offset], e->_permute, Nd))(nu)();
	    e = gStencil_v.GetEntry(4+inc,ss);
	    auto U_mu_xmnu = coalescedReadGeneralPermute(gU_v[e->_offset], e->_permute, Nd)(mu)();
	    e = gStencil_v.GetEntry(5+inc,ss);
	    auto U_nu_xpmu_mnu = coalescedReadGeneralPermute(gU_v[e->_offset], e->_permute, Nd)(nu)();

	    tmp()() = tmp()() + Udag_nu_xmnu * U_mu_xmnu * U_nu_xpmu_mnu;
            
	  }
	}
	coalescedWrite(tmp_stpl_v[ss],rho*tmp);
      });
    pickCheckerboard(cb,Cup,Ghost.Extract(tmp_stpl));
    
    t+=usecond();
    std::cout << GridLogMessage << " BaseSmear_ghost " << t/1e3 << " ms " << std::endl;
  }

  // Adjoint vector to GaugeField force
  void InsertForce(GaugeField &Fdet,AdjVectorField &Fdet_nu,int nu)
  {
    Complex ci(0,1);
    GaugeLinkField Fdet_pol(Fdet.Grid());
    RealD t = 0;

    t-=usecond();
    Fdet_pol=Zero();
    for(int e=0;e<8;e++){
      ColourMatrix te;
      SU3::generator(e, te);
      auto tmp=peekColour(Fdet_nu,e);
      Fdet_pol=Fdet_pol + ci*tmp*te; // but norm of te is different.. why?
    }
    pokeLorentz(Fdet, Fdet_pol, nu);

    t+=usecond();
    std::cout << GridLogPerformance << " InsertForce " << t/1e3 << " ms " << std::endl;
  }

  void Compute_MpInvJx_dNxxdSy(const GaugeLinkField &PlaqL,const GaugeLinkField &PlaqR, AdjMatrixField MpInvJx,AdjVectorField &Fdet2 )
  {
    GRID_TRACE("Compute_MpInvJx_dNxxdSy");
    int cb = PlaqL.Checkerboard();
    GridBase *grid = PlaqL.Grid();
    const int Ngen = SU3Adjoint::Dimension;
    Complex ci(0,1);
    RealD t;

    t=-usecond();
    autoView(Fdet2_v,Fdet2,AcceleratorWrite);
    autoView(PlaqL_v,PlaqL,AcceleratorRead);
    autoView(PlaqR_v,PlaqR,AcceleratorRead);
    autoView(MpInvJx_v,MpInvJx,AcceleratorRead);
    const int nsimd = vAlgebraMatrix::Nsimd();
    accelerator_for(ss,grid->oSites(),nsimd,{
	typedef decltype(coalescedRead(MpInvJx_v[0])) adj_mat;
	typedef decltype(coalescedRead(Fdet2_v[0]))   adj_vec;
	
	adj_mat Dbc;
	adj_vec Fdet;
	ColourMatrix ta;
	
	for(int a=0;a<Ngen;a++) {
	  // Qlat Tb = 2i Tb^Grid
	  SU3::generator(a, ta);
	  ta = 2.0 * ci * ta;
	  auto UtaU = adj(PlaqL_v(ss))*ta*PlaqR_v(ss);
	  SU3::LieAlgebraProject(Dbc,UtaU);
	  Fdet()()(a) = traceProduct(MpInvJx_v(ss),Dbc)()()();
	}
	coalescedWrite(Fdet2_v[ss],Fdet);
      });
    t+=usecond();
    std::cout << GridLogPerformance << " Compute_MpInvJx_dNxxdSy " << t/1e3 <<" ms"<<std::endl;
  }

  void ComputeNxy(const GaugeLinkField &PlaqL,const GaugeLinkField &PlaqR,AdjMatrixField &NxAd)
  {
    GRID_TRACE("ComputeNxy");
    GridBase *grid = PlaqL.Grid();
    RealD t = 0;

    t-=usecond();
    autoView(NxAd_v,NxAd,AcceleratorWrite);
    autoView(PlaqL_v,PlaqL,AcceleratorRead);
    autoView(PlaqR_v,PlaqR,AcceleratorRead);
    const int nsimd = vAlgebraMatrix::Nsimd();
    accelerator_for(ss,grid->oSites(),nsimd,{
        typedef decltype(coalescedRead(NxAd_v[0]))  adj_mat;
        adj_mat NxAd_site;
	SU3::LieAlgebraProject(NxAd_site,PlaqL_v(ss),PlaqR_v(ss));
        coalescedWrite(NxAd_v[ss],NxAd_site);
      });
    t+=usecond();
    std::cout << GridLogPerformance << " ComputeNxy " << t/1e3 <<" ms"<<std::endl;
  }

  void ApplyMask(GaugeField &U,int smr)
  {
    LatticeComplex tmp(U.Grid());
    GaugeLinkField Umu(U.Grid());
    for(int mu=0;mu<Nd;mu++){
      Umu=PeekIndex<LorentzIndex>(U,mu);
      tmp=PeekIndex<LorentzIndex>(masks[smr],mu);
      Umu=Umu*tmp;
      PokeIndex<LorentzIndex>(U, Umu, mu);
    }
  }
public:

  void logDetJacobianForceLevel(const GaugeField &U, GaugeField &force ,int smr)
  {
    GRID_TRACE("logDetJacobianForceLevel");
    GridBase* grid = U.Grid();
    GridBase* hgrid = UrbGrid; // For now, assume masking is based on red-black checkerboarding
    assert(grid==UGrid);
    //GaugeField C(grid);
    GaugeField Umsk(grid);
    std::vector<GaugeLinkField> Umu(Nd,grid);
    GaugeLinkField Cmu(hgrid); // U and staple; C contains factor of epsilon
    GaugeLinkField Zx(hgrid);  // U times Staple, contains factor of epsilon
    GaugeLinkField Utmp(grid);
    GaugeLinkField Ueo(hgrid);
    GaugeLinkField PlaqL(hgrid);
    GaugeLinkField PlaqR(hgrid);
    const int Ngen = SU3Adjoint::Dimension;
    ColourMatrix Ident;
    
    AdjVectorField  dJdXe_nMpInv(hgrid); 
    AdjVectorField  dJdXe_nMpInv_y(hgrid); 
    AdjMatrixField  MpAd(hgrid);    // Mprime luchang's notes
    AdjMatrixField  MpAdInv(hgrid); // Mprime inverse
    AdjMatrixField  NxxAd(hgrid);    // Nxx in adjoint space
    AdjMatrixField  JxAd(hgrid);     
    AdjMatrixField  ZxAd(hgrid);

    RealD t0 = usecond();
    Ident = ComplexD(1.0);
    for(int d=0;d<Nd;d++){
      Umu[d] = peekLorentz(U, d);
    }
    int mu= (smr/2) %Nd;

    ////////////////////////////////////////////////////////////////////////////////
    // Mask the gauge field
    ////////////////////////////////////////////////////////////////////////////////
    int cb = cbs[smr];
    auto mask=PeekIndex<LorentzIndex>(masks[smr],mu); // the cb mask

    Umsk = U;
    ApplyMask(Umsk,smr);
    Utmp = peekLorentz(Umsk,mu);
    pickCheckerboard(cb,Ueo,Utmp);
    
    GaugeField gU(grid);
    std::vector<GaugeLinkField> gUmu(Nd,grid);
    {
      GRID_TRACE("ExchangePeriodic");
      gU = Ghost.ExchangePeriodic(U);
      for(int d=0; d<Nd;d++)
	gUmu[d] = peekLorentz(gU, d);
    }
    GridBase       *ggrid = gUmu[0].Grid();
    
    Cmu.Checkerboard() = cb;
    Zx.Checkerboard() = cb;
    Ueo.Checkerboard() = cb;
    PlaqL.Checkerboard() = cb;
    PlaqR.Checkerboard() = cb;
    MpAd.Checkerboard() = cb;
    MpAdInv.Checkerboard() = cb;
    dJdXe_nMpInv.Checkerboard() = cb;
    dJdXe_nMpInv_y.Checkerboard() = cb;
    NxxAd.Checkerboard() = cb;
    JxAd.Checkerboard() = cb;
    ZxAd.Checkerboard() = cb;
    
    ////////////////////////////////////////////////////////////////////////////////
    // Retrieve the eps/rho parameter(s) -- could allow all different but not so far
    ////////////////////////////////////////////////////////////////////////////////
    double rho=this->StoutSmearing->SmearRho[1];
    int idx=0;
    for(int mu=0;mu<4;mu++){
      for(int nu=0;nu<4;nu++){
	if ( mu!=nu) assert(this->StoutSmearing->SmearRho[idx]==rho);
	else         assert(this->StoutSmearing->SmearRho[idx]==0.0);
	idx++;
      }}
    
    //////////////////////////////////////////////////////////////////
    // Assemble the N matrix
    //////////////////////////////////////////////////////////////////
    // Computes ALL the staples -- could compute one only and do it here
    RealD time;
    time=-usecond();
    BaseSmear_ghost(Cmu, gU, mu, rho);
#if 0 // DEBUG
    GaugeLinkField Cmu2(hgrid);
    Cmu2.Checkerboard() = cb;
    BaseSmear_cb(Cmu2, U, mu, rho);
    std::cout << GridLogMessage << " DEBUG: BaseSmear " <<smr<<" "<<mu<<" "<<cb<<" "<<" simd "<<AdjMatrix::Nsimd()<<" "<<norm2(Cmu-Cmu2)<<std::endl;
#endif
    
    //////////////////////////////////////////////////////////////////
    // Assemble Luscher exp diff map J matrix 
    //////////////////////////////////////////////////////////////////
    // Ta so Z lives in Lie algabra
    {GRID_TRACE("Zx");
    Zx  = Ta(Cmu * adj(Ueo));
    }
    time+=usecond();
    std::cout << GridLogMessage << "Z took "<<time<< " us"<<std::endl;
    
    // Move Z to the Adjoint Rep == make_adjoint_representation
    time=-usecond();
    {GRID_TRACE("ZxAd");
    SU3Adjoint::make_adjoint_rep(ZxAd, Zx);
    }
    time+=usecond();
    std::cout << GridLogPerformance << "ZxAd took "<<time<< " us"<<std::endl;
#if 0 //DEBUG
    AdjMatrixField  ZxAd2(hgrid); ZxAd2.Checkerboard() = cb;
    LatticeComplex  cplx(hgrid);  cplx.Checkerboard() = cb;
    AdjMatrix TRb;
    Complex ci(0,1);
    ColourMatrix   tb;

    ZxAd2 = Zero();
    for(int b=0;b<8;b++) {
      // Adj group sets traceless antihermitian T's -- Guido, really????
      SU3::generator(b, tb);         // Fund group sets traceless hermitian T's
      SU3Adjoint::generator(b,TRb);
      TRb=-TRb;
      cplx = 2.0*trace(ci*tb*Zx); // my convention 1/2 delta ba
      ZxAd2 = ZxAd2 + cplx * TRb; // is this right? YES - Guido used Anti herm Ta's and with bloody wrong sign.
      
    }
    std::cout << GridLogMessage << " DEBUG: ZxAd " <<smr<<" "<<mu<<" "<<cb<<" "<<" simd "<<AdjMatrix::Nsimd()<<" "<<norm2(ZxAd-ZxAd2)<<" "<<norm2(ZxAd)<<" "<<norm2(ZxAd2)<<std::endl;
#endif
    
    //////////////////////////////////////
    // J(x) = 1 + Sum_k=1..N (-Zac)^k/(k+1)!
    //////////////////////////////////////
    time=-usecond();
    {GRID_TRACE("JxAd");
    autoView(JxAd_v,JxAd,AcceleratorWrite);
    autoView(ZxAd_v,ZxAd,AcceleratorRead);
    const int nsimd = vAlgebraMatrix::Nsimd();
    accelerator_for(ss,hgrid->oSites(),nsimd,{
        typedef decltype(coalescedRead(JxAd_v[0])) adj_mat;
        adj_mat X, JxAd_site;
	RealD kpfac = 1;
	
	X=1.0;
	JxAd_site = X;
	for(int k=1;k<12;k++){
	  X=-X*ZxAd_v(ss);
	  kpfac = kpfac /(k+1);
	  JxAd_site = JxAd_site + X * kpfac;
	}
        coalescedWrite(JxAd_v[ss],JxAd_site);
      });
    }
    time+=usecond();
    std::cout << GridLogMessage << "Jx took "<<time<< " us"<<std::endl;
    
    /////////////////////////////////////////////////////////////////
    // NxxAd
    /////////////////////////////////////////////////////////////////
    time=-usecond();
    {
      GRID_TRACE("Plaq_for_Nxy");
    PlaqL = Ident;
    PlaqR = Ueo*adj(Cmu);
    }
    ComputeNxy(PlaqL,PlaqR,NxxAd);
    time+=usecond();
    std::cout << GridLogMessage << "ComputeNxy took "<<time<< " us"<<std::endl;
#if 0 // DEBUG
    AdjMatrixField  NxxAd2(hgrid); NxxAd2.Checkerboard() = cb;
    ComputeNxy(0,PlaqL,PlaqR,NxxAd2);
    std::cout << GridLogMessage << " DEBUG: NxxAd " <<smr<<" "<<mu<<" "<<cb<<" "<<" simd "<<AdjMatrix::Nsimd()<<" "<<norm2(NxxAd-NxxAd2)<<std::endl;
#endif
    ////////////////////////////
    // Mab
    ////////////////////////////
    MpAd = Complex(1.0,0.0);
    MpAd = MpAd - JxAd * NxxAd; 

    /////////////////////////
    // invert the 8x8
    /////////////////////////
    {GRID_TRACE("MpAdInv");
    time=-usecond();
    MpAdInv = Inverse_RealPart(MpAd);
    MpAdInv.Checkerboard() = cb; //inside, it calls Lattice(GridBase *grid,ViewMode mode=AcceleratorWriteDiscard) in Lattice_base.h & sets checkerboard to 0
    time+=usecond();
    }
    std::cout << GridLogPerformance << "MpAdInv took "<<time<< " us"<<std::endl;
    
    /////////////////////////////////////////////////////////////////
    // M'^{-1} J(X)_ad
    /////////////////////////////////////////////////////////////////
    AdjMatrixField MpInvJx(hgrid);      MpInvJx.Checkerboard() = cb; 
    AdjMatrixField MpInvJx_nu(hgrid);   MpInvJx_nu.Checkerboard() = cb;
    {GRID_TRACE("Mp_inv_Jx");
    MpInvJx = (-1.0)*MpAdInv * JxAd;// rho is on the plaq factor
    }

    RealD t3a = usecond();
    /////////////////////////////////////////////////////////////////
    // dJ(x)/dxe N M'^{-1}  
    /////////////////////////////////////////////////////////////////
    time=-usecond();
    {GRID_TRACE("dJdX_nMpinv_combined");

    autoView(dJdXe_nMpInv_v,dJdXe_nMpInv,AcceleratorWrite);
    autoView(ZxAd_v,ZxAd,AcceleratorRead);
    autoView(NxxAd_v,NxxAd,AcceleratorRead);
    autoView(MpAdInv_v,MpAdInv,AcceleratorRead);
    const int nsimd = vAlgebraMatrix::Nsimd();
    accelerator_for(ss,hgrid->oSites(),nsimd,{
	typedef decltype(coalescedRead(ZxAd_v[0]))         adj_mat;
	typedef decltype(coalescedRead(dJdXe_nMpInv_v[0])) adj_vec;
      	adj_mat X, t2, dt2, t3, dt3, aunit, nMpInv_site;
	adj_vec dJdXe_nMpInv_site;
	iVector<adj_mat,Ngen> iTas;
	iVector<adj_mat,Ngen> dJdX_b;

	for(int b=0;b<Ngen;b++){
	  SU3Adjoint::generator(b, iTas(b));
	  dJdX_b(b) = iTas(b);
	}
	aunit = ComplexD(1.0);
	X  = (-1.0)*ZxAd_v(ss);
	t2 = X;
	for (int j = 12; j > 1; --j) {
	  t3  = t2*(1.0 / (j + 1))  + aunit;
	  t2  = X * t3;
	  for(int b=0;b<Ngen;b++){
	    dJdX_b(b)= iTas(b) * t3 + X * dJdX_b(b)*(1.0 / (j + 1));
	  }
	}
	nMpInv_site= NxxAd_v(ss) * MpAdInv_v(ss);
	for(int e=0;e<Ngen;e++){
	  dJdXe_nMpInv_site()()(e) = traceProduct((-1.0)*dJdX_b(e),nMpInv_site)()()();
	}
	coalescedWrite(dJdXe_nMpInv_v[ss],dJdXe_nMpInv_site);
      });
    }//make sure view object is closed
    time += usecond();
    std::cout << GridLogMessage << "dJdX_nMpinv_combined took "<<time<< " us"<<std::endl;

    RealD t3b = usecond();
    
    AdjVectorField  Fdet1_mu(grid);
    AdjVectorField  Fdet2_mu(grid);
    AdjVectorField  Fdet1_nu(grid);
    AdjVectorField  Fdet2_nu(grid);

    AdjVectorField  FdetV(hgrid);       FdetV.Checkerboard() = cb;
    AdjVectorField  Fdet1_nu_eo(hgrid); Fdet1_nu_eo.Checkerboard() = cb;       
    AdjVectorField  Fdet1_nu_oe(hgrid); Fdet1_nu_oe.Checkerboard() = (cb+1)%2; 
    AdjVectorField  Fdet2_nu_eo(hgrid); Fdet2_nu_eo.Checkerboard() = cb;       
    AdjVectorField  Fdet2_nu_oe(hgrid); Fdet1_nu_oe.Checkerboard() = (cb+1)%2;
    
    AdjVectorField  Fdet1_mu_oe(hgrid); Fdet1_mu_oe.Checkerboard() = (cb+1)%2; Fdet1_mu_oe = Zero();
    AdjVectorField  Fdet2_mu_oe(hgrid); Fdet2_mu_oe.Checkerboard() = (cb+1)%2; Fdet2_mu_oe = Zero();

    // Set even part of Fdet1_mu & Fdet2_mu <- cb is referred to as even here, regardless of actual parity corredponding to cb
    setCheckerboard(Fdet1_mu, (AdjVectorField) (transpose(NxxAd)*dJdXe_nMpInv)); 
    Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx,FdetV);
    setCheckerboard(Fdet2_mu,FdetV);

    RealD t3c = usecond();
    
    //    dJdXe_nMpInv needs to multiply:
    //       Nxx_mu (site local)                           (1)
    //       Nxy_mu one site forward  in each nu direction (3)
    //       Nxy_mu one site backward in each nu direction (3)
    //       Nxy_nu 0,0  ; +mu,0; 0,-nu; +mu-nu   [ 3x4 = 12]
    // 19 terms.
    AdjMatrixField Nxy(hgrid);
    // force = Fdet1 + Fdet2
    GaugeField Fdet1(grid);
    GaugeField Fdet2(grid);

    RealD t4 = usecond(), tLR = 0, tNxy = 0, tMJx = 0, t_ins=0, t_ck = 0, t_stencil=0, t_cshift=0;

    GaugeLinkField  gPlaqL(ggrid), gPlaqR(ggrid);
    
    autoView( gPlaqL_v , gPlaqL, AcceleratorWrite);
    autoView( gPlaqR_v , gPlaqR, AcceleratorWrite);
    autoView( gU_mu_v , gUmu[mu], AcceleratorRead);
    
    for(int nu=0;nu<Nd;nu++){
      
      if (nu!=mu) {
	GRID_TRACE("MuNuLoopBody");
	autoView( gU_nu_v , gUmu[nu], AcceleratorRead);


	///////////////// +ve nu /////////////////
	//     __
	//    |  |
	//    x==    // nu polarisation -- clockwise

	time=-usecond(); tLR -= usecond();

	PlaqL=Ident;
	{
	  GRID_TRACE("Staple");
	autoView( gStencil_v  , gStencils[mu*(Nd-1)*6+(nu-(mu<=nu))*6], AcceleratorRead);
	accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
          GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
          auto U_nu_x = coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
          e = gStencil_v.GetEntry(1,ss);
          auto U_mu_xpnu = coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd);
          e = gStencil_v.GetEntry(2,ss);
          auto Udag_nu_xpmu = adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
	  e = gStencil_v.GetEntry(3,ss);
          auto Udag_mu_x = adj(coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd));

          auto stencil_ss = (-rho) * U_nu_x * U_mu_xpnu * Udag_nu_xpmu * Udag_mu_x;

          coalescedWrite(gPlaqR_v[ss],stencil_ss);
	}
        );
	t_ck -= usecond();
	pickCheckerboard(cb,PlaqR,Ghost.Extract(gPlaqR));
	t_ck += usecond();
	}
	time+=usecond(); tLR += usecond();
	std::cout << GridLogMessage << "PlaqLR took "<<time<< " us "<<" checkerboard_extract "<<t_ck<<" us"<<std::endl;

	time=-usecond(); tNxy -= usecond();
	PlaqL.Checkerboard() = cb; Nxy.Checkerboard() = cb; FdetV.Checkerboard() = cb;
	
	dJdXe_nMpInv_y = dJdXe_nMpInv;
	ComputeNxy(PlaqL,PlaqR,Nxy);
	
	Fdet1_nu_eo = transpose(Nxy)*dJdXe_nMpInv_y;
	time+=usecond(); tNxy += usecond();
	std::cout << GridLogMessage << "ComputeNxy (occurs 6x) took "<<time<< " us"<<std::endl;

	time=-usecond(); tMJx -= usecond();
	PlaqR=(-1.0)*PlaqR;
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx,FdetV);
	Fdet2_nu_eo = FdetV;
	time+=usecond(); tMJx += usecond();
	std::cout << GridLogMessage << "Compute_MpInvJx_dNxxSy (occurs 6x) took "<<time<< " us"<<std::endl;
	
	//    x==
	//    |  |
	//    .__|    // nu polarisation -- anticlockwise

	tLR -= usecond();
	{
	  GRID_TRACE("Staple");
	autoView( gStencil_v  , gStencils[mu*(Nd-1)*6+(nu-(mu<=nu))*6+1], AcceleratorRead);
        accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
          GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
          auto U_nu_y = coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
          e = gStencil_v.GetEntry(1,ss);
          auto Udag_mu_ymmupnu = adj(coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd));
          e = gStencil_v.GetEntry(2,ss);
          auto Udag_nu_ymmu = adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
	  auto stencil_ss = (rho) * U_nu_y * Udag_mu_ymmupnu * Udag_nu_ymmu;
	  coalescedWrite(gPlaqR_v[ss],stencil_ss);
	  
          e = gStencil_v.GetEntry(3,ss);
          auto Udag_mu_ymmu = adj(coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd));
	  stencil_ss = Udag_mu_ymmu;
          coalescedWrite(gPlaqL_v[ss],stencil_ss);
        }
        );
	pickCheckerboard((cb+1)%2,PlaqR,Ghost.Extract(gPlaqR));
	pickCheckerboard((cb+1)%2,PlaqL,Ghost.Extract(gPlaqL));
	}
	tLR += usecond();

	tNxy -= usecond();
	Nxy.Checkerboard() = (cb+1)%2; 	FdetV.Checkerboard() = (cb+1)%2;
	t_cshift-= usecond();
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,mu,-1);
	t_cshift += usecond();
	ComputeNxy(PlaqL, PlaqR,Nxy);
	Fdet1_nu_oe = transpose(Nxy)*dJdXe_nMpInv_y;
	tNxy += usecond();

	tMJx -= usecond();
	MpInvJx_nu = Cshift(MpInvJx,mu,-1);
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_nu_oe = FdetV;
	tMJx += usecond();
	
	///////////////// -ve nu /////////////////
	//  __
	// |  |
	// x==          // nu polarisation -- clockwise
	tLR -= usecond();
	{ GRID_TRACE("Staple");
	autoView( gStencil_v  , gStencils[mu*(Nd-1)*6+(nu-(mu<=nu))*6+2], AcceleratorRead);
        accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
          GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
          auto U_mu_y = coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd);
          e = gStencil_v.GetEntry(1,ss);
          auto U_nu_ypmu = coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
          e = gStencil_v.GetEntry(2,ss);
          auto Udag_mu_ypnu = adj(coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd));

          auto stencil_ss = (rho) * U_mu_y * U_nu_ypmu * Udag_mu_ypnu;
          coalescedWrite(gPlaqL_v[ss],stencil_ss);
        }
        );
	pickCheckerboard((cb+1)%2,PlaqL,Ghost.Extract(gPlaqL)); 
        pickCheckerboard((cb+1)%2,PlaqR,Umu[nu]);
	}
	tLR += usecond();

	tNxy -= usecond();
	t_cshift-=usecond();
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,1);
	t_cshift+=usecond();
	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_nu_oe = Fdet1_nu_oe + transpose(Nxy)*dJdXe_nMpInv_y;
	tNxy += usecond();

	tMJx -= usecond();
	t_cshift-=usecond();
	MpInvJx_nu = Cshift(MpInvJx,nu,1);
	t_cshift+=usecond();
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_nu_oe = Fdet2_nu_oe+FdetV;
	tMJx += usecond();
      
	// x==
	// |  |
	// |__|         // nu polarisation
	tLR -= usecond();
	{
	  GRID_TRACE("Staple");
	autoView( gStencil_v  , gStencils[mu*(Nd-1)*6+(nu-(mu<=nu))*6+3], AcceleratorRead);
        accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
          GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
          auto U_nu_y = coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
          e = gStencil_v.GetEntry(1,ss);
          auto Udag_mu_ymmupnu = adj(coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd));
          auto stencil_ss = (-rho) * U_nu_y * Udag_mu_ymmupnu;
          coalescedWrite(gPlaqL_v[ss],stencil_ss);
	  
          e = gStencil_v.GetEntry(2,ss);
          auto Udag_mu_ymmu = adj(coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd));
          e = gStencil_v.GetEntry(3,ss);
          auto U_nu_ymmu = coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
          stencil_ss = Udag_mu_ymmu * U_nu_ymmu;
          coalescedWrite(gPlaqR_v[ss],stencil_ss);
        }
        );
	pickCheckerboard(cb,PlaqL,Ghost.Extract(gPlaqL));
        pickCheckerboard(cb,PlaqR,Ghost.Extract(gPlaqR));
	}
	tLR += usecond();

	tNxy -= usecond();
	Nxy.Checkerboard() = cb;  FdetV.Checkerboard() = cb;
	t_cshift-=usecond();
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,mu,-1);
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv_y,nu,1);
	t_cshift+=usecond();
	
	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_nu_eo = Fdet1_nu_eo + transpose(Nxy)*dJdXe_nMpInv_y;
	tNxy += usecond();

	tMJx -= usecond();
	t_cshift-=usecond();
	MpInvJx_nu = Cshift(MpInvJx,mu,-1);
	MpInvJx_nu = Cshift(MpInvJx_nu,nu,1);
	t_cshift+=usecond();
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_nu_eo = Fdet2_nu_eo+FdetV;
	tMJx += usecond();
	
	/////////////////////////////////////////////////////////////////////
	// Set up the determinant force contribution in 3x3 algebra basis
	/////////////////////////////////////////////////////////////////////
	t_ins -= usecond();
	setCheckerboard(Fdet1_nu, Fdet1_nu_eo);
	setCheckerboard(Fdet1_nu, Fdet1_nu_oe);
	InsertForce(Fdet1,Fdet1_nu,nu);
	setCheckerboard(Fdet2_nu, Fdet2_nu_eo);
        setCheckerboard(Fdet2_nu, Fdet2_nu_oe);
	InsertForce(Fdet2,Fdet2_nu,nu);
	t_ins+= usecond();
	
	//////////////////////////////////////////////////
	// Parallel direction terms
	//////////////////////////////////////////////////

	Nxy.Checkerboard() = (cb+1)%2; FdetV.Checkerboard() = (cb+1)%2;
	
        //     __
	//    |  "
	//    |__"x    // mu polarisation
	tLR -= usecond();
	{
	  GRID_TRACE("Staple");
	autoView( gStencil_v  , gStencils[mu*(Nd-1)*6+(nu-(mu<=nu))*6+4], AcceleratorRead);
        accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
          GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
          auto U_mu_y = coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd);
          e = gStencil_v.GetEntry(1,ss);
          auto Udag_nu_ypmumnu = adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
          e = gStencil_v.GetEntry(2,ss);
          auto Udag_mu_ymnu = adj(coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd));
          auto stencil_ss = (-rho) * U_mu_y * Udag_nu_ypmumnu * Udag_mu_ymnu;
          coalescedWrite(gPlaqL_v[ss],stencil_ss);

          e = gStencil_v.GetEntry(3,ss);
          auto Udag_nu_ymnu = adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
          stencil_ss = Udag_nu_ymnu;
          coalescedWrite(gPlaqR_v[ss],stencil_ss);
        }
        );
	pickCheckerboard((cb+1)%2,PlaqL,Ghost.Extract(gPlaqL));
        pickCheckerboard((cb+1)%2,PlaqR,Ghost.Extract(gPlaqR));
	}
	tLR += usecond();

	tNxy -= usecond();
	t_cshift-=usecond();
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,-1);
	t_cshift+=usecond();
	
	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_mu_oe = Fdet1_mu_oe + transpose(Nxy)*dJdXe_nMpInv_y;
	tNxy += usecond();

	tMJx -= usecond();
	t_cshift-=usecond();
	MpInvJx_nu = Cshift(MpInvJx,nu,-1);
	t_cshift+=usecond();
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_mu_oe = Fdet2_mu_oe+FdetV;
	tMJx += usecond();

	//  __
	// "  |
	// x__|          // mu polarisation
	tLR -= usecond();
	{
	  GRID_TRACE("Staple");
	autoView( gStencil_v  , gStencils[mu*(Nd-1)*6+(nu-(mu<=nu))*6+5], AcceleratorRead);
        accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
          GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
          auto U_mu_y = coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd);
          e = gStencil_v.GetEntry(1,ss);
          auto U_nu_ypmu = coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
          e = gStencil_v.GetEntry(2,ss);
          auto Udag_mu_ypnu = adj(coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd));

          auto stencil_ss = (-rho) * U_mu_y * U_nu_ypmu * Udag_mu_ypnu;
          coalescedWrite(gPlaqL_v[ss],stencil_ss);
        }
        );
        pickCheckerboard((cb+1)%2,PlaqL,Ghost.Extract(gPlaqL));
	pickCheckerboard((cb+1)%2,PlaqR,Umu[nu]);
	}
	tLR += usecond();

	tNxy -= usecond();
	t_cshift-=usecond();
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,1);
	t_cshift+=usecond();
	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_mu_oe = Fdet1_mu_oe + transpose(Nxy)*dJdXe_nMpInv_y;
	tNxy += usecond();

	tMJx -= usecond();
	t_cshift-=usecond();
	MpInvJx_nu = Cshift(MpInvJx,nu,1);
	t_cshift+=usecond();
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_mu_oe = Fdet2_mu_oe+FdetV;
	tMJx += usecond();
	
      }
    }
    RealD t5 = usecond();
    setCheckerboard(Fdet1_mu, Fdet1_mu_oe);
    setCheckerboard(Fdet2_mu, Fdet2_mu_oe);
    InsertForce(Fdet1,Fdet1_mu,mu);
    InsertForce(Fdet2,Fdet2_mu,mu);

    force= (-0.5)*( Fdet1 + Fdet2);
    RealD t1 = usecond();
    std::cout << GridLogMessage << " logDetJacobianForce t3-t0 "<<t3a-t0<<" us "<<std::endl;
    std::cout << GridLogMessage << " logDetJacobianForce t4-t3b dJdXe_nMpInv+even part of Fdet1_mu & Fdet2_mu "<<t4-t3b<<" us "<<std::endl;
    std::cout << GridLogMessage << " logDetJacobianForce t3b-t3a dJdXe_nMpInv  "<<t3b-t3a<<" us "<<std::endl;
    std::cout << GridLogMessage << " logDetJacobianForce t3c-t3b even part of Fdet1_mu & Fdet2_mu  "<<t3c-t3b<<" us "<<std::endl;
    std::cout << GridLogMessage << " logDetJacobianForce t5-t4 mu nu loop "<<t5-t4<<" us Plaq "
	      <<tLR/1e3<<" ms Nxy "<<tNxy/1e3<<" ms MpInvJx_dNxxdSy "<<tMJx/1e3<<" ms "<<" insert_force "<<t_ins/1e3<< " ms Stencil "
	      <<t_stencil/1e3<<" ms Cshift "<<t_cshift/1e3<<" ms"<<std::endl;
    std::cout << GridLogMessage << " logDetJacobianForce t1-t5 "<<t1-t5<<" us "<<std::endl; // turn adj vec to SU3 force
    std::cout << GridLogMessage << " logDetJacobianForce level took "<<t1-t0<<" us "<<std::endl;
  }

  //Assume: masking is based on red-black checkerboarding
  RealD logDetJacobianLevel(const GaugeField &U,int smr)
  {
    GRID_TRACE("logDetJacobianLevel");
    GridBase* grid = U.Grid();
    GridBase* hgrid = UrbGrid; // For now, assume masking is based on red-black checkerboarding
    assert(grid==UGrid);
    GaugeField gU(grid);
    std::vector<GaugeLinkField> gUmu(Nd,grid);
    GaugeLinkField PlaqL(hgrid);
    GaugeLinkField Z(hgrid);
    GaugeLinkField Umu(hgrid), Cmu(hgrid);
    LatticeComplex ln_det(hgrid); 
    typedef typename SU3Adjoint::LatticeAdjMatrix  AdjMatrixField;
    AdjMatrixField  Ncb(hgrid);
    AdjMatrixField  Zac(hgrid); 
    ColourMatrix Ident;
    
    RealD time=0, tN=0, tZ=0, t_M=0, tJ_lnDet=0, t_det=0, t_ln=0;
    int mu= (smr/2) %Nd;

    time -= usecond();
    Ident = ComplexD(1.0);

    int cb = cbs[smr]; //Assume: the applied mask is the cb mask

    Z.Checkerboard() = cb;
    Cmu.Checkerboard() = cb;
    PlaqL.Checkerboard() = cb;
    Ncb.Checkerboard() = cb;
    ln_det.Checkerboard() = cb; 

    pickCheckerboard(cb,Umu,peekLorentz(U,mu));

    //////////////////////////////////////////////////////////////////
    // Assemble the N matrix
    //////////////////////////////////////////////////////////////////

    tN -= usecond();
    {GRID_TRACE("ExchangePeriodic");
    gU = Ghost.ExchangePeriodic(U);
    }
    double rho=this->StoutSmearing->SmearRho[1];
    PlaqL = Ident;
    BaseSmear_ghost(Cmu,gU,mu,rho);
    ComputeNxy(PlaqL,Umu * adj(Cmu),Ncb);//we can expand function body and form a bigger accelerated loop
    tN += usecond();

    //////////////////////////////////////////////////////////////////
    // Assemble Luscher exp diff map J matrix 
    //////////////////////////////////////////////////////////////////
    tZ -= usecond();
    {GRID_TRACE("Z");
    // Ta so Z lives in Lie algabra
    Z  = Ta(Cmu * adj(Umu));
    // Move Z to the Adjoint Rep == make_adjoint_representation
    SU3Adjoint::make_adjoint_rep(Zac, Z);
    }
    tZ += usecond();

    tJ_lnDet -= usecond();
    t_M -= usecond();
    {GRID_TRACE("Mab");
    autoView(ln_det_v,ln_det,AcceleratorWrite);
    autoView(Zac_v,Zac,AcceleratorRead);
    autoView(Ncb_v,Ncb,AcceleratorRead);
    accelerator_for(ss,hgrid->oSites(),hgrid->Nsimd(),{
	typedef decltype(coalescedRead(Zac_v(0)))    adj_mat;

	adj_mat X, Jac, Mab_ss;
	RealD kpfac = 1;
	
	//////////////////////////////////////
	// J(x) = 1 + Sum_k=1..N (-Zac)^k/(k+1)!
	//////////////////////////////////////
	X=1.0; 
	Jac = X;
	for(int k=1;k<12;k++){
	  X=(-1.0)*X*Zac_v(ss);
	  kpfac = kpfac /((RealD) (k+1));
	  Jac = Jac + X * kpfac;
	}
	
	////////////////////////////
	// Mab
	////////////////////////////
	Mab_ss = Complex(1.0,0.0);
	Mab_ss = Mab_ss - Jac * Ncb_v(ss);


	////////////////////////////
	// ln det
	////////////////////////////
	auto detD = Determinant(Mab_ss);
	coalescedWrite(ln_det_v[ss],log(detD));
      });
    }
    t_M+=usecond();

    tJ_lnDet += usecond();
    Complex result = sum(ln_det);
    time += usecond();
    std::cout << GridLogMessage << " logDetJacobianLevel " << time/1e3 <<" ms N "<<tN/1e3<<" ms Z "<<tZ/1e3<<" ms J_lnDetM " <<tJ_lnDet/1e3
	      <<" ms Det " <<t_det/1e3<<" ms log "<<t_ln/1e3<<" ms M "<<t_M/1e3<<std::endl;

    return result.real();
  }

  /*------------------    old implementation for testing      #########################   -------------------------------- */
  void Compute_MpInvJx_dNxxdSy(int old, const GaugeLinkField &PlaqL,const GaugeLinkField &PlaqR, AdjMatrixField MpInvJx,AdjVectorField &Fdet2 )
  {
    GRID_TRACE("Compute_MpInvJx_dNxxdSy_old");
    int cb = PlaqL.Checkerboard();
    GaugeLinkField UtaU(PlaqL.Grid());         UtaU.Checkerboard() = cb;
    GaugeLinkField D(PlaqL.Grid());            D.Checkerboard() = cb;
    AdjMatrixField Dbc(PlaqL.Grid());          Dbc.Checkerboard() = cb;
    AdjMatrixField Dbc_opt(PlaqL.Grid());      Dbc_opt.Checkerboard() = cb;
    LatticeComplex tmp(PlaqL.Grid());          tmp.Checkerboard() = cb;
    const int Ngen = SU3Adjoint::Dimension;
    Complex ci(0,1);
    ColourMatrix   ta,tb,tc;
    RealD t=0, time;
    RealD tp=0, tpl=0;
    RealD tta=0;
    RealD tpk=0;
    t-=usecond();
    for(int a=0;a<Ngen;a++) {
      tta-=usecond();
      SU3::generator(a, ta);
      ta = 2.0 * ci * ta;
      // Qlat Tb = 2i Tb^Grid
      {
	GRID_TRACE("UtaU");
	UtaU= adj(PlaqL)*ta*PlaqR; // 6ms
      }
      tta+=usecond();
      ////////////////////////////////////////////
      // Could add this entire C-loop to a projection routine
      // for performance. Could also pick checkerboard on UtaU
      // and set checkerboard on result for 2x perf
      ////////////////////////////////////////////
      for(int c=0;c<Ngen;c++) {
	SU3::generator(c, tc);
	tc = 2.0*ci*tc;
	tp-=usecond(); 
	{
	  GRID_TRACE("D = Ta tc*UtaU");
	  D = Ta( tc *UtaU); // 2ms
	}
#if 1
	{
	  GRID_TRACE("LieAlgebraProject");
	  time=-usecond();
	  SU3::LieAlgebraProject(Dbc_opt,D,c); // 5.5ms
	  time+=usecond();
	  tpl+=time;
	  //std::cout << GridLogMessage << " LieAlgebraProject_in_Compute_MpInvJx " << a << " " << c << " "<<time/1e3<<" ms"<<std::endl;
	}
#else
	for(int b=0;b<Ngen;b++){
	  SU3::generator(b, tb);
	  tmp =-trace(ci*tb*D); 
	  PokeIndex<ColourIndex>(Dbc,tmp,b,c);  // Adjoint rep
	}
#endif
	tp+=usecond();
      }
      //Dump(Dbc_opt,"Dbc_opt");
      //Dump(Dbc,"Dbc");
      tpk-=usecond();
      {
	GRID_TRACE("traceMpDbc");
	tmp = trace(MpInvJx * Dbc_opt);
      }
#if 0
      // Tesing trace_prod
      if (a==0){
	GRID_TRACE("traceMpDbc2");
	LatticeComplex tmp2(PlaqL.Grid()); tmp2.Checkerboard() = cb;
	SU3::trace_product(tmp2,MpInvJx,Dbc_opt);
	std::cout << GridLogMessage <<  "DEBUG: Compute_MpInvJx_dNxxdSy trace " << norm2(tmp2-tmp)<<std::endl;
      }
	
#endif	
      {
	GRID_TRACE("pokeIndecx");
	PokeIndex<ColourIndex>(Fdet2,tmp,a);
      }
      tpk+=usecond();
    }
    t+=usecond();
    std::cout << GridLogPerformance << " Compute_MpInvJx_dNxxdSy_old " << t/1e3 << " ms  proj "<<tp/1e3<< " ms"
	      << " ta "<<tta/1e3<<" ms" << " poke "<<tpk/1e3<< " ms LieAlgebraProject "<<tpl/1e3<<" ms"<<std::endl;
  }
  
  void ComputeNxy(int old, const GaugeLinkField &PlaqL,const GaugeLinkField &PlaqR,AdjMatrixField &NxAd)
  {
    GRID_TRACE("ComputeNxy_old");
    GaugeLinkField Nx(PlaqL.Grid());
    const int Ngen = SU3Adjoint::Dimension;
    Complex ci(0,1);
    ColourMatrix   tb;
    ColourMatrix   tc;
    RealD t = 0, tta = 0, tp = 0, tgen = 0;

    t-=usecond();
    for(int b=0;b<Ngen;b++) {
      tgen-=usecond();
      SU3::generator(b, tb);
      tgen+=usecond();
      tb = 2.0 * ci * tb;
      tta-=usecond();
      {
        GRID_TRACE("UtaU");
	Nx = Ta( adj(PlaqL)*tb * PlaqR );
      }
      tta+=usecond();
      tp-=usecond();
#if 1
      {
	GRID_TRACE("LieAlgebraProject");
	SU3::LieAlgebraProject(NxAd,Nx,b);
      }
#else
      for(int c=0;c<Ngen;c++) {
	SU3::generator(c, tc);
	auto tmp =closure( -trace(ci*tc*Nx));
	{
	  GRID_TRACE("pokeIndecx");
	  PokeIndex<ColourIndex>(NxAd,tmp,c,b);
	}
      }
#endif
      tp+=usecond();
    }
    t+=usecond();
    std::cout << GridLogPerformance << " ComputeNxy_old " << t/1e3 << " ms  proj "<<tp/1e3<< " ms"
              << " ta "<<tta/1e3<<" ms tgen "<< tgen/1e3 << std::endl;
  }
  void logDetJacobianForceLevel(int old, const GaugeField &U, GaugeField &force ,int smr)
  {
    GRID_TRACE("logDetJacobianForceLevel_old");
    GridBase* grid = U.Grid();
    ColourMatrix   tb;
    ColourMatrix   tc;
    ColourMatrix   ta;
    GaugeField C(grid);
    GaugeField Umsk(grid);
    std::vector<GaugeLinkField> Umu(Nd,grid);
    GaugeLinkField Cmu(grid); // U and staple; C contains factor of epsilon
    GaugeLinkField Zx(grid);  // U times Staple, contains factor of epsilon
    GaugeLinkField Nxx(grid);  // Nxx fundamental space
    GaugeLinkField Utmp(grid);
    GaugeLinkField PlaqL(grid);
    GaugeLinkField PlaqR(grid);
    const int Ngen = SU3Adjoint::Dimension;
    AdjMatrix TRb;
    ColourMatrix Ident;
    LatticeComplex  cplx(grid);
    
    AdjVectorField  dJdXe_nMpInv(grid); 
    AdjVectorField  dJdXe_nMpInv_y(grid); 
    AdjMatrixField  MpAd(grid);    // Mprime luchang's notes
    AdjMatrixField  MpAdInv(grid); // Mprime inverse
    AdjMatrixField  NxxAd(grid);    // Nxx in adjoint space
    AdjMatrixField  JxAd(grid);     
    AdjMatrixField  ZxAd(grid);
    AdjMatrixField  mZxAd(grid);
    AdjMatrixField  X(grid);
    Complex ci(0,1);

    RealD t0 = usecond();
    Ident = ComplexD(1.0);
    for(int d=0;d<Nd;d++){
      Umu[d] = peekLorentz(U, d);
    }
    int mu= (smr/2) %Nd;

    ////////////////////////////////////////////////////////////////////////////////
    // Mask the gauge field
    ////////////////////////////////////////////////////////////////////////////////
    auto mask=PeekIndex<LorentzIndex>(masks[smr],mu); // the cb mask

    Umsk = U;
    ApplyMask(Umsk,smr);
    Utmp = peekLorentz(Umsk,mu);

    ////////////////////////////////////////////////////////////////////////////////
    // Retrieve the eps/rho parameter(s) -- could allow all different but not so far
    ////////////////////////////////////////////////////////////////////////////////
    double rho=this->StoutSmearing->SmearRho[1];
    int idx=0;
    for(int mu=0;mu<4;mu++){
    for(int nu=0;nu<4;nu++){
      if ( mu!=nu) assert(this->StoutSmearing->SmearRho[idx]==rho);
      else         assert(this->StoutSmearing->SmearRho[idx]==0.0);
      idx++;
    }}
    //////////////////////////////////////////////////////////////////
    // Assemble the N matrix
    //////////////////////////////////////////////////////////////////
    // Computes ALL the staples -- could compute one only and do it here
    RealD time;
    time=-usecond();
    BaseSmear(Cmu, U,mu,rho);

    //////////////////////////////////////////////////////////////////
    // Assemble Luscher exp diff map J matrix 
    //////////////////////////////////////////////////////////////////
    // Ta so Z lives in Lie algabra
    {GRID_TRACE("Zx");
    Zx  = Ta(Cmu * adj(Umu[mu]));
    time+=usecond();
    }
    std::cout << GridLogMessage << "Full: Z took "<<time<< " us"<<std::endl;
    {GRID_TRACE("ZxAd");
    time=-usecond();
    // Move Z to the Adjoint Rep == make_adjoint_representation
    ZxAd = Zero();
    for(int b=0;b<8;b++) {
      // Adj group sets traceless antihermitian T's -- Guido, really????
      SU3::generator(b, tb);         // Fund group sets traceless hermitian T's
      SU3Adjoint::generator(b,TRb);
      TRb=-TRb;
      cplx = 2.0*trace(ci*tb*Zx); // my convention 1/2 delta ba
      ZxAd = ZxAd + cplx * TRb; // is this right? YES - Guido used Anti herm Ta's and with bloody wrong sign.
    }
    time+=usecond();
    }
    std::cout << GridLogMessage << "Full: ZxAd took "<<time<< " us"<<std::endl;

    //////////////////////////////////////
    // J(x) = 1 + Sum_k=1..N (-Zac)^k/(k+1)!
    //////////////////////////////////////
    {GRID_TRACE("Jx");
    time=-usecond();
    X=1.0; 
    JxAd = X;
    mZxAd = (-1.0)*ZxAd; 
    RealD kpfac = 1;
    for(int k=1;k<12;k++){
      X=X*mZxAd;
      kpfac = kpfac /(k+1);
      JxAd = JxAd + X * kpfac;
    }
    time+=usecond();
    }
    std::cout << GridLogMessage << "Full: Jx took "<<time<< " us"<<std::endl;

    //////////////////////////////////////
    // dJ(x)/dxe
    //////////////////////////////////////
    time=-usecond();
#if 1
    std::vector<AdjMatrixField>  dJdX;    dJdX.resize(8,grid);
    std::vector<AdjMatrix> TRb_s; TRb_s.resize(8);
    AdjMatrixField tbXn(grid);
    AdjMatrixField sumXtbX(grid);
    AdjMatrixField t2(grid);
    AdjMatrixField dt2(grid);
    AdjMatrixField t3(grid);
    AdjMatrixField dt3(grid);
    AdjMatrixField aunit(grid);
    
    {GRID_TRACE("dJdX");
    for(int b=0;b<8;b++){
      SU3Adjoint::generator(b, TRb_s[b]);
      dJdX[b] = TRb_s[b];
    }
    aunit = ComplexD(1.0);
    // Could put into an accelerator_for
    X  = (-1.0)*ZxAd; 
    t2 = X;
    for (int j = 12; j > 1; --j) {
      t3  = t2*(1.0 / (j + 1))  + aunit;
      t2  = X * t3;
      for(int b=0;b<8;b++){
	dJdX[b]= TRb_s[b] * t3 + X * dJdX[b]*(1.0 / (j + 1));
      }
    }
    for(int b=0;b<8;b++){
      dJdX[b] = -dJdX[b];
    }
    }
#else
    std::vector<AdjMatrixField>  dJdX;    dJdX.resize(8,grid);
    AdjMatrixField tbXn(grid);
    AdjMatrixField sumXtbX(grid);
    AdjMatrixField t2(grid);
    AdjMatrixField dt2(grid);
    AdjMatrixField t3(grid);
    AdjMatrixField dt3(grid);
    AdjMatrixField aunit(grid);
    for(int b=0;b<8;b++){
      aunit = ComplexD(1.0);
      SU3Adjoint::generator(b, TRb); //dt2

      X  = (-1.0)*ZxAd; 
      t2 = X;
      dt2 = TRb;
      for (int j = 12; j > 1; --j) {
	t3  = t2*(1.0 / (j + 1))  + aunit;
	dt3 = dt2*(1.0 / (j + 1));
	t2 = X * t3;
	dt2 = TRb * t3 + X * dt3;
      }
      dJdX[b] = -dt2; 
    }
#endif  
    time+=usecond();
    std::cout << GridLogMessage << "Full: dJx took "<<time<< " us"<<std::endl;
    /////////////////////////////////////////////////////////////////
    // Mask Umu for this link
    /////////////////////////////////////////////////////////////////
    time=-usecond();
    PlaqL = Ident;
    PlaqR = Utmp*adj(Cmu);
    ComputeNxy(old,PlaqL,PlaqR,NxxAd);
    time+=usecond();
    std::cout << GridLogMessage << "Full: ComputeNxy took "<<time<< " us"<<std::endl;
    
    ////////////////////////////
    // Mab
    ////////////////////////////
    MpAd = Complex(1.0,0.0);
    MpAd = MpAd - JxAd * NxxAd;

    /////////////////////////
    // invert the 8x8
    /////////////////////////
    time=-usecond();
    MpAdInv = Inverse(MpAd);
    time+=usecond();
    std::cout << GridLogMessage << "Full: MpAdInv took "<<time<< " us"<<std::endl;
    
    RealD t3a = usecond();
    /////////////////////////////////////////////////////////////////
    // Nxx Mp^-1
    /////////////////////////////////////////////////////////////////
    AdjVectorField  FdetV(grid);
    AdjVectorField  Fdet1_nu(grid);
    AdjVectorField  Fdet2_nu(grid);
    AdjVectorField  Fdet2_mu(grid);
    AdjVectorField  Fdet1_mu(grid);

    AdjMatrixField nMpInv(grid);
    nMpInv= NxxAd *MpAdInv;

    AdjMatrixField MpInvJx(grid);
    AdjMatrixField MpInvJx_nu(grid);
    MpInvJx = (-1.0)*MpAdInv * JxAd;// rho is on the plaq factor

    Compute_MpInvJx_dNxxdSy(old,PlaqL,PlaqR,MpInvJx,FdetV);
    Fdet2_mu=FdetV;
    Fdet1_mu=Zero();
    
    for(int e =0 ; e<8 ; e++){
      LatticeComplexD tr(grid);
      //      ColourMatrix te;
      //      SU3::generator(e, te);
      tr = trace(dJdX[e] * nMpInv);
      pokeColour(dJdXe_nMpInv,tr,e);
    }
    ///////////////////////////////
    // Mask it off
    ///////////////////////////////
    auto tmp=PeekIndex<LorentzIndex>(masks[smr],mu);
    dJdXe_nMpInv = dJdXe_nMpInv*tmp;
    
    //    dJdXe_nMpInv needs to multiply:
    //       Nxx_mu (site local)                           (1)
    //       Nxy_mu one site forward  in each nu direction (3)
    //       Nxy_mu one site backward in each nu direction (3)
    //       Nxy_nu 0,0  ; +mu,0; 0,-nu; +mu-nu   [ 3x4 = 12]
    // 19 terms.
    AdjMatrixField Nxy(grid);

    GaugeField Fdet1(grid);
    GaugeField Fdet2(grid);
    GaugeLinkField Fdet_pol(grid); // one polarisation

    RealD t4 = usecond(), tLR = 0, tNxy = 0, tMJx = 0;
    for(int nu=0;nu<Nd;nu++){

      if (nu!=mu) {
	///////////////// +ve nu /////////////////
	//     __
	//    |  |
	//    x==    // nu polarisation -- clockwise

	time=-usecond(); tLR -= usecond();
	PlaqL=Ident;

	PlaqR=(-rho)*Gimpl::CovShiftForward(Umu[nu], nu,
 	       Gimpl::CovShiftForward(Umu[mu], mu,
	         Gimpl::CovShiftBackward(Umu[nu], nu,
		   Gimpl::CovShiftIdentityBackward(Utmp, mu))));
	time+=usecond(); tLR += usecond();
	std::cout << GridLogMessage << "Full: PlaqLR took "<<time<< " us"<<std::endl;

	time=-usecond(); tNxy -= usecond();
	dJdXe_nMpInv_y =   dJdXe_nMpInv;
	ComputeNxy(old,PlaqL,PlaqR,Nxy);
	Fdet1_nu = transpose(Nxy)*dJdXe_nMpInv_y;
	time+=usecond(); tNxy += usecond();
	std::cout << GridLogMessage << "Full: ComputeNxy (occurs 6x) took "<<time<< " us"<<std::endl;

	time=-usecond(); tMJx -= usecond();
	PlaqR=(-1.0)*PlaqR;
	Compute_MpInvJx_dNxxdSy(old,PlaqL,PlaqR,MpInvJx,FdetV);
	Fdet2_nu = FdetV;
	time+=usecond(); tMJx += usecond();
	std::cout << GridLogMessage << "Full: Compute_MpInvJx_dNxxSy (occurs 6x) took "<<time<< " us"<<std::endl;
	
	//    x==
	//    |  |
	//    .__|    // nu polarisation -- anticlockwise

	tLR -= usecond();
	PlaqR=(rho)*Gimpl::CovShiftForward(Umu[nu], nu,
		      Gimpl::CovShiftBackward(Umu[mu], mu,
    	 	        Gimpl::CovShiftIdentityBackward(Umu[nu], nu)));

	PlaqL=Gimpl::CovShiftIdentityBackward(Utmp, mu);
	tLR += usecond();

	tNxy -= usecond();
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,mu,-1);
	ComputeNxy(old,PlaqL, PlaqR,Nxy);
	Fdet1_nu = Fdet1_nu+transpose(Nxy)*dJdXe_nMpInv_y;
	tNxy += usecond();

	tMJx -= usecond();
	MpInvJx_nu = Cshift(MpInvJx,mu,-1);
	Compute_MpInvJx_dNxxdSy(old,PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_nu = Fdet2_nu+FdetV;
	tMJx += usecond();
	
	///////////////// -ve nu /////////////////
	//  __
	// |  |
	// x==          // nu polarisation -- clockwise

	tLR -= usecond();
	PlaqL=(rho)* Gimpl::CovShiftForward(Umu[mu], mu,
		       Gimpl::CovShiftForward(Umu[nu], nu,
			 Gimpl::CovShiftIdentityBackward(Utmp, mu)));

        PlaqR = Gimpl::CovShiftIdentityForward(Umu[nu], nu);
	tLR += usecond();

	tNxy -= usecond();
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,1);
	ComputeNxy(old,PlaqL,PlaqR,Nxy);
	Fdet1_nu = Fdet1_nu + transpose(Nxy)*dJdXe_nMpInv_y;
	tNxy += usecond();

	tMJx -= usecond();
	MpInvJx_nu = Cshift(MpInvJx,nu,1);
	Compute_MpInvJx_dNxxdSy(old,PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_nu = Fdet2_nu+FdetV;
	tMJx += usecond();
	
	// x==
	// |  |
	// |__|         // nu polarisation

	tLR -= usecond();
	PlaqL=(-rho)*Gimpl::CovShiftForward(Umu[nu], nu,
 	        Gimpl::CovShiftIdentityBackward(Utmp, mu));

	PlaqR=Gimpl::CovShiftBackward(Umu[mu], mu,
	        Gimpl::CovShiftIdentityForward(Umu[nu], nu));
	tLR += usecond();

	tNxy -= usecond();
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,mu,-1);
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv_y,nu,1);

	ComputeNxy(old,PlaqL,PlaqR,Nxy);
	Fdet1_nu = Fdet1_nu + transpose(Nxy)*dJdXe_nMpInv_y;
	tNxy += usecond();

	tMJx -= usecond();
	MpInvJx_nu = Cshift(MpInvJx,mu,-1);
	MpInvJx_nu = Cshift(MpInvJx_nu,nu,1);
	Compute_MpInvJx_dNxxdSy(old,PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_nu = Fdet2_nu+FdetV;
	tMJx += usecond();
	
	/////////////////////////////////////////////////////////////////////
	// Set up the determinant force contribution in 3x3 algebra basis
	/////////////////////////////////////////////////////////////////////
	InsertForce(Fdet1,Fdet1_nu,nu);
	InsertForce(Fdet2,Fdet2_nu,nu);
	
	//////////////////////////////////////////////////
	// Parallel direction terms
	//////////////////////////////////////////////////

        //     __
	//    |  "
	//    |__"x    // mu polarisation
	tLR -= usecond();
	PlaqL=(-rho)*Gimpl::CovShiftForward(Umu[mu], mu,
		      Gimpl::CovShiftBackward(Umu[nu], nu,
   		        Gimpl::CovShiftIdentityBackward(Utmp, mu)));

	PlaqR=Gimpl::CovShiftIdentityBackward(Umu[nu], nu);
	tLR += usecond();

	tNxy -= usecond();
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,-1);

	ComputeNxy(old,PlaqL,PlaqR,Nxy);
	Fdet1_mu = Fdet1_mu + transpose(Nxy)*dJdXe_nMpInv_y;
	tNxy += usecond();

	tMJx -= usecond();
	MpInvJx_nu = Cshift(MpInvJx,nu,-1);

	Compute_MpInvJx_dNxxdSy(old,PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_mu = Fdet2_mu+FdetV;
	tMJx += usecond();

	//  __
	// "  |
	// x__|          // mu polarisation
	tLR -= usecond();
	PlaqL=(-rho)*Gimpl::CovShiftForward(Umu[mu], mu,
		       Gimpl::CovShiftForward(Umu[nu], nu,
		 	 Gimpl::CovShiftIdentityBackward(Utmp, mu)));

        PlaqR=Gimpl::CovShiftIdentityForward(Umu[nu], nu);
	tLR += usecond();

	tNxy -= usecond();
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,1);

	ComputeNxy(old,PlaqL,PlaqR,Nxy);
	Fdet1_mu = Fdet1_mu + transpose(Nxy)*dJdXe_nMpInv_y;
	tNxy += usecond();

	tMJx -= usecond();
	MpInvJx_nu = Cshift(MpInvJx,nu,1);

	Compute_MpInvJx_dNxxdSy(old,PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_mu = Fdet2_mu+FdetV;
	tMJx += usecond();
      }
    }
    RealD t5 = usecond();

    Fdet1_mu = Fdet1_mu + transpose(NxxAd)*dJdXe_nMpInv;

    InsertForce(Fdet1,Fdet1_mu,mu);
    InsertForce(Fdet2,Fdet2_mu,mu);

    force= (-0.5)*( Fdet1 + Fdet2);
    RealD t1 = usecond();
    std::cout << GridLogMessage << " Full: logDetJacobianForce t3-t0 "<<t3a-t0<<" us "<<std::endl;
    std::cout << GridLogMessage << " Full: logDetJacobianForce t4-t3 dJdXe_nMpInv "<<t4-t3a<<" us "<<std::endl;
    std::cout << GridLogMessage << " Full: logDetJacobianForce t5-t4 mu nu loop "<<t5-t4<<" us Plaq "
	      <<tLR/1e3<<" ms Nxy "<<tNxy/1e3<<" ms MpInvJx_dNxxdSy "<<tMJx/1e3<<" ms"<<std::endl;
    std::cout << GridLogMessage << " Full: logDetJacobianForce t1-t5 "<<t1-t5<<" us "<<std::endl; // turn adj vec to SU3 force
    std::cout << GridLogMessage << " Full: logDetJacobianForce level took "<<t1-t0<<" us "<<std::endl;
  }
  RealD logDetJacobianLevel(int old, const GaugeField &U,int smr)
  {
    GridBase* grid = U.Grid();
    GaugeField C(grid);
    GaugeLinkField Nb(grid);
    GaugeLinkField Z(grid);
    GaugeLinkField Umu(grid), Cmu(grid);
    ColourMatrix   Tb;
    ColourMatrix   Tc;
    typedef typename SU3Adjoint::AMatrix AdjMatrix;
    typedef typename SU3Adjoint::LatticeAdjMatrix  AdjMatrixField;
    typedef typename SU3Adjoint::LatticeAdjVector  AdjVectorField;
    const int Ngen = SU3Adjoint::Dimension;
    AdjMatrix TRb;
    LatticeComplex       cplx(grid); 
    AdjVectorField  AlgV(grid); 
    AdjMatrixField  Mab(grid);
    AdjMatrixField  Ncb(grid), Ncb_opt(grid);
    AdjMatrixField  Jac(grid);
    AdjMatrixField  Zac(grid);
    AdjMatrixField  mZac(grid);
    AdjMatrixField  X(grid);

    RealD time=0, tta=0, tpk=0, tN=0, tZ=0, tJ=0, tlnDetM=0;
    int mu= (smr/2) %Nd;

    time -= usecond();
    auto mask=PeekIndex<LorentzIndex>(masks[smr],mu); // the cb mask

    //////////////////////////////////////////////////////////////////
    // Assemble the N matrix
    //////////////////////////////////////////////////////////////////

    tN -= usecond();
    double rho=this->StoutSmearing->SmearRho[1];
    BaseSmear(Cmu, U,mu,rho);

    Umu = peekLorentz(U, mu);
    Complex ci(0,1);
#if 1
    for(int b=0;b<Ngen;b++) {
      SU3::generator(b, Tb);
      // Qlat Tb = 2i Tb^Grid
      tta -= usecond();
      Nb = (2.0)*Ta( ci*Tb * Umu * adj(Cmu));
      tta += usecond();
      // FIXME -- replace this with LieAlgebraProject
#if 0
      // Fixed it
      SU3::LieAlgebraProject(Ncb_opt,Nb,b);
#else
      for(int c=0;c<Ngen;c++) {
	SU3::generator(c, Tc);
	auto tmp = -trace(ci*Tc*Nb); // Luchang's norm: (2Tc) (2Td) N^db = -2 delta cd N^db // - was important
	tpk -= usecond();
	PokeIndex<ColourIndex>(Ncb,tmp,c,b);
	tpk += usecond();
      }

#endif
    }
#else
      autoView(NxAd_v,NxAd,AcceleratorWrite);
    autoView(PlaqL_v,PlaqL,AcceleratorRead);
    autoView(PlaqR_v,PlaqR,AcceleratorRead);
    const int nsimd = vAlgebraMatrix::Nsimd();
    accelerator_for(ss,grid->oSites(),nsimd,{
        typedef decltype(coalescedRead(PlaqL_v[0])) SU3_mat;
        typedef decltype(coalescedRead(NxAd_v[0]))  adj_mat;
        SU3_mat Nx;
        adj_mat NxAd_site;
        for(int b=0;b<Ngen;b++) {
          SU3::generator(b, tb);
          tb = 2.0 * ci * tb;
          //auto Nx =Ta( adj(PlaqL_v(ss)) * tb * PlaqR_v(ss) );
          auto Nx = 0.5*( adj(PlaqL_v(ss))*tb*PlaqR_v(ss) + adj(PlaqR_v(ss))*tb*PlaqL_v(ss) );//-tr(*) part does not contribute
          SU3::LieAlgebraProject(NxAd_site,Nx,b);
        }
        coalescedWrite(NxAd_v[ss],NxAd_site);
      });
#endif
    //Dump(Ncb_opt,"Ncb_opt");
    //Dump(Ncb,"Ncb");
    tN += usecond();
    
    //////////////////////////////////////////////////////////////////
    // Assemble Luscher exp diff map J matrix 
    //////////////////////////////////////////////////////////////////
    tZ -= usecond();
    // Ta so Z lives in Lie algabra
    Z  = Ta(Cmu * adj(Umu));

    // Move Z to the Adjoint Rep == make_adjoint_representation
    Zac = Zero();
    for(int b=0;b<8;b++) {
      // Adj group sets traceless antihermitian T's -- Guido, really????
      // Is the mapping of these the same? Same structure constants
      // Might never have been checked.
      SU3::generator(b, Tb);         // Fund group sets traceless hermitian T's
      SU3Adjoint::generator(b,TRb);
      TRb=-TRb;
      cplx = 2.0*trace(ci*Tb*Z); // my convention 1/2 delta ba
      Zac = Zac + cplx * TRb; // is this right? YES - Guido used Anti herm Ta's and with bloody wrong sign.
    }
    tZ += usecond();
    
    //////////////////////////////////////
    // J(x) = 1 + Sum_k=1..N (-Zac)^k/(k+1)!
    //////////////////////////////////////
    tJ -= usecond();
    X=1.0; 
    Jac = X;
    mZac = (-1.0)*Zac; 
    RealD kpfac = 1;
    for(int k=1;k<12;k++){
      X=X*mZac;
      kpfac = kpfac /(k+1);
      Jac = Jac + X * kpfac;
    }
    tJ += usecond();

    tlnDetM -= usecond();
    ////////////////////////////
    // Mab
    ////////////////////////////
    Mab = Complex(1.0,0.0);
    Mab = Mab - Jac * Ncb;

    ////////////////////////////
    // det
    ////////////////////////////
    LatticeComplex       det(grid); 
    det = Determinant(Mab);

    ////////////////////////////
    // ln det
    ////////////////////////////
    LatticeComplex       ln_det(grid); 
    ln_det = log(det);

    ////////////////////////////
    // Masked sum
    ////////////////////////////
    ln_det = ln_det * mask;
    tlnDetM += usecond();
    time += usecond();
    std::cout << GridLogMessage << " Full: logDetJacobianLevel " << time/1e3 << " ms ta "<<tta/1e3<<" ms" << " poke "<<tpk/1e3<< " ms"
	      <<" N "<<tN/1e3<<" ms Z "<<tZ/1e3<<" ms J " <<tJ/1e3<<" ms lnDetM "<<tlnDetM/1e3<<" ms" <<std::endl;
    
    Complex result = sum(ln_det);
    return result.real();
  }
  
public:
  RealD logDetJacobian(int old)
  {
    RealD ln_det = 0;
    if (this->smearingLevels > 0)
    {
      double start = usecond();
      for (int ismr = this->smearingLevels - 1; ismr > 0; --ismr) {
	ln_det+= logDetJacobianLevel(old,this->get_smeared_conf(ismr-1),ismr);
      }
      ln_det +=logDetJacobianLevel(old,*(this->ThinLinks),0);

      double end = usecond();
      double time = (end - start)/ 1e3;
      std::cout << GridLogMessage << "Full: GaugeConfigurationMasked: logDetJacobian took " << time << " ms" << std::endl;  
    }
    std::cout << GridLogMessage << " DEBUG: logDetJacobian Full " << std::endl;
    return ln_det;
  }
  void logDetJacobianForce(int old, GaugeField &force)
  {
    GRID_TRACE("logDetJacobianForce");
    force =Zero();
    GaugeField force_det(force.Grid());

    if (this->smearingLevels > 0)
    {
      double start = usecond(), tas = 0;

      GaugeLinkField tmp_mu(force.Grid());

      for (int ismr = this->smearingLevels - 1; ismr > 0; --ismr) {

	// remove U in UdSdU...
	for (int mu = 0; mu < Nd; mu++) {
	  tmp_mu = adj(peekLorentz(this->get_smeared_conf(ismr), mu)) * peekLorentz(force, mu);
	  pokeLorentz(force, tmp_mu, mu);
	}

	tas -= usecond();
      	// Propagate existing force
        force = this->AnalyticSmearedForce(force, this->get_smeared_conf(ismr - 1), ismr);
	tas += usecond();
	
	// Add back U in UdSdU...
	for (int mu = 0; mu < Nd; mu++) {
	  tmp_mu = peekLorentz(this->get_smeared_conf(ismr - 1), mu) * peekLorentz(force, mu);
	  pokeLorentz(force, tmp_mu, mu);
	}
    	
	// Get this levels determinant force
	force_det = Zero();
	logDetJacobianForceLevel(old, this->get_smeared_conf(ismr-1),force_det,ismr);

	// Sum the contributions
	force = force + force_det;
      }
    
      // remove U in UdSdU...
      for (int mu = 0; mu < Nd; mu++) {
	tmp_mu = adj(peekLorentz(this->get_smeared_conf(0), mu)) * peekLorentz(force, mu);
	pokeLorentz(force, tmp_mu, mu);
      }

      tas -= usecond();
      force = this->AnalyticSmearedForce(force, *this->ThinLinks,0);
      tas += usecond();
      
      for (int mu = 0; mu < Nd; mu++) {
	tmp_mu = peekLorentz(*this->ThinLinks, mu) * peekLorentz(force, mu);
	pokeLorentz(force, tmp_mu, mu);
      }

      force_det = Zero();

      logDetJacobianForceLevel(old, *this->ThinLinks,force_det,0);

      force = force + force_det;

      force=Ta(force); // Ta
      
    }  // if smearingLevels = 0 do nothing
    std::cout << GridLogMessage << " DEBUG: logDetJacobianForce Full " << std::endl;
  }
  /*------------------------------- OLD IMPLEMENTATION ----------------------------------------------------------*/
  
  RealD logDetJacobian(void)
  {
    RealD ln_det = 0;
    if (this->smearingLevels > 0)
    {
      double start = usecond();
      for (int ismr = this->smearingLevels - 1; ismr > 0; --ismr) {
	ln_det+= logDetJacobianLevel(this->get_smeared_conf(ismr-1),ismr);
      }
      ln_det +=logDetJacobianLevel(*(this->ThinLinks),0);

      double end = usecond();
      double time = (end - start)/ 1e3;
      std::cout << GridLogMessage << "GaugeConfigurationMasked: logDetJacobian took " << time << " ms" << std::endl;
#if 0 //DEBUG
      RealD ln_det2 = logDetJacobian(1);
      std::cout << GridLogMessage << " DEBUG: logDetJacobian_diff " << abs(ln_det2-ln_det) <<":"<<ln_det<<" "<<ln_det2<<std::endl;
#endif
    }
    return ln_det;
  }
  void logDetJacobianForce(GaugeField &force)
  {
    GRID_TRACE("logDetJacobianForce");
    force =Zero();
    GaugeField force_det(force.Grid());

    if (this->smearingLevels > 0)
    {
      double start = usecond(), tas = 0;

      GaugeLinkField tmp_mu(force.Grid());

      for (int ismr = this->smearingLevels - 1; ismr > 0; --ismr) {

	// remove U in UdSdU...
	for (int mu = 0; mu < Nd; mu++) {
	  tmp_mu = adj(peekLorentz(this->get_smeared_conf(ismr), mu)) * peekLorentz(force, mu);
	  pokeLorentz(force, tmp_mu, mu);
	}

	tas -= usecond();
      	// Propagate existing force
        force = this->AnalyticSmearedForce(force, this->get_smeared_conf(ismr - 1), ismr);
	tas += usecond();
	
	// Add back U in UdSdU...
	for (int mu = 0; mu < Nd; mu++) {
	  tmp_mu = peekLorentz(this->get_smeared_conf(ismr - 1), mu) * peekLorentz(force, mu);
	  pokeLorentz(force, tmp_mu, mu);
	}
    	
	// Get this levels determinant force
	force_det = Zero();
	logDetJacobianForceLevel(this->get_smeared_conf(ismr-1),force_det,ismr);

	// Sum the contributions
	force = force + force_det;
      }
    
      // remove U in UdSdU...
      for (int mu = 0; mu < Nd; mu++) {
	tmp_mu = adj(peekLorentz(this->get_smeared_conf(0), mu)) * peekLorentz(force, mu);
	pokeLorentz(force, tmp_mu, mu);
      }

      tas -= usecond();
      force = this->AnalyticSmearedForce(force, *this->ThinLinks,0);
      tas += usecond();
      
      for (int mu = 0; mu < Nd; mu++) {
	tmp_mu = peekLorentz(*this->ThinLinks, mu) * peekLorentz(force, mu);
	pokeLorentz(force, tmp_mu, mu);
      }

      force_det = Zero();

      logDetJacobianForceLevel(*this->ThinLinks,force_det,0);

      force = force + force_det;

      force=Ta(force); // Ta
      
#if 0 // debug
      GaugeField force_debug(force.Grid()); 
      logDetJacobianForce(1,force_debug);
      std::cout << GridLogMessage << " DEBUG: logDetJacobianForce_diff " << norm2(force-force_debug) <<":"<<norm2(force)<<" "<<norm2(force_debug)<< std::endl;
#endif
      double end = usecond();
      double time = (end - start)/ 1e3;
      std::cout << GridLogMessage << "GaugeConfigurationMasked: AnalyticSmearedForce in lnDetJacobianForce took " << tas/1e3 << " ms" << std::endl;
      std::cout << GridLogMessage << "GaugeConfigurationMasked: lnDetJacobianForce took " << time << " ms" << std::endl;  
    }  // if smearingLevels = 0 do nothing
  }

private:
  //====================================================================
  // Override base clas here to mask it
  virtual void fill_smearedSet(GaugeField &U)
  {
    this->ThinLinks = &U;  // attach the smearing routine to the field U

    // check the pointer is not null
    if (this->ThinLinks == NULL)
      std::cout << GridLogError << "[SmearedConfigurationMasked] Error in ThinLinks pointer\n";

    if (this->smearingLevels > 0)
    {
      std::cout << GridLogMessage << "[SmearedConfigurationMasked] Filling SmearedSet\n";
      GaugeField previous_u(this->ThinLinks->Grid());

      GaugeField smeared_A(this->ThinLinks->Grid());
      GaugeField smeared_B(this->ThinLinks->Grid());

      previous_u = *this->ThinLinks;
      double start = usecond();
      for (int smearLvl = 0; smearLvl < this->smearingLevels; ++smearLvl)
      {
        this->StoutSmearing->smear(smeared_A, previous_u);
	ApplyMask(smeared_A,smearLvl);
	smeared_B = previous_u;
	ApplyMask(smeared_B,smearLvl);
	// Replace only the masked portion
	this->SmearedSet[smearLvl] = previous_u-smeared_B + smeared_A;
        previous_u = this->SmearedSet[smearLvl];

        // For debug purposes
        RealD impl_plaq = WilsonLoops<Gimpl>::avgPlaquette(previous_u);
        std::cout << GridLogMessage << "[SmearedConfigurationMasked] smeared Plaq: " << impl_plaq << std::endl;
      }
      double end = usecond();
      double time = (end - start)/ 1e3;
      std::cout << GridLogMessage << "GaugeConfigurationMasked: Link smearing took " << time << " ms" << std::endl;  
    }
  }
  //====================================================================
  // Override base to add masking
  virtual GaugeField AnalyticSmearedForce(const GaugeField& SigmaKPrime,
					  const GaugeField& GaugeK,int level) 
  {
    GRID_TRACE("AnalyticSmearedForce");
    GridBase* grid = GaugeK.Grid();
    GaugeField SigmaK(grid), iLambda(grid);
    GaugeField SigmaKPrimeA(grid);
    GaugeField SigmaKPrimeB(grid);
    GaugeLinkField iLambda_mu(grid);
    GaugeLinkField iQ(grid), e_iQ(grid);
    GaugeLinkField SigmaKPrime_mu(grid);
    GaugeLinkField GaugeKmu(grid), Cmu(grid);

    int mmu= (level/2) %Nd;
    int cb= (level%2);
    double rho=this->StoutSmearing->SmearRho[1];

    RealD time = 0;

    time -= usecond();
    // Can override this to do one direction only.
    SigmaK = Zero();
    iLambda = Zero();

    SigmaKPrimeA = SigmaKPrime;
    ApplyMask(SigmaKPrimeA,level);
    SigmaKPrimeB = SigmaKPrime - SigmaKPrimeA;
    
    // Could get away with computing only one polarisation here
    // int mu= (smr/2) %Nd;
    // SigmaKprime_A has only one component
#if 0
    BaseSmear(Cmu, GaugeK,mu,rho);
    GaugeKmu = peekLorentz(GaugeK, mu);
    SigmaKPrime_mu = peekLorentz(SigmaKPrimeA, mu);
    iQ = Ta(Cmu * adj(GaugeKmu));
    this->set_iLambda(iLambda_mu, e_iQ, iQ, SigmaKPrime_mu, GaugeKmu);
    pokeLorentz(SigmaK, SigmaKPrime_mu * e_iQ + adj(Cmu) * iLambda_mu, mu);
    pokeLorentz(iLambda, iLambda_mu, mu);
    BaseSmearDerivative(SigmaK, iLambda,GaugeK,mu,rho);  // derivative of SmearBase
#else
    //    GaugeField C(grid);
    //    this->StoutSmearing->BaseSmear(C, GaugeK);
    //    for (int mu = 0; mu < Nd; mu++)
    int mu =mmu;
    BaseSmear(Cmu, GaugeK,mu,rho);
    {
      // Cmu = peekLorentz(C, mu);
      GaugeKmu = peekLorentz(GaugeK, mu);
      SigmaKPrime_mu = peekLorentz(SigmaKPrimeA, mu);
      iQ = Ta(Cmu * adj(GaugeKmu));
      this->set_iLambda(iLambda_mu, e_iQ, iQ, SigmaKPrime_mu, GaugeKmu);
      pokeLorentz(SigmaK, SigmaKPrime_mu * e_iQ + adj(Cmu) * iLambda_mu, mu);
      pokeLorentz(iLambda, iLambda_mu, mu);
      std::cout << " mu "<<mu<<" SigmaKPrime_mu "<<norm2(SigmaKPrime_mu)<< " iLambda_mu " <<norm2(iLambda_mu)<<std::endl;
    }
    //    GaugeField SigmaKcopy(grid);
    //    SigmaKcopy = SigmaK;
    BaseSmearDerivative(SigmaK, iLambda,GaugeK,mu,rho);  // derivative of SmearBase
    //    this->StoutSmearing->derivative(SigmaK, iLambda,GaugeK);  // derivative of SmearBase
    //    SigmaKcopy = SigmaKcopy - SigmaK;
    //    std::cout << " BaseSmearDerivative fast path error" <<norm2(SigmaKcopy)<<std::endl;
#endif
    ////////////////////////////////////////////////////////////////////////////////////
    // propagate the rest of the force as identity map, just add back
    ////////////////////////////////////////////////////////////////////////////////////
    SigmaK = SigmaK+SigmaKPrimeB;
    time += usecond();

    std::cout << GridLogMessage << "GaugeConfigurationMasked: Analytic smearing force  took " << time/1e3 << " ms" << std::endl;
    return SigmaK;
  }

public:

  /* Standard constructor */
  virtual ~SmearedConfigurationMasked()
  {
    delete UrbGrid;
    gStencils.clear();
    gStencils_smear.clear();
  }
  SmearedConfigurationMasked(GridCartesian* _UGrid, unsigned int Nsmear, Smear_Stout<Gimpl>& Stout)
    : SmearedConfiguration<Gimpl>(_UGrid, Nsmear,Stout),
      UGrid(_UGrid),
      Ghost(1,_UGrid)
  {
    assert(Nsmear%(2*Nd)==0); // Or multiply by 8??

    // was resized in base class
    assert(this->SmearedSet.size()==Nsmear);
    
    UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(_UGrid);
    
    LatticeComplex one(_UGrid); one = ComplexD(1.0,0.0);
    LatticeComplex tmp(_UGrid);
    for (unsigned int i = 0; i < this->smearingLevels; ++i) {

      masks.push_back(*(new LatticeLorentzComplex(_UGrid)));

      int mu= (i/2) %Nd;
      int cb= (i%2);
      LatticeComplex tmpcb(UrbGrid);

      cbs.push_back(cb);
      
      masks[i]=Zero();
      ////////////////////
      // Setup the mask
      ////////////////////
      tmp = Zero();
      pickCheckerboard(cb,tmpcb,one);
      setCheckerboard(tmp,tmpcb);
      PokeIndex<LorentzIndex>(masks[i],tmp, mu);
	
    }
    ///////////////////////////////
    // Setup stencils for staples
    ///////////////////////////////
    int depth            = 1;
    Coordinate local     =UGrid->LocalDimensions();
    Coordinate simd      =UGrid->_simd_layout;
    Coordinate processors=UGrid->_processors;
    Coordinate plocal    =UGrid->LocalDimensions();
    Coordinate global(Nd);

    // pre-stencil calculation for force calculation
    GridBase *ggrid = Ghost.grids[Nd-1];
    std::vector<Coordinate> shifts;

    gStencils.clear();
    Coordinate shift_0(Nd,0);
    for(int mu=0;mu<Nd;mu++){
      Coordinate shift_mu(Nd,0);  shift_mu[mu]=1;
      Coordinate shift_mmu(Nd,0); shift_mmu[mu]=-1;
      for(int nu=0;nu<Nd;nu++){
	if (nu!=mu) {
	  Coordinate shift_nu(Nd,0);  shift_nu[nu]=1;
	  Coordinate shift_mnu(Nd,0); shift_mnu[nu]=-1;
	  Coordinate shift_pmu_pnu(Nd,0); shift_pmu_pnu[mu]= 1;  shift_pmu_pnu[nu]= 1;
	  Coordinate shift_pmu_mnu(Nd,0); shift_pmu_mnu[mu]= 1;  shift_pmu_mnu[nu]=-1;
	  Coordinate shift_mmu_pnu(Nd,0); shift_mmu_pnu[mu]=-1;  shift_mmu_pnu[nu]= 1;
	  shifts.clear();
	  shifts.push_back(shift_0);shifts.push_back(shift_nu);shifts.push_back(shift_mu);shifts.push_back(shift_0);
	  gStencils.push_back(GeneralLocalStencil(ggrid,shifts));

	  shifts.clear();
	  shifts.push_back(shift_0);shifts.push_back(shift_mmu_pnu);shifts.push_back(shift_mmu);shifts.push_back(shift_mmu);
	  gStencils.push_back(GeneralLocalStencil(ggrid,shifts));

	  shifts.clear();
	  shifts.push_back(shift_0);shifts.push_back(shift_mu);shifts.push_back(shift_nu);shifts.push_back(shift_0);
	  gStencils.push_back(GeneralLocalStencil(ggrid,shifts));
	  
	  shifts.clear();
	  shifts.push_back(shift_0);shifts.push_back(shift_mmu_pnu);shifts.push_back(shift_mmu);shifts.push_back(shift_mmu);
	  gStencils.push_back(GeneralLocalStencil(ggrid,shifts));
	  
	  shifts.clear();
	  shifts.push_back(shift_0);shifts.push_back(shift_pmu_mnu);shifts.push_back(shift_mnu);shifts.push_back(shift_mnu);
	  gStencils.push_back(GeneralLocalStencil(ggrid,shifts));
	  
	  shifts.clear();
	  shifts.push_back(shift_0);shifts.push_back(shift_mu);shifts.push_back(shift_nu);shifts.push_back(shift_0);
	  gStencils.push_back(GeneralLocalStencil(ggrid,shifts));
	}
      }
    }
    // pre-stencil calculation for BaseSmear
    gStencils_smear.clear();
    for(int mu=0;mu<Nd;mu++){
      Coordinate shift_mu(Nd,0);  shift_mu[mu]=1;
      shifts.clear();
      for(int nu=0;nu<Nd;nu++){
        if (nu!=mu) {
	  Coordinate shift_nu(Nd,0);  shift_nu[nu]=1;
          Coordinate shift_mnu(Nd,0); shift_mnu[nu]=-1;
	  Coordinate shift_pmu_mnu(Nd,0); shift_pmu_mnu[mu]=1; shift_pmu_mnu[nu]=-1;

	  //U_nu(x+mu)U^dag_mu(x+nu) U^dag_nu(x)
	  shifts.push_back(shift_0); shifts.push_back(shift_nu); shifts.push_back(shift_mu);
	  //U_nu^dag(x-nu+mu) U_mu^dag(x-nu) U_nu(x-nu)
	  shifts.push_back(shift_mnu); shifts.push_back(shift_mnu); shifts.push_back(shift_pmu_mnu);
	}
      }
      gStencils_smear.push_back(GeneralLocalStencil(ggrid,shifts));
    }
  }
  
  virtual void smeared_force(GaugeField &SigmaTilde) 
  {
    GRID_TRACE("smeared_force");
    if (this->smearingLevels > 0)
    {
      double start = usecond();
      GaugeField force = SigmaTilde; // actually = U*SigmaTilde
      GaugeLinkField tmp_mu(SigmaTilde.Grid());

      // Remove U from UdSdU
      for (int mu = 0; mu < Nd; mu++)
      {
        // to get just SigmaTilde
        tmp_mu = adj(peekLorentz(this->SmearedSet[this->smearingLevels - 1], mu)) * peekLorentz(force, mu);
        pokeLorentz(force, tmp_mu, mu);
      }

      for (int ismr = this->smearingLevels - 1; ismr > 0; --ismr) {
        force = this->AnalyticSmearedForce(force, this->get_smeared_conf(ismr - 1),ismr);
      }
      
      force = this->AnalyticSmearedForce(force, *this->ThinLinks,0);

      // Add U to UdSdU
      for (int mu = 0; mu < Nd; mu++)
      {
        tmp_mu = peekLorentz(*this->ThinLinks, mu) * peekLorentz(force, mu);
        pokeLorentz(SigmaTilde, tmp_mu, mu);
      }


      double end = usecond();
      double time = (end - start)/ 1e3;
      std::cout << GridLogMessage << " GaugeConfigurationMasked: Smeared Force chain rule took " << time << " ms" << std::endl;

    }  // if smearingLevels = 0 do nothing
    SigmaTilde=Gimpl::projectForce(SigmaTilde); // Ta
  }

};

NAMESPACE_END(Grid);

