/*!
  @file GaugeConfiguration.h
  @brief Declares the GaugeConfiguration class
*/
#pragma once

NAMESPACE_BEGIN(Grid);

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
  
  std::vector<LatticeLorentzComplex> masks;

  typedef typename SU3Adjoint::AMatrix AdjMatrix;
  typedef typename SU3Adjoint::LatticeAdjMatrix  AdjMatrixField;
  typedef typename SU3Adjoint::LatticeAdjVector  AdjVectorField;

  // Adjoint vector to GaugeField force
  void InsertForce(GaugeField &Fdet,AdjVectorField &Fdet_nu,int nu)
  {
    Complex ci(0,1);
    GaugeLinkField Fdet_pol(Fdet.Grid());
    Fdet_pol=Zero();
    for(int e=0;e<8;e++){
      ColourMatrix te;
      SU3::generator(e, te);
      auto tmp=peekColour(Fdet_nu,e);
      Fdet_pol=Fdet_pol + ci*tmp*te; // but norm of te is different.. why?
    }
    pokeLorentz(Fdet, Fdet_pol, nu);
  }
  void Compute_MpInvJx_dNxxdSy(const GaugeLinkField &PlaqL,const GaugeLinkField &PlaqR, AdjMatrixField MpInvJx,AdjVectorField &Fdet2 )
  {
    GaugeLinkField UtaU(PlaqL.Grid());
    GaugeLinkField D(PlaqL.Grid());
    AdjMatrixField Dbc(PlaqL.Grid());
    LatticeComplex tmp(PlaqL.Grid());
    const int Ngen = SU3Adjoint::Dimension;
    Complex ci(0,1);
    ColourMatrix   ta,tb,tc;
    
    for(int a=0;a<Ngen;a++) {
      SU3::generator(a, ta);
      // Qlat Tb = 2i Tb^Grid
      UtaU= 2.0*ci*adj(PlaqL)*ta*PlaqR;
      for(int c=0;c<Ngen;c++) {
	SU3::generator(c, tc);
	D = Ta( (2.0)*ci*tc *UtaU);
	for(int b=0;b<Ngen;b++){
	  SU3::generator(b, tb);
	  tmp =-trace(ci*tb*D); 
	  PokeIndex<ColourIndex>(Dbc,tmp,b,c);  // Adjoint rep
	}
      }
      tmp = trace(MpInvJx * Dbc);
      PokeIndex<ColourIndex>(Fdet2,tmp,a);
    }
  }
  
  void ComputeNxy(const GaugeLinkField &PlaqL,const GaugeLinkField &PlaqR,AdjMatrixField &NxAd)
  {
    GaugeLinkField Nx(PlaqL.Grid());
    const int Ngen = SU3Adjoint::Dimension;
    Complex ci(0,1);
    ColourMatrix   tb;
    ColourMatrix   tc;
    for(int b=0;b<Ngen;b++) {
      SU3::generator(b, tb);
      Nx = (2.0)*Ta( adj(PlaqL)*ci*tb * PlaqR );
      for(int c=0;c<Ngen;c++) {
	SU3::generator(c, tc);
	auto tmp =closure( -trace(ci*tc*Nx)); 
	PokeIndex<ColourIndex>(NxAd,tmp,c,b); 
      }
    }
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
    this->StoutSmearing->BaseSmear(C, U);
    Cmu = peekLorentz(C, mu);

    //////////////////////////////////////////////////////////////////
    // Assemble Luscher exp diff map J matrix 
    //////////////////////////////////////////////////////////////////
    // Ta so Z lives in Lie algabra
    Zx  = Ta(Cmu * adj(Umu[mu]));
    time+=usecond();
    std::cout << GridLogMessage << "Z took "<<time<< " us"<<std::endl;

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
    std::cout << GridLogMessage << "ZxAd took "<<time<< " us"<<std::endl;

    //////////////////////////////////////
    // J(x) = 1 + Sum_k=1..N (-Zac)^k/(k+1)!
    //////////////////////////////////////
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
    std::cout << GridLogMessage << "Jx took "<<time<< " us"<<std::endl;

    //////////////////////////////////////
    // dJ(x)/dxe
    //////////////////////////////////////
    time=-usecond();
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
      for (int j = 20; j > 1; --j) {
	t3 = t2*(1.0 / (j + 1))  + aunit;
	dt3 = dt2*(1.0 / (j + 1));
	t2 = X * t3;
	dt2 = TRb * t3 + X * dt3;
      }
      dJdX[b] = -dt2; 
    }
    time+=usecond();
    std::cout << GridLogMessage << "dJx took "<<time<< " us"<<std::endl;
    /////////////////////////////////////////////////////////////////
    // Mask Umu for this link
    /////////////////////////////////////////////////////////////////
    time=-usecond();
    PlaqL = Ident;
    PlaqR = Utmp*adj(Cmu);
    ComputeNxy(PlaqL,PlaqR,NxxAd);
    time+=usecond();
    std::cout << GridLogMessage << "ComputeNxy took "<<time<< " us"<<std::endl;
    
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
    std::cout << GridLogMessage << "MpAdInv took "<<time<< " us"<<std::endl;
    
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

    Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx,FdetV);
    Fdet2_mu=FdetV;
    Fdet1_mu=Zero();
    
    for(int e =0 ; e<8 ; e++){
      LatticeComplexD tr(grid);
      ColourMatrix te;
      SU3::generator(e, te);
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

    RealD t4 = usecond();
    for(int nu=0;nu<Nd;nu++){

      if (nu!=mu) {
	///////////////// +ve nu /////////////////
	//     __
	//    |  |
	//    x==    // nu polarisation -- clockwise

	time=-usecond();
	PlaqL=Ident;

	PlaqR=(-rho)*Gimpl::CovShiftForward(Umu[nu], nu,
 	       Gimpl::CovShiftForward(Umu[mu], mu,
	         Gimpl::CovShiftBackward(Umu[nu], nu,
		   Gimpl::CovShiftIdentityBackward(Utmp, mu))));
	time+=usecond();
	std::cout << GridLogMessage << "PlaqLR took "<<time<< " us"<<std::endl;

	time=-usecond();
	dJdXe_nMpInv_y =   dJdXe_nMpInv;
	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_nu = transpose(Nxy)*dJdXe_nMpInv_y;
	time+=usecond();
	std::cout << GridLogMessage << "ComputeNxy (occurs 6x) took "<<time<< " us"<<std::endl;

	time=-usecond();
	PlaqR=(-1.0)*PlaqR;
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx,FdetV);
	Fdet2_nu = FdetV;
	time+=usecond();
	std::cout << GridLogMessage << "Compute_MpInvJx_dNxxSy (occurs 6x) took "<<time<< " us"<<std::endl;
	
	//    x==
	//    |  |
	//    .__|    // nu polarisation -- anticlockwise

	PlaqR=(rho)*Gimpl::CovShiftForward(Umu[nu], nu,
		      Gimpl::CovShiftBackward(Umu[mu], mu,
    	 	        Gimpl::CovShiftIdentityBackward(Umu[nu], nu)));

	PlaqL=Gimpl::CovShiftIdentityBackward(Utmp, mu);

	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,mu,-1);
	ComputeNxy(PlaqL, PlaqR,Nxy);
	Fdet1_nu = Fdet1_nu+transpose(Nxy)*dJdXe_nMpInv_y;
	

	MpInvJx_nu = Cshift(MpInvJx,mu,-1);
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_nu = Fdet2_nu+FdetV;
	
	///////////////// -ve nu /////////////////
	//  __
	// |  |
	// x==          // nu polarisation -- clockwise

	PlaqL=(rho)* Gimpl::CovShiftForward(Umu[mu], mu,
		       Gimpl::CovShiftForward(Umu[nu], nu,
			 Gimpl::CovShiftIdentityBackward(Utmp, mu)));

        PlaqR = Gimpl::CovShiftIdentityForward(Umu[nu], nu);

	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,1);
	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_nu = Fdet1_nu + transpose(Nxy)*dJdXe_nMpInv_y;

	MpInvJx_nu = Cshift(MpInvJx,nu,1);
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_nu = Fdet2_nu+FdetV;
	
	// x==
	// |  |
	// |__|         // nu polarisation

	PlaqL=(-rho)*Gimpl::CovShiftForward(Umu[nu], nu,
 	        Gimpl::CovShiftIdentityBackward(Utmp, mu));

	PlaqR=Gimpl::CovShiftBackward(Umu[mu], mu,
	        Gimpl::CovShiftIdentityForward(Umu[nu], nu));

	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,mu,-1);
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv_y,nu,1);

	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_nu = Fdet1_nu + transpose(Nxy)*dJdXe_nMpInv_y;

	MpInvJx_nu = Cshift(MpInvJx,mu,-1);
	MpInvJx_nu = Cshift(MpInvJx_nu,nu,1);
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_nu = Fdet2_nu+FdetV;

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
	PlaqL=(-rho)*Gimpl::CovShiftForward(Umu[mu], mu,
		      Gimpl::CovShiftBackward(Umu[nu], nu,
   		        Gimpl::CovShiftIdentityBackward(Utmp, mu)));

	PlaqR=Gimpl::CovShiftIdentityBackward(Umu[nu], nu);
	
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,-1);

	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_mu = Fdet1_mu + transpose(Nxy)*dJdXe_nMpInv_y;

	MpInvJx_nu = Cshift(MpInvJx,nu,-1);

	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_mu = Fdet2_mu+FdetV;

	//  __
	// "  |
	// x__|          // mu polarisation

	PlaqL=(-rho)*Gimpl::CovShiftForward(Umu[mu], mu,
		       Gimpl::CovShiftForward(Umu[nu], nu,
		 	 Gimpl::CovShiftIdentityBackward(Utmp, mu)));

        PlaqR=Gimpl::CovShiftIdentityForward(Umu[nu], nu);

	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,1);

	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_mu = Fdet1_mu + transpose(Nxy)*dJdXe_nMpInv_y;

	MpInvJx_nu = Cshift(MpInvJx,nu,1);

	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_mu = Fdet2_mu+FdetV;
	
      }
    }
    RealD t5 = usecond();

    Fdet1_mu = Fdet1_mu + transpose(NxxAd)*dJdXe_nMpInv;

    InsertForce(Fdet1,Fdet1_mu,mu);
    InsertForce(Fdet2,Fdet2_mu,mu);

    force= (-0.5)*( Fdet1 + Fdet2);
    RealD t1 = usecond();
    std::cout << GridLogMessage << " logDetJacobianForce level took "<<t1-t0<<" us "<<std::endl;
    std::cout << GridLogMessage << " logDetJacobianForce t3-t0 "<<t3a-t0<<" us "<<std::endl;
    std::cout << GridLogMessage << " logDetJacobianForce t4-t3 dJdXe_nMpInv "<<t4-t3a<<" us "<<std::endl;
    std::cout << GridLogMessage << " logDetJacobianForce t5-t4 mu nu loop "<<t5-t4<<" us "<<std::endl;
    std::cout << GridLogMessage << " logDetJacobianForce t1-t5 "<<t1-t5<<" us "<<std::endl;
    std::cout << GridLogMessage << " logDetJacobianForce level took "<<t1-t0<<" us "<<std::endl;
  }
  RealD logDetJacobianLevel(const GaugeField &U,int smr)
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
    AdjMatrixField  Ncb(grid);
    AdjMatrixField  Jac(grid);
    AdjMatrixField  Zac(grid);
    AdjMatrixField  mZac(grid);
    AdjMatrixField  X(grid);

    int mu= (smr/2) %Nd;

    auto mask=PeekIndex<LorentzIndex>(masks[smr],mu); // the cb mask

    //////////////////////////////////////////////////////////////////
    // Assemble the N matrix
    //////////////////////////////////////////////////////////////////
    // Computes ALL the staples -- could compute one only here
    this->StoutSmearing->BaseSmear(C, U);
    Cmu = peekLorentz(C, mu);
    Umu = peekLorentz(U, mu);
    Complex ci(0,1);
    for(int b=0;b<Ngen;b++) {
      SU3::generator(b, Tb);
      // Qlat Tb = 2i Tb^Grid
      Nb = (2.0)*Ta( ci*Tb * Umu * adj(Cmu));
      for(int c=0;c<Ngen;c++) {
	SU3::generator(c, Tc);
	auto tmp = -trace(ci*Tc*Nb); // Luchang's norm: (2Tc) (2Td) N^db = -2 delta cd N^db // - was important
	PokeIndex<ColourIndex>(Ncb,tmp,c,b); 
      }
    }      

    //////////////////////////////////////////////////////////////////
    // Assemble Luscher exp diff map J matrix 
    //////////////////////////////////////////////////////////////////
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

    //////////////////////////////////////
    // J(x) = 1 + Sum_k=1..N (-Zac)^k/(k+1)!
    //////////////////////////////////////
    X=1.0; 
    Jac = X;
    mZac = (-1.0)*Zac; 
    RealD kpfac = 1;
    for(int k=1;k<12;k++){
      X=X*mZac;
      kpfac = kpfac /(k+1);
      Jac = Jac + X * kpfac;
    }

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
    Complex result = sum(ln_det);
    return result.real();
  }
public:
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
    }
    return ln_det;
  }
  void logDetJacobianForce(GaugeField &force)
  {
    force =Zero();
    GaugeField force_det(force.Grid());

    if (this->smearingLevels > 0)
    {
      double start = usecond();

      GaugeLinkField tmp_mu(force.Grid());

      for (int ismr = this->smearingLevels - 1; ismr > 0; --ismr) {

	// remove U in UdSdU...
	for (int mu = 0; mu < Nd; mu++) {
	  tmp_mu = adj(peekLorentz(this->get_smeared_conf(ismr), mu)) * peekLorentz(force, mu);
	  pokeLorentz(force, tmp_mu, mu);
	}
	
      	// Propagate existing force
        force = this->AnalyticSmearedForce(force, this->get_smeared_conf(ismr - 1), ismr);

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

      force = this->AnalyticSmearedForce(force, *this->ThinLinks,0);

      for (int mu = 0; mu < Nd; mu++) {
	tmp_mu = peekLorentz(*this->ThinLinks, mu) * peekLorentz(force, mu);
	pokeLorentz(force, tmp_mu, mu);
      }

      force_det = Zero();

      logDetJacobianForceLevel(*this->ThinLinks,force_det,0);

      force = force + force_det;

      force=Ta(force); // Ta
      
      double end = usecond();
      double time = (end - start)/ 1e3;
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
    GridBase* grid = GaugeK.Grid();
    GaugeField C(grid), SigmaK(grid), iLambda(grid);
    GaugeField SigmaKPrimeA(grid);
    GaugeField SigmaKPrimeB(grid);
    GaugeLinkField iLambda_mu(grid);
    GaugeLinkField iQ(grid), e_iQ(grid);
    GaugeLinkField SigmaKPrime_mu(grid);
    GaugeLinkField GaugeKmu(grid), Cmu(grid);
    
    this->StoutSmearing->BaseSmear(C, GaugeK);
    SigmaK = Zero();
    iLambda = Zero();

    SigmaKPrimeA = SigmaKPrime;
    ApplyMask(SigmaKPrimeA,level);
    SigmaKPrimeB = SigmaKPrime - SigmaKPrimeA;
    
    // Could get away with computing only one polarisation here
    // int mu= (smr/2) %Nd;
    // SigmaKprime_A has only one component
    for (int mu = 0; mu < Nd; mu++)
    {
      Cmu = peekLorentz(C, mu);
      GaugeKmu = peekLorentz(GaugeK, mu);
      SigmaKPrime_mu = peekLorentz(SigmaKPrimeA, mu);
      iQ = Ta(Cmu * adj(GaugeKmu));
      this->set_iLambda(iLambda_mu, e_iQ, iQ, SigmaKPrime_mu, GaugeKmu);
      pokeLorentz(SigmaK, SigmaKPrime_mu * e_iQ + adj(Cmu) * iLambda_mu, mu);
      pokeLorentz(iLambda, iLambda_mu, mu);
    }
    this->StoutSmearing->derivative(SigmaK, iLambda,GaugeK);  // derivative of SmearBase

    ////////////////////////////////////////////////////////////////////////////////////
    // propagate the rest of the force as identity map, just add back
    ////////////////////////////////////////////////////////////////////////////////////
    SigmaK = SigmaK+SigmaKPrimeB;

    return SigmaK;
  }

public:

  /* Standard constructor */
  SmearedConfigurationMasked(GridCartesian* _UGrid, unsigned int Nsmear, Smear_Stout<Gimpl>& Stout)
    : SmearedConfiguration<Gimpl>(_UGrid, Nsmear,Stout)
  {
    assert(Nsmear%(2*Nd)==0); // Or multiply by 8??

    // was resized in base class
    assert(this->SmearedSet.size()==Nsmear);
    
    GridRedBlackCartesian * UrbGrid;
    UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(_UGrid);
    LatticeComplex one(_UGrid); one = ComplexD(1.0,0.0);
    LatticeComplex tmp(_UGrid);

    for (unsigned int i = 0; i < this->smearingLevels; ++i) {

      masks.push_back(*(new LatticeLorentzComplex(_UGrid)));

      int mu= (i/2) %Nd;
      int cb= (i%2);
      LatticeComplex tmpcb(UrbGrid);
	
      masks[i]=Zero();
      ////////////////////
      // Setup the mask
      ////////////////////
      tmp = Zero();
      pickCheckerboard(cb,tmpcb,one);
      setCheckerboard(tmp,tmpcb);
      PokeIndex<LorentzIndex>(masks[i],tmp, mu);
	
    }
    delete UrbGrid;
  }
  
  virtual void smeared_force(GaugeField &SigmaTilde) 
  {
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

