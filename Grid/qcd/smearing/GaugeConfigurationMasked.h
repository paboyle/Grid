/*!
  @file GaugeConfiguration.h
  @brief Declares the GaugeConfiguration class
*/
#pragma once

NAMESPACE_BEGIN(Grid);


template<class T> void Dump(Lattice<T> & lat, std::string s,Coordinate site = Coordinate({0,0,0,0}))
{
  typename T::scalar_object tmp;
  peekSite(tmp,lat,site);
  std::cout << " GRID "<<s<<" "<<tmp<<std::endl;
}
/*!
  @brief Smeared configuration masked container
  Modified for a multi-subset smearing (aka Luscher Flowed HMC)
*/
template <class Gimpl>
class SmearedConfigurationMasked
{
public:
  INHERIT_GIMPL_TYPES(Gimpl);

private:
  const unsigned int smearingLevels;
  Smear_Stout<Gimpl> *StoutSmearing;
  std::vector<GaugeField> SmearedSet;
  std::vector<LatticeLorentzComplex> masks;

  ///////////////////////////////////////////////////
  // Additions
  // Could just set masks to 1 in case of oldsmearing
  // and share the code. Pass in a mask creation object with alternate constructor -- SmearMaskMaker
  ///////////////////////////////////////////////////
  typedef typename SU3Adjoint::AMatrix AdjMatrix;
  typedef typename SU3Adjoint::LatticeAdjMatrix  AdjMatrixField;
  typedef typename SU3Adjoint::LatticeAdjVector  AdjVectorField;
  /*
      set_zero(f_det_basis);
      for (int a = 0; a < 8; ++a) {
        const ColorMatrix uc = uc_pre * ts[a] * uc_post;
        for (int c = 0; c < 8; ++c) {
          const ColorMatrix d_n = make_tr_less_anti_herm_matrix(ts[c] * uc);
          const array<double, 8> d_n_b =
              basis_projection_anti_hermitian_matrix(d_n);
          for (int b = 0; b < 8; ++b) {
	    f_det_basis[a] += n_e_mp_inv_j_x_mat(c, b) * d_n_b[b];
          }
        }
      }
  */
  void AdjointDeriv2(const GaugeLinkField &PlaqL,const GaugeLinkField &PlaqR, AdjMatrixField MpInvJx,AdjVectorField &Fdet2 )
  {
    GaugeLinkField UtaU(PlaqL.Grid());
    GaugeLinkField D(PlaqL.Grid());
    GaugeLinkField aa(PlaqL.Grid());
    AdjMatrixField Dbc(PlaqL.Grid());
    LatticeComplex tmp(PlaqL.Grid());
    const int Ngen = SU3Adjoint::Dimension;
    Complex ci(0,1);
    ColourMatrix   ta,tb,tc;
    
    for(int a=0;a<Ngen;a++) {
      SU3::generatorQlat(a, ta);
      // Qlat Tb = 2i Tb^Grid
      UtaU= 2.0*ci*adj(PlaqL)*ta*PlaqR;
      for(int c=0;c<Ngen;c++) {
	SU3::generatorQlat(c, tc);
	D = Ta( (2.0)*ci*tc *UtaU);
	//	Dump(D," (-4.0)*Ta( tc *UtaU)");
	for(int b=0;b<Ngen;b++){
	  SU3::generatorQlat(b, tb);
	  tmp =-trace(ci*tb*D); 
	  PokeIndex<ColourIndex>(Dbc,tmp,b,c);  // Adjoint rep
	  //	  Dump(tmp," -trace(ci*tb*D)");
	}
      }
      //      Dump(Dbc," Dbc ");
      //      Dump(MpInvJx," MpInvJx ");
      tmp = trace(MpInvJx * Dbc);
      //      Dump(tmp,"  trace(MpInvJx * Dbc) ");
      PokeIndex<ColourIndex>(Fdet2,tmp,a);
    }
    Dump(Fdet2," Fdet2 ");
  }
  
  void AdjointDeriv(const GaugeLinkField &PlaqL,const GaugeLinkField &PlaqR,AdjMatrixField &NxAd)
  {
    GaugeLinkField Nx(PlaqL.Grid());
    const int Ngen = SU3Adjoint::Dimension;
    Complex ci(0,1);
    ColourMatrix   tb;
    ColourMatrix   tc;
    for(int b=0;b<Ngen;b++) {
      SU3::generatorQlat(b, tb);
      // Qlat Tb = 2i Tb^Grid
      Nx = (2.0)*Ta( adj(PlaqL)*ci*tb * PlaqR );
      for(int c=0;c<Ngen;c++) {
	SU3::generatorQlat(c, tc);
	auto tmp =closure( -trace(ci*tc*Nx)); // Luchang's norm: (2Tc) (2Td) N^db = -2 delta cd N^db // - was important
	PokeIndex<ColourIndex>(NxAd,tmp,c,b);  // Adjoint rep
      }
    }
  }
  void ApplyMask(GaugeField &U,int smr)
  {
    LatticeComplex tmp(UGrid);
    GaugeLinkField Umu(UGrid);
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
    conformable(grid,UGrid);
    ColourMatrix   tb;
    ColourMatrix   tc;
    ColourMatrix   ta;
    GaugeField C(grid);
    GaugeField Umsk(grid);
    std::vector<GaugeLinkField> Umu(Nd,grid);
    GaugeLinkField Cmu(grid); // U and staple; C contains factor of epsilon
    GaugeLinkField Zx(grid);  // U times Staple, contains factor of epsilon
    GaugeLinkField Nx(grid);  // Nxx fundamental space
    GaugeLinkField Utmp(grid);
    GaugeLinkField PlaqL(grid);
    GaugeLinkField PlaqR(grid);
    const int Ngen = SU3Adjoint::Dimension;
    AdjMatrix TRb;
    ColourMatrix Ident;
    LatticeComplex  cplx(grid); 
    AdjVectorField  AlgV(grid); 
    AdjVectorField  AlgVtmp(grid); 
    AdjMatrixField  MpAd(grid);    // Mprime luchang's notes
    AdjMatrixField  MpAdInv(grid); // Mprime inverse
    AdjMatrixField  NxAd(grid);    // Nxx in adjoint space
    AdjMatrixField  JxAd(grid);     
    AdjMatrixField  ZxAd(grid);
    AdjMatrixField  mZxAd(grid);
    AdjMatrixField  X(grid);
    Complex ci(0,1);

    Ident = ComplexD(1.0);
    for(int d=0;d<Nd;d++){
      Umu[d] = peekLorentz(U, d);
    }
    int mu= (smr/2) %Nd;

    auto mask=PeekIndex<LorentzIndex>(masks[smr],mu); // the cb mask

    // Mask the gauge field
    Umsk = U;
    ApplyMask(Umsk,smr);
    Utmp = peekLorentz(Umsk,mu);

    //////////////////////////////////////////////////////////////////
    // Assemble the N matrix
    //////////////////////////////////////////////////////////////////
    // Computes ALL the staples -- could compute one only here
    StoutSmearing->BaseSmear(C, U);
    double rho=0.1;
    Cmu = peekLorentz(C, mu);
    Dump(Cmu,std::string(" Cmu "));

    //////////////////////////////////////////////////////////////////
    // Assemble Luscher exp diff map J matrix 
    //////////////////////////////////////////////////////////////////
    // Ta so Z lives in Lie algabra
    Zx  = Ta(Cmu * adj(Umu[mu]));
    //    Dump(Zx,std::string("Zx"));

    // Move Z to the Adjoint Rep == make_adjoint_representation
    ZxAd = Zero();
    for(int b=0;b<8;b++) {
      // Adj group sets traceless antihermitian T's -- Guido, really????
      // Is the mapping of these the same? Same structure constants
      // Might never have been checked.
      SU3::generatorQlat(b, tb);         // Fund group sets traceless hermitian T's
      SU3Adjoint::generatorQlat(b,TRb);
      TRb=-TRb;
      cplx = 2.0*trace(ci*tb*Zx); // my convention 1/2 delta ba
      ZxAd = ZxAd + cplx * TRb; // is this right? YES - Guido used Anti herm Ta's and with bloody wrong sign.
    }
    Dump(ZxAd,std::string("ZxAd"));

    //////////////////////////////////////
    // J(x) = 1 + Sum_k=1..N (-Zac)^k/(k+1)!
    //////////////////////////////////////
    X=1.0; 
    JxAd = X;
    mZxAd = (-1.0)*ZxAd; 
    RealD kpfac = 1;
    for(int k=1;k<12;k++){
      X=X*mZxAd;
      kpfac = kpfac /(k+1);
      JxAd = JxAd + X * kpfac;
    }
    Dump(JxAd,std::string("JxAd"));

    //////////////////////////////////////
    // dJ(x)/dxe = d/dxe Sum_k=1..N X^k/(k+1)!
    //           = 1/2! te
    //           + 1/3! (te x + x te ) )
    //           + 1/4! (te x x + x te x + x x te )
    //           + 1/5! (te x x x + x te x x + x x te x + x x x te )
    // Iterate:
    // teX_n = teX_{n-1} x
    // S_n = x S_{n-1} + teX_{n}
    //////////////////////////////////////
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
      SU3Adjoint::generatorQlat(b, TRb); //dt2
      Dump(ZxAd,std::string("ZxAd"));
      X  = (-1.0)*ZxAd; 
      t2 = X;
      dt2 = TRb;
      for (int j = 20; j > 1; --j) {
	t3 = t2*(1.0 / (j + 1))  + aunit;
	dt3 = dt2*(1.0 / (j + 1));
	t2 = X * t3;
	dt2 = TRb * t3 + X * dt3;
      }
      //      dt3 = .5 * dt2;
      dJdX[b] = -dt2; // sign and 2x ?
      Dump(dJdX[b],std::string("dJdX"));
    }
      /*
      X = (-1.0)*ZxAd; //x=t2
      // n=1 point
      tbXn=TRb;
      sumXtbX= TRb;
      RealD kpfac = 1.0/2.0;
      dJdX[b] = sumXtbX *kpfac;
      for(int k=2;k<12;k++){
	kpfac = kpfac /(k+1);
	tbXn = tbXn*X;
	sumXtbX = X*sumXtbX + tbXn;
	dJdX[b] = dJdX[b] + sumXtbX *kpfac;
      }
      */
    /////////////////////////////////////////////////////////////////
    // Mask Umu for this link
    /////////////////////////////////////////////////////////////////
    //      Nx = (2.0)*Ta( adj(PlaqL)*ci*tb * PlaqR );
    PlaqL = Ident;
    PlaqR = Utmp*adj(Cmu);
    AdjointDeriv(PlaqL,PlaqR,NxAd);
    Dump(NxAd,std::string("NxAd"));

    
    ////////////////////////////
    // Mab
    ////////////////////////////
    //    Mab = Complex(1.0,0.0);
    //    Mab = Mab - Jac * Ncb;

    MpAd = Complex(1.0,0.0);
    MpAd = MpAd - JxAd * NxAd;
    Dump(MpAd,std::string("MpAd"));

    /////////////////////////
    // invert the 8x8
    /////////////////////////
    MpAdInv = Inverse(MpAd);
    Dump(MpAdInv,std::string("MpAdInv"));
    
    /////////////////////////////////////////////////////////////////
    // Alternate way of creating
    // May need to keep the +nu and -nu plaq fields distinct to stop collision
    /////////////////////////////////////////////////////////////////

    AdjMatrixField nMpInv(grid);
    nMpInv= NxAd *MpAdInv;
    Dump(nMpInv," nMpInv ");

    AdjMatrixField MpInvJx(grid);
    AdjMatrixField MpInvJx_nu(grid);
    MpInvJx = (-1.0)*MpAdInv * JxAd;// rho is on the plaq factor
    Dump(MpInvJx," MpInvJx ");

    AdjVectorField  FdetV(grid);
    AdjVectorField  FdetV2(grid);
    AdjVectorField  FdetV2_mu(grid);
    // First shot at the Fdet2
    AdjointDeriv2(PlaqL,PlaqR,MpInvJx,FdetV);
    FdetV2_mu=FdetV;
    Dump(FdetV,std::string(" FdetV2xx_mu "));
    
    for(int e =0 ; e<8 ; e++){
      LatticeComplexD tr(grid);
      ColourMatrix te;
      SU3::generatorQlat(e, te);
      tr = trace(dJdX[e] * nMpInv);
      pokeColour(AlgV,tr,e);
      //      std::cout << " ***** tr " <<e<<std::endl;
      //      Dump(tr,std::string("tr"));
    }
    // Mask it off
    auto tmp=PeekIndex<LorentzIndex>(masks[smr],mu);
    AlgV = AlgV*tmp;
    Dump(AlgV,std::string("AlgV"));
    //    AlgV needs to multiply:
    //       NxAd (site local)                                            (1)
    //       NfmuPlus (site local) one site forward  in each nu direction (3)
    //       NfmuMinus (site local) one site backward in each nu direction(3)
    //       Nfnu (site local) 0,0 ; +mu,0; 0,-nu; +mu-nu   [ 3x4 = 12]
    // 19 terms.
    AdjVectorField  AlgVmu_p(grid); AlgVmu_p=Zero();
    AdjVectorField  AlgVmu_m(grid); AlgVmu_m=Zero();
    AdjVectorField  AlgVnu(grid);   AlgVnu=Zero();

    GaugeLinkField FmuPlus(grid);
    GaugeLinkField FmuMinus(grid);
    GaugeLinkField Fnu(grid);
    GaugeLinkField Fnumu(grid);
    std::vector<AdjMatrixField> Nfperp(Nd,grid); // Why needed vs nu?
    AdjMatrixField NfmuPlus(grid);
    AdjMatrixField NfmuMinus(grid);
    for(int d=0;d<Nd;d++){
      Nfperp[d] = Zero();
    }
    NfmuPlus=Zero();
    NfmuMinus=Zero();
    //////////////////////////////////////////////////////
    // six force inserts x 3 directions from OTHER links y!=x || mu!=nu
    //
    // To avoid collison of x-y pairs need to keep +ve nu and -ve nu separate, at least for mu force
    // FmuPlus, FmuMinus and Fnu/Fperp
    //////////////////////////////////////////////////////
    GaugeField Fdet1(grid);
    GaugeField Fdet2(grid);
    GaugeLinkField Fdet1_pol(grid);
    GaugeLinkField Fdet2_pol(grid);
    
    AdjVectorField  FdetV_acc(grid); 
    AdjVectorField  FdetV2_acc(grid); 
    AdjVectorField  FdetV_mu(grid);  FdetV_mu=Zero();
    for(int nu=0;nu<Nd;nu++){

      if (nu!=mu) {
	LatticeComplexD tr(grid);

	int Lnu = grid->GlobalDimensions()[nu];
	Coordinate coormu({0,0,0,0});	coormu[mu]=1;
	Coordinate coornu({0,0,0,0});	coornu[nu]=1;
	Coordinate coornnu({0,0,0,0});	coornnu[nu]=Lnu-1;
	Coordinate coormunu({0,0,0,0});	 coormunu[mu]=1; coormunu[nu]=1;
	Coordinate coormunnu({0,0,0,0}); coormunnu[mu]=1; coormunnu[nu]=Lnu-1;

	///////////////// +ve nu /////////////////
	//     __
	//    |  |
	//    x==    // nu polarisation -- clockwise
	std::cout << " ********************* "<<std::endl;
	std::cout << "   nu+ = "<<nu<<std::endl;
	std::cout << " ********************* "<<std::endl;
	PlaqL=Ident;
	PlaqR=(-rho)*Gimpl::CovShiftForward(Umu[nu], nu,
 	       Gimpl::CovShiftForward(Umu[mu], mu,
	         Gimpl::CovShiftBackward(Umu[nu], nu,
		   Gimpl::CovShiftIdentityBackward(Utmp, mu))));

	//      Nx = -2.0*ci*trace(tc*Ta( adj(PlaqL)*ci*tb * PlaqR )));
	AdjointDeriv(PlaqL,PlaqR,Nfperp[nu]);

	AlgVnu =   AlgV;

	PlaqR=(-1.0)*PlaqR;
	FdetV = transpose(Nfperp[nu])*AlgVnu;
	Dump(FdetV,std::string("FdetVxy_nu ; y=x ")); // OK

	FdetV_acc = FdetV;

	AdjointDeriv2(PlaqL,PlaqR,MpInvJx,FdetV2);
	FdetV2_acc = FdetV2;
	Dump(FdetV2,std::string("FdetV2xy_nu ; y=x ")); 
	
	//    x==
	//    |  |
	//    .__|    // nu polarisation -- anticlockwise

	//      Nx = (2.0)*Ta( adj(PlaqL)*ci*tb * PlaqR );
	PlaqR=(rho)*Gimpl::CovShiftForward(Umu[nu], nu,
		      Gimpl::CovShiftBackward(Umu[mu], mu,
    	 	        Gimpl::CovShiftIdentityBackward(Umu[nu], nu)));
	PlaqL=Gimpl::CovShiftIdentityBackward(Utmp, mu);

	AdjointDeriv(PlaqL, PlaqR,Nfperp[nu]);

	AlgVnu = Cshift(AlgV,mu,-1);
	
	FdetV = transpose(Nfperp[nu])*AlgVnu;
	Dump(FdetV,std::string("FdetVxy_nu ; +mu "),coormu); // OK

	
	FdetV_acc = FdetV_acc + FdetV;

	MpInvJx_nu = Cshift(MpInvJx,mu,-1);
	AdjointDeriv2(PlaqL,PlaqR,MpInvJx_nu,FdetV2);
	FdetV2_acc = FdetV2_acc+FdetV2;
	Dump(FdetV2,std::string("FdetV2xy_nu ; +mu "),coormu); 
	
	///////////////// -ve nu /////////////////
	//  __
	// |  |
	// x==          // nu polarisation -- clockwise
	std::cout << " ********************* "<<std::endl;
	std::cout << "   nu- = "<<nu<<std::endl;
	std::cout << " ********************* "<<std::endl;

	PlaqL=(rho)* Gimpl::CovShiftForward(Umu[mu], mu,
		       Gimpl::CovShiftForward(Umu[nu], nu,
			 Gimpl::CovShiftIdentityBackward(Utmp, mu)));
        PlaqR = Gimpl::CovShiftIdentityForward(Umu[nu], nu);
	AdjointDeriv(PlaqL,PlaqR,Nfperp[nu]);
	
	AlgVnu = Cshift(AlgV,nu,1);

	FdetV = transpose(Nfperp[nu])*AlgVnu;

	Dump(FdetV,std::string("FdetVxy_nu ; -nu"),coornnu); 

	FdetV_acc = FdetV_acc + FdetV;

	MpInvJx_nu = Cshift(MpInvJx,nu,1);
	AdjointDeriv2(PlaqL,PlaqR,MpInvJx_nu,FdetV2);
	FdetV2_acc = FdetV2_acc+FdetV2;
	Dump(FdetV2,std::string("FdetV2xy_nu ; -nu "),coornnu); 
	
	// x==
	// |  |
	// |__|         // nu polarisation

	PlaqL=(-rho)*Gimpl::CovShiftForward(Umu[nu], nu,
 	        Gimpl::CovShiftIdentityBackward(Utmp, mu));
	PlaqR=Gimpl::CovShiftBackward(Umu[mu], mu,
	        Gimpl::CovShiftIdentityForward(Umu[nu], nu));

	AdjointDeriv(PlaqL,PlaqR,Nfperp[nu]);

	AlgVnu = Cshift(AlgV,mu,-1);
	AlgVnu = Cshift(AlgVnu,nu,1);

	FdetV = transpose(Nfperp[nu])*AlgVnu;
	Dump(FdetV,std::string("FdetVxy_nu; -nu +mu"),coormunnu);

	FdetV_acc = FdetV_acc + FdetV;

	MpInvJx_nu = Cshift(MpInvJx,mu,-1);
	MpInvJx_nu = Cshift(MpInvJx_nu,nu,1);
	AdjointDeriv2(PlaqL,PlaqR,MpInvJx_nu,FdetV2);
	FdetV2_acc = FdetV2_acc+FdetV2;
	Dump(FdetV2,std::string("FdetV2xy_nu ; -nu +mu "),coormunu); 
	
	// Set up the determinant force contribution in 3x3 algebra basis
	Fdet1_pol=Zero();
	for(int e=0;e<8;e++){
	  ColourMatrix te;
	  SU3::generatorQlat(e, te);
	  auto tmp=peekColour(FdetV_acc,e);
	  Fdet1_pol=Fdet1_pol + ci*tmp*te; // but norm of te is different.. why?
	}
	pokeLorentz(Fdet1, Fdet1_pol, nu);

	Fdet2_pol=Zero();
	for(int e=0;e<8;e++){
	  ColourMatrix te;
	  SU3::generatorQlat(e, te);
	  auto tmp=peekColour(FdetV2_acc,e);
	  Fdet2_pol=Fdet2_pol + ci*tmp*te; // but norm of te is different.. why?
	}
	pokeLorentz(Fdet2, Fdet2_pol, nu);
	
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
	
	AdjointDeriv(PlaqL,PlaqR,NfmuPlus);

	AlgVnu = Cshift(AlgV,nu,-1);

	FdetV = transpose(NfmuPlus)*AlgVnu;
	Dump(FdetV,std::string("FdetV mu ; +nu "),coornu);

	FdetV_mu = FdetV_mu + FdetV;

	MpInvJx_nu = Cshift(MpInvJx,nu,-1);
	AdjointDeriv2(PlaqL,PlaqR,MpInvJx_nu,FdetV2);
	FdetV2_mu = FdetV2_mu+FdetV2;
	Dump(FdetV2,std::string("FdetV2xy_mu ; +nu "),coornu); 

	//  __
	// "  |
	// x__|          // mu polarisation

	PlaqL=(-rho)*Gimpl::CovShiftForward(Umu[mu], mu,
		       Gimpl::CovShiftForward(Umu[nu], nu,
		 	 Gimpl::CovShiftIdentityBackward(Utmp, mu)));
        PlaqR=Gimpl::CovShiftIdentityForward(Umu[nu], nu);

	AdjointDeriv(PlaqL,PlaqR,NfmuMinus);

	AlgVnu = Cshift(AlgV,nu,1);

	FdetV = transpose(NfmuMinus)*AlgVnu;

	Dump(FdetV,std::string("FdetV_xy mu ; -nu "),coornnu);

	FdetV_mu = FdetV_mu + FdetV;

	MpInvJx_nu = Cshift(MpInvJx,nu,1);
	AdjointDeriv2(PlaqL,PlaqR,MpInvJx_nu,FdetV2);
	FdetV2_mu = FdetV2_mu+FdetV2;
	Dump(FdetV2,std::string("FdetV2xy_mu ; -nu "),coornnu); 
	
      }
    }

    FdetV = transpose(NxAd)*AlgV;
    Dump(FdetV,std::string(" FdetVxx_mu "));
    FdetV_mu = FdetV_mu + FdetV;
    Fdet1_pol=Zero();
    for(int e=0;e<8;e++){
      ColourMatrix te;
      SU3::generatorQlat(e, te);
      auto tmp=peekColour(FdetV_mu,e);
      Fdet1_pol=Fdet1_pol + ci*tmp*te; // but norm of te is different.. why?
    }
    pokeLorentz(Fdet1, Fdet1_pol, mu);

    Fdet2_pol=Zero();
    for(int e=0;e<8;e++){
      ColourMatrix te;
      SU3::generatorQlat(e, te);
      auto tmp=peekColour(FdetV2_mu,e);
      Fdet2_pol=Fdet2_pol + ci*tmp*te; // but norm of te is different.. why?
    }
    pokeLorentz(Fdet2, Fdet2_pol, mu);

    force = Fdet1 + Fdet2;
    //    force = Fdet1;
    //    force = Fdet2;
    
  }
  RealD logDetJacobianLevel(const GaugeField &U,int smr)
  {
    GridBase* grid = U.Grid();
    conformable(grid,UGrid);
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
    StoutSmearing->BaseSmear(C, U);
    Cmu = peekLorentz(C, mu);
    Umu = peekLorentz(U, mu);
    Complex ci(0,1);
    for(int b=0;b<Ngen;b++) {
      SU3::generatorQlat(b, Tb);
      // Qlat Tb = 2i Tb^Grid
      Nb = (2.0)*Ta( ci*Tb * Umu * adj(Cmu));
      for(int c=0;c<Ngen;c++) {
	SU3::generatorQlat(c, Tc);
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
      SU3::generatorQlat(b, Tb);         // Fund group sets traceless hermitian T's
      SU3Adjoint::generatorQlat(b,TRb);
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
    if (smearingLevels > 0)
    {
      for (int ismr = smearingLevels - 1; ismr > 0; --ismr) {
	ln_det+= logDetJacobianLevel(get_smeared_conf(ismr-1),ismr);
      }
      ln_det +=logDetJacobianLevel(*ThinLinks,0);
    }
    return ln_det;
  }
  void logDetJacobianForce(GaugeField &force)
  {
    RealD ln_det = 0;
    if (smearingLevels > 0)
    {
      for (int ismr = smearingLevels - 1; ismr > 0; --ismr) {
	ln_det+= logDetJacobianForceLevel(get_smeared_conf(ismr-1),force,ismr);
      }
      ln_det +=logDetJacobianForeceLevel(*ThinLinks,force,0);
    }
  }

private:
  // Member functions
  //====================================================================
  void fill_smearedSet(GaugeField &U)
  {
    ThinLinks = &U;  // attach the smearing routine to the field U

    // check the pointer is not null
    if (ThinLinks == NULL)
      std::cout << GridLogError << "[SmearedConfigurationMasked] Error in ThinLinks pointer\n";

    if (smearingLevels > 0)
    {
      std::cout << GridLogMessage << "[SmearedConfigurationMasked] Filling SmearedSet\n";
      GaugeField previous_u(ThinLinks->Grid());

      GaugeField smeared_A(ThinLinks->Grid());
      GaugeField smeared_B(ThinLinks->Grid());

      previous_u = *ThinLinks;
      for (int smearLvl = 0; smearLvl < smearingLevels; ++smearLvl)
      {
        StoutSmearing->smear(smeared_A, previous_u);
	ApplyMask(smeared_A,smearLvl);
	smeared_B = previous_u;
	ApplyMask(smeared_B,smearLvl);
	// Replace only the masked portion
	SmearedSet[smearLvl] = previous_u-smeared_B + smeared_A;
        previous_u = SmearedSet[smearLvl];

        // For debug purposes
        RealD impl_plaq = WilsonLoops<Gimpl>::avgPlaquette(previous_u);
        std::cout << GridLogMessage << "[SmearedConfigurationMasked] Plaq: " << impl_plaq << std::endl;
      }
    }
  }
  //====================================================================
  GaugeField AnalyticSmearedForce(const GaugeField& SigmaKPrime,
                                  const GaugeField& GaugeK,int level) 
  {
    GridBase* grid = GaugeK.Grid();
    conformable(grid,UGrid);
    GaugeField C(grid), SigmaK(grid), iLambda(grid);
    GaugeField SigmaKPrimeA(grid);
    GaugeField SigmaKPrimeB(grid);
    GaugeLinkField iLambda_mu(grid);
    GaugeLinkField iQ(grid), e_iQ(grid);
    GaugeLinkField SigmaKPrime_mu(grid);
    GaugeLinkField GaugeKmu(grid), Cmu(grid);

    StoutSmearing->BaseSmear(C, GaugeK);
    SigmaK = Zero();
    iLambda = Zero();

    SigmaK = Zero();

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
      set_iLambda(iLambda_mu, e_iQ, iQ, SigmaKPrime_mu, GaugeKmu);
      pokeLorentz(SigmaK, SigmaKPrime_mu * e_iQ + adj(Cmu) * iLambda_mu, mu);
      pokeLorentz(iLambda, iLambda_mu, mu);
    }
    StoutSmearing->derivative(SigmaK, iLambda,GaugeK);  // derivative of SmearBase

    ////////////////////////////////////////////////////////////////////////////////////
    // propagate the rest of the force as identity map, just add back
    ////////////////////////////////////////////////////////////////////////////////////
    SigmaK = SigmaK+SigmaKPrimeB;
    return SigmaK;
  }


  ////////////////////////////////////////
  // INHERIT THESE
  ////////////////////////////////////////
  
  /*! @brief Returns smeared configuration at level 'Level' */
  const GaugeField &get_smeared_conf(int Level) const
  {
    return SmearedSet[Level];
  }

  // Duplicates code that is in GaugeConfiguration.h
  // Should inherit or share.
  //====================================================================
  void set_iLambda(GaugeLinkField& iLambda, GaugeLinkField& e_iQ,
                   const GaugeLinkField& iQ, const GaugeLinkField& Sigmap,
                   const GaugeLinkField& GaugeK) const 
  {
    GridBase* grid = iQ.Grid();
    GaugeLinkField iQ2(grid), iQ3(grid), B1(grid), B2(grid), USigmap(grid);
    GaugeLinkField unity(grid);
    unity = 1.0;

    LatticeComplex u(grid), w(grid);
    LatticeComplex f0(grid), f1(grid), f2(grid);
    LatticeComplex xi0(grid), xi1(grid), tmp(grid);
    LatticeComplex u2(grid), w2(grid), cosw(grid);
    LatticeComplex emiu(grid), e2iu(grid), qt(grid), fden(grid);
    LatticeComplex r01(grid), r11(grid), r21(grid), r02(grid), r12(grid);
    LatticeComplex r22(grid), tr1(grid), tr2(grid);
    LatticeComplex b10(grid), b11(grid), b12(grid), b20(grid), b21(grid),
      b22(grid);
    LatticeComplex LatticeUnitComplex(grid);

    LatticeUnitComplex = 1.0;

    // Exponential
    iQ2 = iQ * iQ;
    iQ3 = iQ * iQ2;
    StoutSmearing->set_uw(u, w, iQ2, iQ3);
    StoutSmearing->set_fj(f0, f1, f2, u, w);
    e_iQ = f0 * unity + timesMinusI(f1) * iQ - f2 * iQ2;

    // Getting B1, B2, Gamma and Lambda
    // simplify this part, reduntant calculations in set_fj
    xi0 = StoutSmearing->func_xi0(w);
    xi1 = StoutSmearing->func_xi1(w);
    u2 = u * u;
    w2 = w * w;
    cosw = cos(w);

    emiu = cos(u) - timesI(sin(u));
    e2iu = cos(2.0 * u) + timesI(sin(2.0 * u));

    r01 = (2.0 * u + timesI(2.0 * (u2 - w2))) * e2iu +
      emiu * ((16.0 * u * cosw + 2.0 * u * (3.0 * u2 + w2) * xi0) +
	      timesI(-8.0 * u2 * cosw + 2.0 * (9.0 * u2 + w2) * xi0));

    r11 = (2.0 * LatticeUnitComplex + timesI(4.0 * u)) * e2iu +
      emiu * ((-2.0 * cosw + (3.0 * u2 - w2) * xi0) +
	      timesI((2.0 * u * cosw + 6.0 * u * xi0)));

    r21 =
      2.0 * timesI(e2iu) + emiu * (-3.0 * u * xi0 + timesI(cosw - 3.0 * xi0));

    r02 = -2.0 * e2iu +
      emiu * (-8.0 * u2 * xi0 +
	      timesI(2.0 * u * (cosw + xi0 + 3.0 * u2 * xi1)));

    r12 = emiu * (2.0 * u * xi0 + timesI(-cosw - xi0 + 3.0 * u2 * xi1));

    r22 = emiu * (xi0 - timesI(3.0 * u * xi1));

    fden = LatticeUnitComplex / (2.0 * (9.0 * u2 - w2) * (9.0 * u2 - w2));

    b10 = 2.0 * u * r01 + (3.0 * u2 - w2) * r02 - (30.0 * u2 + 2.0 * w2) * f0;
    b11 = 2.0 * u * r11 + (3.0 * u2 - w2) * r12 - (30.0 * u2 + 2.0 * w2) * f1;
    b12 = 2.0 * u * r21 + (3.0 * u2 - w2) * r22 - (30.0 * u2 + 2.0 * w2) * f2;

    b20 = r01 - (3.0 * u) * r02 - (24.0 * u) * f0;
    b21 = r11 - (3.0 * u) * r12 - (24.0 * u) * f1;
    b22 = r21 - (3.0 * u) * r22 - (24.0 * u) * f2;

    b10 *= fden;
    b11 *= fden;
    b12 *= fden;
    b20 *= fden;
    b21 *= fden;
    b22 *= fden;

    B1 = b10 * unity + timesMinusI(b11) * iQ - b12 * iQ2;
    B2 = b20 * unity + timesMinusI(b21) * iQ - b22 * iQ2;
    USigmap = GaugeK * Sigmap;

    tr1 = trace(USigmap * B1);
    tr2 = trace(USigmap * B2);

    GaugeLinkField QUS = iQ * USigmap;
    GaugeLinkField USQ = USigmap * iQ;

    GaugeLinkField iGamma = tr1 * iQ - timesI(tr2) * iQ2 +
      timesI(f1) * USigmap + f2 * QUS + f2 * USQ;

    iLambda = Ta(iGamma);
  }

  //====================================================================
public:
  GaugeField* ThinLinks; /* Pointer to the thin links configuration */
  ////////////////////////
  // Derived class
  ////////////////////////
  GridRedBlackCartesian * UrbGrid;
  GridCartesian         * UGrid;
  /* Standard constructor */
  SmearedConfigurationMasked(GridCartesian* _UGrid, unsigned int Nsmear, Smear_Stout<Gimpl>& Stout,bool domask=false)
    : UGrid(_UGrid), smearingLevels(Nsmear), StoutSmearing(&Stout), ThinLinks(NULL)
  {
    if(domask) assert(Nsmear%(2*Nd)==0); // Or multiply by 8??
    
    UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(_UGrid);
    LatticeComplex one(UGrid); one = ComplexD(1.0,0.0);
    LatticeComplex tmp(UGrid);

    for (unsigned int i = 0; i < smearingLevels; ++i) {
      SmearedSet.push_back(*(new GaugeField(UGrid)));
      masks.push_back(*(new LatticeLorentzComplex(UGrid)));
      if (domask) {

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
	
      } else {
	for(int mu=0;mu<Nd;mu++){
	  PokeIndex<LorentzIndex>(masks[i],one, mu);
	}
      }
    }
    delete UrbGrid;
  }

  /*! For just thin links */
  //  SmearedConfigurationMasked()
  //    : smearingLevels(0), StoutSmearing(nullptr), SmearedSet(), ThinLinks(NULL), UGrid(NULL), UrbGrid(NULL), masks() {}

  // attach the smeared routines to the thin links U and fill the smeared set
  void set_Field(GaugeField &U)
  {
    double start = usecond();
    fill_smearedSet(U);
    double end = usecond();
    double time = (end - start)/ 1e3;
    std::cout << GridLogMessage << "Smearing in " << time << " ms" << std::endl;  
  }

  //====================================================================
  void smeared_force(GaugeField &SigmaTilde) 
  {
    if (smearingLevels > 0)
    {
      double start = usecond();
      GaugeField force = SigmaTilde; // actually = U*SigmaTilde
      GaugeLinkField tmp_mu(SigmaTilde.Grid());

      for (int mu = 0; mu < Nd; mu++)
      {
        // to get just SigmaTilde
        tmp_mu = adj(peekLorentz(SmearedSet[smearingLevels - 1], mu)) * peekLorentz(force, mu);
        pokeLorentz(force, tmp_mu, mu);
      }

      for (int ismr = smearingLevels - 1; ismr > 0; --ismr) {
        force = AnalyticSmearedForce(force, get_smeared_conf(ismr - 1),ismr);
      }
      
      force = AnalyticSmearedForce(force, *ThinLinks,0);

      for (int mu = 0; mu < Nd; mu++)
      {
        tmp_mu = peekLorentz(*ThinLinks, mu) * peekLorentz(force, mu);
        pokeLorentz(SigmaTilde, tmp_mu, mu);
      }
      double end = usecond();
      double time = (end - start)/ 1e3;
      std::cout << GridLogMessage << "Smearing force in " << time << " ms" << std::endl;  
    }  // if smearingLevels = 0 do nothing
  }
  //====================================================================

  GaugeField& get_SmearedU() { return SmearedSet[smearingLevels - 1]; }

  GaugeField &get_U(bool smeared = false)
  {
    // get the config, thin links by default
    if (smeared)
    {
      if (smearingLevels)
      {
        RealD impl_plaq =
	  WilsonLoops<Gimpl>::avgPlaquette(SmearedSet[smearingLevels - 1]);
        std::cout << GridLogDebug << "getting Usmr Plaq: " << impl_plaq
                  << std::endl;
        return get_SmearedU();
      }
      else
      {
        RealD impl_plaq = WilsonLoops<Gimpl>::avgPlaquette(*ThinLinks);
        std::cout << GridLogDebug << "getting Thin Plaq: " << impl_plaq
                  << std::endl;
        return *ThinLinks;
      }
    }
    else
    {
      RealD impl_plaq = WilsonLoops<Gimpl>::avgPlaquette(*ThinLinks);
      std::cout << GridLogDebug << "getting Thin Plaq: " << impl_plaq
                << std::endl;
      return *ThinLinks;
    }
  }
};

NAMESPACE_END(Grid);

