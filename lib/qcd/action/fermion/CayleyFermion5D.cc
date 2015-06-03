#include <Grid.h>
namespace Grid {
namespace QCD {

 CayleyFermion5D::CayleyFermion5D(LatticeGaugeField &_Umu,
				  GridCartesian         &FiveDimGrid,
				  GridRedBlackCartesian &FiveDimRedBlackGrid,
				  GridCartesian         &FourDimGrid,
				  GridRedBlackCartesian &FourDimRedBlackGrid,
				  RealD _mass,RealD _M5) :
   WilsonFermion5D(_Umu,
		   FiveDimGrid,
		   FiveDimRedBlackGrid,
		   FourDimGrid,
		   FourDimRedBlackGrid,_M5),
   mass(_mass)
 {
 }

  // override multiply
  RealD CayleyFermion5D::M    (const LatticeFermion &psi, LatticeFermion &chi)
  {
    LatticeFermion Din(psi._grid);

    // Assemble Din
    for(int s=0;s<Ls;s++){
      if ( s==0 ) {
	//	Din = bs psi[s] + cs[s] psi[s+1}
	axpby_ssp_pminus(Din,bs[s],psi,cs[s],psi,s,s+1);
	//      Din+= -mass*cs[s] psi[s+1}
	axpby_ssp_pplus (Din,1.0,Din,-mass*cs[s],psi,s,Ls-1);
      } else if ( s==(Ls-1)) { 
	axpby_ssp_pminus(Din,bs[s],psi,-mass*cs[s],psi,s,0);
	axpby_ssp_pplus (Din,1.0,Din,cs[s],psi,s,s-1);
      } else {
	axpby_ssp_pminus(Din,bs[s],psi,cs[s],psi,s,s+1);
	axpby_ssp_pplus(Din,1.0,Din,cs[s],psi,s,s-1);
      }
    }

    DW(Din,chi,DaggerNo);
    // ((b D_W + D_w hop terms +1) on s-diag
    axpby(chi,1.0,1.0,chi,psi); 

    for(int s=0;s<Ls;s++){
      if ( s==0 ){
	axpby_ssp_pminus(chi,1.0,chi,-1.0,psi,s,s+1);
	axpby_ssp_pplus (chi,1.0,chi,mass,psi,s,Ls-1);
      } else if ( s==(Ls-1)) {
	axpby_ssp_pminus(chi,1.0,chi,mass,psi,s,0);
	axpby_ssp_pplus (chi,1.0,chi,-1.0,psi,s,s-1);
      } else {
	axpby_ssp_pminus(chi,1.0,chi,-1.0,psi,s,s+1);
	axpby_ssp_pplus (chi,1.0,chi,-1.0,psi,s,s-1);
      }
    }
    return norm2(chi);
  }

  RealD CayleyFermion5D::Mdag (const LatticeFermion &psi, LatticeFermion &chi)
  {
    // Under adjoint
    //D1+        D1- P-    ->   D1+^dag   P+ D2-^dag
    //D2- P+     D2+            P-D1-^dag D2+dag

    LatticeFermion Din(psi._grid);
    // Apply Dw
    DW(psi,Din,DaggerYes); 

    for(int s=0;s<Ls;s++){
      // Collect the terms in DW
      //	Chi = bs Din[s] + cs[s] Din[s+1}
      //    Chi+= -mass*cs[s] psi[s+1}
      if ( s==0 ) {
	axpby_ssp_pplus (chi,bs[s],Din,cs[s+1],Din,s,s+1);
	axpby_ssp_pminus(chi,1.0,chi,-mass*cs[Ls-1],Din,s,Ls-1);
      } else if ( s==(Ls-1)) { 
	axpby_ssp_pplus (chi,bs[s],Din,-mass*cs[0],Din,s,0);
	axpby_ssp_pminus(chi,1.0,chi,cs[s-1],Din,s,s-1);
      } else {
	axpby_ssp_pplus (chi,bs[s],Din,cs[s+1],Din,s,s+1);
	axpby_ssp_pminus(chi,1.0,chi,cs[s-1],Din,s,s-1);
      }
      // Collect the terms indept of DW
      if ( s==0 ){
	axpby_ssp_pplus (chi,1.0,chi,-1.0,psi,s,s+1);
	axpby_ssp_pminus(chi,1.0,chi,mass,psi,s,Ls-1);
      } else if ( s==(Ls-1)) {
	axpby_ssp_pplus (chi,1.0,chi,mass,psi,s,0);
	axpby_ssp_pminus(chi,1.0,chi,-1.0,psi,s,s-1);
      } else {
	axpby_ssp_pplus(chi,1.0,chi,-1.0,psi,s,s+1);
	axpby_ssp_pminus(chi,1.0,chi,-1.0,psi,s,s-1);
      }
    }
    // ((b D_W + D_w hop terms +1) on s-diag
    axpby (chi,1.0,1.0,chi,psi); 
    return norm2(chi);
  }

  // half checkerboard operations
  void CayleyFermion5D::Meooe       (const LatticeFermion &psi, LatticeFermion &chi)
  {
    LatticeFermion tmp(psi._grid);
    // Assemble the 5d matrix
    for(int s=0;s<Ls;s++){
      if ( s==0 ) {
	//	tmp = bs psi[s] + cs[s] psi[s+1}
	//      tmp+= -mass*cs[s] psi[s+1}
	axpby_ssp_pminus(tmp,beo[s],psi,-ceo[s],psi ,s, s+1);
	axpby_ssp_pplus(tmp,1.0,tmp,mass*ceo[s],psi,s,Ls-1);
      } else if ( s==(Ls-1)) { 
	axpby_ssp_pminus(tmp,beo[s],psi,mass*ceo[s],psi,s,0);
	axpby_ssp_pplus(tmp,1.0,tmp,-ceo[s],psi,s,s-1);
      } else {
	axpby_ssp_pminus(tmp,beo[s],psi,-ceo[s],psi,s,s+1);
	axpby_ssp_pplus (tmp,1.0,tmp,-ceo[s],psi,s,s-1);
      }
    }
    // Apply 4d dslash
    if ( psi.checkerboard == Odd ) {
      DhopEO(tmp,chi,DaggerNo);
    } else {
      DhopOE(tmp,chi,DaggerNo);
    }
  }

  void CayleyFermion5D::MeooeDag    (const LatticeFermion &psi, LatticeFermion &chi)
  {
    LatticeFermion tmp(psi._grid);
    // Apply 4d dslash
    if ( psi.checkerboard == Odd ) {
      DhopEO(psi,tmp,DaggerYes);
    } else {
      DhopOE(psi,tmp,DaggerYes);
    }
    // Assemble the 5d matrix
    for(int s=0;s<Ls;s++){
      if ( s==0 ) {
	axpby_ssp_pplus(chi,beo[s],tmp,   -ceo[s+1]  ,tmp,s,s+1);
	axpby_ssp_pminus(chi,   1.0,chi,mass*ceo[Ls-1],tmp,s,Ls-1);
      } else if ( s==(Ls-1)) { 
	axpby_ssp_pplus(chi,beo[s],tmp,mass*ceo[0],tmp,s,0);
	axpby_ssp_pminus(chi,1.0,chi,-ceo[s-1],tmp,s,s-1);
      } else {
	axpby_ssp_pplus(chi,beo[s],tmp,-ceo[s+1],tmp,s,s+1);
	axpby_ssp_pminus(chi,1.0   ,chi,-ceo[s-1],tmp,s,s-1);
      }
    }
  }

  void CayleyFermion5D::Mooee       (const LatticeFermion &psi, LatticeFermion &chi)
  {
    for (int s=0;s<Ls;s++){
      if ( s==0 ) {
	axpby_ssp_pminus(chi,bee[s],psi ,-cee[s],psi,s,s+1);
	axpby_ssp_pplus (chi,1.0,chi,mass*cee[s],psi,s,Ls-1);
      } else if ( s==(Ls-1)) { 
	axpby_ssp_pminus(chi,bee[s],psi,mass*cee[s],psi,s,0);
	axpby_ssp_pplus (chi,1.0,chi,-cee[s],psi,s,s-1);
      } else {
	axpby_ssp_pminus(chi,bee[s],psi,-cee[s],psi,s,s+1);
	axpby_ssp_pplus (chi,1.0,chi,-cee[s],psi,s,s-1);
      }
    }
  }

  void CayleyFermion5D::MooeeDag    (const LatticeFermion &psi, LatticeFermion &chi)
  {
    for (int s=0;s<Ls;s++){
      // Assemble the 5d matrix
      if ( s==0 ) {
	axpby_ssp_pplus(chi,bee[s],psi,-cee[s+1]  ,psi,s,s+1);
	axpby_ssp_pminus(chi,1.0,chi,mass*cee[Ls-1],psi,s,Ls-1);
      } else if ( s==(Ls-1)) { 
	axpby_ssp_pplus(chi,bee[s],psi,mass*cee[0],psi,s,0);
	axpby_ssp_pminus(chi,1.0,chi,-cee[s-1],psi,s,s-1);
      } else {
	axpby_ssp_pplus(chi,bee[s],psi,-cee[s+1],psi,s,s+1);
	axpby_ssp_pminus(chi,1.0   ,chi,-cee[s-1],psi,s,s-1);
      }
    }
  }

  void CayleyFermion5D::MooeeInv    (const LatticeFermion &psi, LatticeFermion &chi)
  {
    // Apply (L^{\prime})^{-1}
    axpby_ssp (chi,1.0,psi,     0.0,psi,0,0);      // chi[0]=psi[0]
    for (int s=1;s<Ls;s++){
      axpby_ssp_pplus(chi,1.0,psi,-lee[s-1],chi,s,s-1);// recursion Psi[s] -lee P_+ chi[s-1]
    }
    // L_m^{-1} 
    for (int s=0;s<Ls-1;s++){ // Chi[ee] = 1 - sum[s<Ls-1] -leem[s]P_- chi
      axpby_ssp_pminus(chi,1.0,chi,-leem[s],chi,Ls-1,s);
    }
    // U_m^{-1} D^{-1}
    for (int s=0;s<Ls-1;s++){
      // Chi[s] + 1/d chi[s] 
      axpby_ssp_pplus(chi,1.0/dee[s],chi,-ueem[s]/dee[Ls-1],chi,s,Ls-1);
    }	
    axpby_ssp(chi,1.0/dee[Ls-1],chi,0.0,chi,Ls-1,Ls-1); // Modest avoidable 
    
    // Apply U^{-1}
    for (int s=Ls-2;s>=0;s--){
      axpby_ssp_pminus (chi,1.0,chi,-uee[s],chi,s,s+1);  // chi[Ls]
    }
  }

  void CayleyFermion5D::MooeeInvDag (const LatticeFermion &psi, LatticeFermion &chi)
  {
    // Apply (U^{\prime})^{-dagger}
    axpby_ssp (chi,1.0,psi,     0.0,psi,0,0);      // chi[0]=psi[0]
    for (int s=1;s<Ls;s++){
      axpby_ssp_pminus(chi,1.0,psi,-uee[s-1],chi,s,s-1);
    }
    // U_m^{-\dagger} 
    for (int s=0;s<Ls-1;s++){
      axpby_ssp_pplus(chi,1.0,chi,-ueem[s],chi,Ls-1,s);
    }
    // L_m^{-\dagger} D^{-dagger}
    for (int s=0;s<Ls-1;s++){
      axpby_ssp_pminus(chi,1.0/dee[s],chi,-leem[s]/dee[Ls-1],chi,s,Ls-1);
    }	
    axpby_ssp(chi,1.0/dee[Ls-1],chi,0.0,chi,Ls-1,Ls-1); // Modest avoidable 
    
    // Apply L^{-dagger}
    for (int s=Ls-2;s>=0;s--){
      axpby_ssp_pplus (chi,1.0,chi,-lee[s],chi,s,s+1);  // chi[Ls]
    }
  }
  
  // Tanh
  void CayleyFermion5D::SetCoefficientsTanh(Approx::zolotarev_data *zdata,RealD b,RealD c)
  {
    SetCoefficientsZolotarev(1.0,zdata,b,c);

  }
  //Zolo
  void CayleyFermion5D::SetCoefficientsZolotarev(RealD zolo_hi,Approx::zolotarev_data *zdata,RealD b,RealD c)
  {

    ///////////////////////////////////////////////////////////
    // The Cayley coeffs (unprec)
    ///////////////////////////////////////////////////////////
    omega.resize(Ls);
    bs.resize(Ls);
    cs.resize(Ls);
    as.resize(Ls);
    
    // 
    // Ts = (    [bs+cs]Dw        )^-1 (    (bs+cs) Dw         )
    //     -(g5  -------       -1 )    ( g5 ---------     + 1  )
    //      (   {2+(bs-cs)Dw}     )    (    2+(bs-cs) Dw       )
    //
    //  bs = 1/2( (1/omega_s + 1)*b + (1/omega - 1)*c ) = 1/2(  1/omega(b+c) + (b-c) )
    //  cs = 1/2( (1/omega_s - 1)*b + (1/omega + 1)*c ) = 1/2(  1/omega(b+c) - (b-c) )
    //
    // bs+cs = 0.5*( 1/omega(b+c) + (b-c) + 1/omega(b+c) - (b-c) ) = 1/omega(b+c)
    // bs-cs = 0.5*( 1/omega(b+c) + (b-c) - 1/omega(b+c) + (b-c) ) = b-c
    //
    // So 
    //
    // Ts = (    [b+c]Dw/omega_s    )^-1 (    (b+c) Dw /omega_s        )
    //     -(g5  -------         -1 )    ( g5 ---------           + 1  )
    //      (   {2+(b-c)Dw}         )    (    2+(b-c) Dw               )
    //
    // Ts = (    [b+c]Dw            )^-1 (    (b+c) Dw                 )
    //     -(g5  -------    -omega_s)    ( g5 ---------      + omega_s )
    //      (   {2+(b-c)Dw}         )    (    2+(b-c) Dw               )
    // 
    
    double bpc = b+c;
    double bmc = b-c;
    for(int i=0; i < Ls; i++){
      as[i] = 1.0;
      omega[i] = ((double)zdata->gamma[i])*zolo_hi; //NB reciprocal relative to Chroma NEF code
      bs[i] = 0.5*(bpc/omega[i] + bmc);
      cs[i] = 0.5*(bpc/omega[i] - bmc);
    }

    ////////////////////////////////////////////////////////
    // Constants for the preconditioned matrix Cayley form
    ////////////////////////////////////////////////////////
    bee.resize(Ls);
    cee.resize(Ls);
    beo.resize(Ls);
    ceo.resize(Ls);
    
    for(int i=0;i<Ls;i++){
      bee[i]=as[i]*(bs[i]*(4.0-M5) +1.0);
      cee[i]=as[i]*(1.0-cs[i]*(4.0-M5));
      beo[i]=as[i]*bs[i];
      ceo[i]=-as[i]*cs[i];
    }

    aee.resize(Ls);
    aeo.resize(Ls);
    for(int i=0;i<Ls;i++){
      aee[i]=cee[i];
      aeo[i]=ceo[i];
    }

    //////////////////////////////////////////
    // LDU decomposition of eeoo
    //////////////////////////////////////////
    dee.resize(Ls);
    lee.resize(Ls);
    leem.resize(Ls);
    uee.resize(Ls);
    ueem.resize(Ls);
    
    for(int i=0;i<Ls;i++){
      
      dee[i] = bee[i];
      
      if ( i < Ls-1 ) {
	
	lee[i] =-cee[i+1]/bee[i]; // sub-diag entry on the ith column
	    
	leem[i]=mass*cee[Ls-1]/bee[0];
	for(int j=0;j<i;j++)  leem[i]*= aee[j]/bee[j+1];
	
	uee[i] =-aee[i]/bee[i];   // up-diag entry on the ith row
	
	ueem[i]=mass;
	for(int j=1;j<=i;j++) ueem[i]*= cee[j]/bee[j];
	ueem[i]*= aee[0]/bee[0];
	    
      } else { 
	lee[i] =0.0;
	leem[i]=0.0;
	uee[i] =0.0;
	ueem[i]=0.0;
      }
    }
	
    { 
      double delta_d=mass*cee[Ls-1];
      for(int j=0;j<Ls-1;j++) delta_d *= cee[j]/bee[j];
      dee[Ls-1] += delta_d;
    }
  }

}}


