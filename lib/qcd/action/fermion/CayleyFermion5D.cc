    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/CayleyFermion5D.cc

    Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: paboyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#include <Grid.h>
namespace Grid {
namespace QCD {

 template<class Impl>
 CayleyFermion5D<Impl>::CayleyFermion5D(GaugeField &_Umu,
					GridCartesian         &FiveDimGrid,
					GridRedBlackCartesian &FiveDimRedBlackGrid,
					GridCartesian         &FourDimGrid,
					GridRedBlackCartesian &FourDimRedBlackGrid,
					RealD _mass,RealD _M5,const ImplParams &p) :
   WilsonFermion5D<Impl>(_Umu,
		   FiveDimGrid,
		   FiveDimRedBlackGrid,
		   FourDimGrid,
 	 	   FourDimRedBlackGrid,_M5,p),
   mass(_mass)
 {
 }

 template<class Impl>
  void CayleyFermion5D<Impl>::Meooe5D    (const FermionField &psi, FermionField &Din)
  {
    // Assemble Din
    int Ls=this->Ls;
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
  }
 template<class Impl>
  void CayleyFermion5D<Impl>::MeooeDag5D    (const FermionField &psi, FermionField &Din)
  {
    int Ls=this->Ls;
    for(int s=0;s<Ls;s++){
      if ( s==0 ) {
	axpby_ssp_pplus (Din,bs[s],psi,cs[s+1],psi,s,s+1);
	axpby_ssp_pminus(Din,1.0,Din,-mass*cs[Ls-1],psi,s,Ls-1);
      } else if ( s==(Ls-1)) { 
	axpby_ssp_pplus (Din,bs[s],psi,-mass*cs[0],psi,s,0);
	axpby_ssp_pminus(Din,1.0,Din,cs[s-1],psi,s,s-1);
      } else {
	axpby_ssp_pplus (Din,bs[s],psi,cs[s+1],psi,s,s+1);
	axpby_ssp_pminus(Din,1.0,Din,cs[s-1],psi,s,s-1);
      }
    }
  }

  // override multiply
 template<class Impl>
  RealD CayleyFermion5D<Impl>::M    (const FermionField &psi, FermionField &chi)
  {
    int Ls=this->Ls;

    FermionField Din(psi._grid);

    // Assemble Din
    /*
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
    */
    Meooe5D(psi,Din);

    this->DW(Din,chi,DaggerNo);
    // ((b D_W + D_w hop terms +1) on s-diag
    axpby(chi,1.0,1.0,chi,psi); 

    // Call Mooee??
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

 template<class Impl>
  RealD CayleyFermion5D<Impl>::Mdag (const FermionField &psi, FermionField &chi)
  {
    // Under adjoint
    //D1+        D1- P-    ->   D1+^dag   P+ D2-^dag
    //D2- P+     D2+            P-D1-^dag D2+dag

    FermionField Din(psi._grid);
    // Apply Dw
    this->DW(psi,Din,DaggerYes); 

    MeooeDag5D(Din,chi);

    int Ls=this->Ls;
    for(int s=0;s<Ls;s++){

      // Collect the terms in DW
      //	Chi = bs Din[s] + cs[s] Din[s+1}
      //    Chi+= -mass*cs[s] psi[s+1}
      /*
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
      */

      // FIXME just call MooeeDag??

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
 template<class Impl>
  void CayleyFermion5D<Impl>::Meooe       (const FermionField &psi, FermionField &chi)
  {
    int Ls=this->Ls;

    FermionField tmp(psi._grid);
    // Assemble the 5d matrix
    Meooe5D(psi,tmp); 
#if 0
    std::cout << "Meooe Test replacement norm2 tmp = " <<norm2(tmp)<<std::endl;
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
    std::cout << "Meooe Test replacement norm2 tmp old = " <<norm2(tmp)<<std::endl;
#endif

    // Apply 4d dslash
    if ( psi.checkerboard == Odd ) {
      this->DhopEO(tmp,chi,DaggerNo);
    } else {
      this->DhopOE(tmp,chi,DaggerNo);
    }
  }

  template<class Impl>
  void CayleyFermion5D<Impl>::MeooeDag    (const FermionField &psi, FermionField &chi)
  {
    FermionField tmp(psi._grid);
    // Apply 4d dslash
    if ( psi.checkerboard == Odd ) {
      this->DhopEO(psi,tmp,DaggerYes);
    } else {
      this->DhopOE(psi,tmp,DaggerYes);
    }

    MeooeDag5D(tmp,chi); 
#if 0
    std::cout << "Meooe Test replacement norm2 chi new = " <<norm2(chi)<<std::endl;
    // Assemble the 5d matrix
    int Ls=this->Ls;
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
    std::cout << "Meooe Test replacement norm2 chi old = " <<norm2(chi)<<std::endl;
#endif

  }

 template<class Impl>
  void CayleyFermion5D<Impl>::Mooee       (const FermionField &psi, FermionField &chi)
  {
    int Ls=this->Ls;
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

 template<class Impl>
  void  CayleyFermion5D<Impl>::Mdir (const FermionField &psi, FermionField &chi,int dir,int disp){
    int Ls=this->Ls;
    FermionField tmp(psi._grid);
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
    // Apply 4d dslash fragment
    this->DhopDir(tmp,chi,dir,disp);
  }

 template<class Impl>
  void CayleyFermion5D<Impl>::MooeeDag    (const FermionField &psi, FermionField &chi)
  {
    int Ls=this->Ls;
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

 template<class Impl>
  void CayleyFermion5D<Impl>::MooeeInv    (const FermionField &psi, FermionField &chi)
  {
    int Ls=this->Ls;
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

 template<class Impl>
  void CayleyFermion5D<Impl>::MooeeInvDag (const FermionField &psi, FermionField &chi)
  {
    int Ls=this->Ls;
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

  // force terms; five routines; default to Dhop on diagonal
  template<class Impl>
  void CayleyFermion5D<Impl>::MDeriv  (GaugeField &mat,const FermionField &U,const FermionField &V,int dag)
  {
    FermionField Din(V._grid);

    if ( dag == DaggerNo ) {
      //      U d/du [D_w D5] V = U d/du DW D5 V
      Meooe5D(V,Din);
      this->DhopDeriv(mat,U,Din,dag);
    } else {
      //      U d/du [D_w D5]^dag V = U D5^dag d/du DW^dag Y // implicit adj on U in call
      Meooe5D(U,Din);
      this->DhopDeriv(mat,Din,V,dag);
    }
  };
 template<class Impl>
  void CayleyFermion5D<Impl>::MoeDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag)
  {
    FermionField Din(V._grid);

    if ( dag == DaggerNo ) {
      //      U d/du [D_w D5] V = U d/du DW D5 V
      Meooe5D(V,Din);
      this->DhopDerivOE(mat,U,Din,dag);
    } else {
      //      U d/du [D_w D5]^dag V = U D5^dag d/du DW^dag Y // implicit adj on U in call
      Meooe5D(U,Din);
      this->DhopDerivOE(mat,Din,V,dag);
    }
  };
 template<class Impl>
  void CayleyFermion5D<Impl>::MeoDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag)
  {
    FermionField Din(V._grid);

    if ( dag == DaggerNo ) {
      //      U d/du [D_w D5] V = U d/du DW D5 V
      Meooe5D(V,Din);
      this->DhopDerivEO(mat,U,Din,dag);
    } else {
      //      U d/du [D_w D5]^dag V = U D5^dag d/du DW^dag Y // implicit adj on U in call
      Meooe5D(U,Din);
      this->DhopDerivEO(mat,Din,V,dag);
    }
  };
  
  // Tanh
 template<class Impl>
  void CayleyFermion5D<Impl>::SetCoefficientsTanh(Approx::zolotarev_data *zdata,RealD b,RealD c)
  {
    SetCoefficientsZolotarev(1.0,zdata,b,c);

  }
  //Zolo
 template<class Impl>
  void CayleyFermion5D<Impl>::SetCoefficientsZolotarev(RealD zolo_hi,Approx::zolotarev_data *zdata,RealD b,RealD c)
  {
    int Ls=this->Ls;

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
      bee[i]=as[i]*(bs[i]*(4.0-this->M5) +1.0);
      cee[i]=as[i]*(1.0-cs[i]*(4.0-this->M5));
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

  FermOpTemplateInstantiate(CayleyFermion5D);

}}


