    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/ContinuedFractionFermion5D.cc

    Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
    void ContinuedFractionFermion5D<Impl>::SetCoefficientsTanh(Approx::zolotarev_data *zdata,RealD scale)
    {
      SetCoefficientsZolotarev(1.0/scale,zdata);
    }
    template<class Impl>
    void ContinuedFractionFermion5D<Impl>::SetCoefficientsZolotarev(RealD zolo_hi,Approx::zolotarev_data *zdata)
    {
      // How to check Ls matches??
      //      std::cout<<GridLogMessage << Ls << " Ls"<<std::endl;
      //      std::cout<<GridLogMessage << zdata->n  << " - n"<<std::endl;
      //      std::cout<<GridLogMessage << zdata->da << " -da "<<std::endl;
      //      std::cout<<GridLogMessage << zdata->db << " -db"<<std::endl;
      //      std::cout<<GridLogMessage << zdata->dn << " -dn"<<std::endl;
      //      std::cout<<GridLogMessage << zdata->dd << " -dd"<<std::endl;
      int Ls = this->Ls;
      assert(zdata->db==Ls);// Beta has Ls coeffs

      R=(1+this->mass)/(1-this->mass);

      Beta.resize(Ls);
      cc.resize(Ls);
      cc_d.resize(Ls);
      sqrt_cc.resize(Ls);
      for(int i=0; i < Ls ; i++){
	Beta[i] = zdata -> beta[i];
	cc[i] = 1.0/Beta[i];
	cc_d[i]=sqrt(cc[i]);
      }
    
      cc_d[Ls-1]=1.0;
      for(int i=0; i < Ls-1 ; i++){
	sqrt_cc[i]= sqrt(cc[i]*cc[i+1]);
      }    
      sqrt_cc[Ls-2]=sqrt(cc[Ls-2]);


      ZoloHiInv =1.0/zolo_hi;
      dw_diag = (4.0-this->M5)*ZoloHiInv;
    
      See.resize(Ls);
      Aee.resize(Ls);
      int sign=1;
      for(int s=0;s<Ls;s++){
	Aee[s] = sign * Beta[s] * dw_diag;
	sign   = - sign;
      }
      Aee[Ls-1] += R;
    
      See[0] = Aee[0];
      for(int s=1;s<Ls;s++){
	See[s] = Aee[s] - 1.0/See[s-1];
      }
      for(int s=0;s<Ls;s++){
	std::cout<<GridLogMessage <<"s = "<<s<<" Beta "<<Beta[s]<<" Aee "<<Aee[s] <<" See "<<See[s] <<std::endl;
      }
    }



    template<class Impl>
    RealD  ContinuedFractionFermion5D<Impl>::M           (const FermionField &psi, FermionField &chi)
    {
      int Ls = this->Ls;

      FermionField D(psi._grid);

      this->DW(psi,D,DaggerNo); 

      int sign=1;
      for(int s=0;s<Ls;s++){
	if ( s==0 ) {
	  ag5xpby_ssp(chi,cc[0]*Beta[0]*sign*ZoloHiInv,D,sqrt_cc[0],psi,s,s+1); // Multiplies Dw by G5 so Hw
	} else if ( s==(Ls-1) ){
	  RealD R=(1.0+mass)/(1.0-mass);
	  ag5xpby_ssp(chi,Beta[s]*ZoloHiInv,D,sqrt_cc[s-1],psi,s,s-1);
	  ag5xpby_ssp(chi,R,psi,1.0,chi,s,s);
	} else {
	  ag5xpby_ssp(chi,cc[s]*Beta[s]*sign*ZoloHiInv,D,sqrt_cc[s],psi,s,s+1);
  	  axpby_ssp(chi,1.0,chi,sqrt_cc[s-1],psi,s,s-1);
	}
	sign=-sign; 
      }
      return norm2(chi);
    }
    template<class Impl>
    RealD  ContinuedFractionFermion5D<Impl>::Mdag        (const FermionField &psi, FermionField &chi)
    {
      // This matrix is already hermitian. (g5 Dw) = Dw dag g5 = (g5 Dw)dag
      // The rest of matrix is symmetric.
      // Can ignore "dag"
      return M(psi,chi);
    }
    template<class Impl>
    void  ContinuedFractionFermion5D<Impl>::Mdir (const FermionField &psi, FermionField &chi,int dir,int disp){
      int Ls = this->Ls;

      this->DhopDir(psi,chi,dir,disp); // Dslash on diagonal. g5 Dslash is hermitian

      int sign=1;
      for(int s=0;s<Ls;s++){
	if ( s==(Ls-1) ){
	  ag5xpby_ssp(chi,Beta[s]*ZoloHiInv,chi,0.0,chi,s,s);
	} else {
	  ag5xpby_ssp(chi,cc[s]*Beta[s]*sign*ZoloHiInv,chi,0.0,chi,s,s);
	}
	sign=-sign; 
      }
    }
    template<class Impl>
    void   ContinuedFractionFermion5D<Impl>::Meooe       (const FermionField &psi, FermionField &chi)
    {
      int Ls = this->Ls;

      // Apply 4d dslash
      if ( psi.checkerboard == Odd ) {
	this->DhopEO(psi,chi,DaggerNo); // Dslash on diagonal. g5 Dslash is hermitian
      } else {
	this->DhopOE(psi,chi,DaggerNo); // Dslash on diagonal. g5 Dslash is hermitian
      }
      
      int sign=1;
      for(int s=0;s<Ls;s++){
	if ( s==(Ls-1) ){
	  ag5xpby_ssp(chi,Beta[s]*ZoloHiInv,chi,0.0,chi,s,s);
	} else {
	  ag5xpby_ssp(chi,cc[s]*Beta[s]*sign*ZoloHiInv,chi,0.0,chi,s,s);
	}
	sign=-sign; 
      }
    }
    template<class Impl>
    void   ContinuedFractionFermion5D<Impl>::MeooeDag    (const FermionField &psi, FermionField &chi)
    {
      this->Meooe(psi,chi);
    }
    template<class Impl>
    void   ContinuedFractionFermion5D<Impl>::Mooee       (const FermionField &psi, FermionField &chi)
    {
      int Ls = this->Ls;

      int sign=1;
      for(int s=0;s<Ls;s++){
	if ( s==0 ) {
	  ag5xpby_ssp(chi,cc[0]*Beta[0]*sign*dw_diag,psi,sqrt_cc[0],psi,s,s+1); // Multiplies Dw by G5 so Hw
	} else if ( s==(Ls-1) ){
	  // Drop the CC here.
	  double R=(1+mass)/(1-mass);
	  ag5xpby_ssp(chi,Beta[s]*dw_diag,psi,sqrt_cc[s-1],psi,s,s-1);
	  ag5xpby_ssp(chi,R,psi,1.0,chi,s,s);
	} else {
	  ag5xpby_ssp(chi,cc[s]*Beta[s]*sign*dw_diag,psi,sqrt_cc[s],psi,s,s+1);
	  axpby_ssp(chi,1.0,chi,sqrt_cc[s-1],psi,s,s-1);
	}
	sign=-sign; 
      }
    }

    template<class Impl>
    void   ContinuedFractionFermion5D<Impl>::MooeeDag    (const FermionField &psi, FermionField &chi)
    {
      this->Mooee(psi,chi);
    }
    template<class Impl>
    void   ContinuedFractionFermion5D<Impl>::MooeeInv    (const FermionField &psi, FermionField &chi)
    {
      int Ls = this->Ls;

      // Apply Linv
      axpby_ssp(chi,1.0/cc_d[0],psi,0.0,psi,0,0); 
      for(int s=1;s<Ls;s++){
	axpbg5y_ssp(chi,1.0/cc_d[s],psi,-1.0/See[s-1],chi,s,s-1);
      }
      // Apply Dinv
      for(int s=0;s<Ls;s++){
	ag5xpby_ssp(chi,1.0/See[s],chi,0.0,chi,s,s); //only appearance of See[0]
      }
      // Apply Uinv = (Linv)^T
      axpby_ssp(chi,1.0/cc_d[Ls-1],chi,0.0,chi,Ls-1,Ls-1);
      for(int s=Ls-2;s>=0;s--){
	axpbg5y_ssp(chi,1.0/cc_d[s],chi,-1.0*cc_d[s+1]/See[s]/cc_d[s],chi,s,s+1);
      }
    }
    template<class Impl>
    void   ContinuedFractionFermion5D<Impl>::MooeeInvDag (const FermionField &psi, FermionField &chi)
    {
      this->MooeeInv(psi,chi);
    }

  // force terms; five routines; default to Dhop on diagonal
    template<class Impl>
   void ContinuedFractionFermion5D<Impl>::MDeriv  (GaugeField &mat,const FermionField &U,const FermionField &V,int dag)
  {
    int Ls = this->Ls;

    FermionField D(V._grid);

    int sign=1;
    for(int s=0;s<Ls;s++){
      if ( s==(Ls-1) ){
	ag5xpby_ssp(D,Beta[s]*ZoloHiInv,U,0.0,U,s,s);
      } else {
	ag5xpby_ssp(D,cc[s]*Beta[s]*sign*ZoloHiInv,U,0.0,U,s,s);
      }
      sign=-sign; 
    }
    this->DhopDeriv(mat,D,V,DaggerNo); 
  };
    template<class Impl>
   void ContinuedFractionFermion5D<Impl>::MoeDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag)
  {
    int Ls = this->Ls;

    FermionField D(V._grid);

    int sign=1;
    for(int s=0;s<Ls;s++){
      if ( s==(Ls-1) ){
	ag5xpby_ssp(D,Beta[s]*ZoloHiInv,U,0.0,U,s,s);
      } else {
	ag5xpby_ssp(D,cc[s]*Beta[s]*sign*ZoloHiInv,U,0.0,U,s,s);
      }
      sign=-sign; 
    }
    this->DhopDerivOE(mat,D,V,DaggerNo); 
  };
  template<class Impl>
  void ContinuedFractionFermion5D<Impl>::MeoDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag)
  {
    int Ls = this->Ls;

    FermionField D(V._grid);

    int sign=1;
    for(int s=0;s<Ls;s++){
      if ( s==(Ls-1) ){
	ag5xpby_ssp(D,Beta[s]*ZoloHiInv,U,0.0,U,s,s);
      } else {
	ag5xpby_ssp(D,cc[s]*Beta[s]*sign*ZoloHiInv,U,0.0,U,s,s);
      }
      sign=-sign; 
    }
    this->DhopDerivEO(mat,D,V,DaggerNo); 
  };
    
    // Constructors
    template<class Impl>
    ContinuedFractionFermion5D<Impl>::ContinuedFractionFermion5D(
							   GaugeField &_Umu,
							   GridCartesian         &FiveDimGrid,
							   GridRedBlackCartesian &FiveDimRedBlackGrid,
							   GridCartesian         &FourDimGrid,
							   GridRedBlackCartesian &FourDimRedBlackGrid,
							   RealD _mass,RealD M5,const ImplParams &p) :
      WilsonFermion5D<Impl>(_Umu,
			    FiveDimGrid, FiveDimRedBlackGrid,
			    FourDimGrid, FourDimRedBlackGrid,M5,p),
      mass(_mass)
    {
      int Ls = this->Ls;
      assert((Ls&0x1)==1); // Odd Ls required
    }

    FermOpTemplateInstantiate(ContinuedFractionFermion5D);

  }
}

