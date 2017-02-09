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

#include <Grid/Eigen/Dense>
#include <Grid/Grid.h>


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
void CayleyFermion5D<Impl>::Dminus(const FermionField &psi, FermionField &chi)
{
  int Ls=this->Ls;

  this->DW(psi,this->tmp(),DaggerNo);

  for(int s=0;s<Ls;s++){
    axpby_ssp(chi,Coeff_t(1.0),psi,-cs[s],this->tmp(),s,s);// chi = (1-c[s] D_W) psi
  }
}


template<class Impl> void CayleyFermion5D<Impl>::CayleyReport(void)
{
  this->Report();
  std::vector<int> latt = GridDefaultLatt();          
  RealD volume = this->Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt[mu];
  RealD NP     = this->_FourDimGrid->_Nprocessors;
  if ( M5Dcalls > 0 ) {
    std::cout << GridLogMessage << "#### M5D calls report " << std::endl;
    std::cout << GridLogMessage << "CayleyFermion5D Number of M5D Calls     : " << M5Dcalls   << std::endl;
    std::cout << GridLogMessage << "CayleyFermion5D ComputeTime/Calls       : " << M5Dtime / M5Dcalls << " us" << std::endl;

    // Flops = 6.0*(Nc*Ns) *Ls*vol
    RealD mflops = 6.0*12*volume*M5Dcalls/M5Dtime/2; // 2 for red black counting
    std::cout << GridLogMessage << "Average mflops/s per call                : " << mflops << std::endl;
    std::cout << GridLogMessage << "Average mflops/s per call per rank       : " << mflops/NP << std::endl;
  }

  if ( MooeeInvCalls > 0 ) {

    std::cout << GridLogMessage << "#### MooeeInv calls report " << std::endl;
    std::cout << GridLogMessage << "CayleyFermion5D Number of MooeeInv Calls     : " << MooeeInvCalls   << std::endl;
    std::cout << GridLogMessage << "CayleyFermion5D ComputeTime/Calls            : " << MooeeInvTime / MooeeInvCalls << " us" << std::endl;

    // Flops = MADD * Ls *Ls *4dvol * spin/colour/complex
    RealD mflops = 2.0*24*this->Ls*volume*MooeeInvCalls/MooeeInvTime/2; // 2 for red black counting
    std::cout << GridLogMessage << "Average mflops/s per call                : " << mflops << std::endl;
    std::cout << GridLogMessage << "Average mflops/s per call per rank       : " << mflops/NP << std::endl;
  }

}
template<class Impl> void CayleyFermion5D<Impl>::CayleyZeroCounters(void)
{
  this->ZeroCounters();
  M5Dflops=0;
  M5Dcalls=0;
  M5Dtime=0;
  MooeeInvFlops=0;
  MooeeInvCalls=0;
  MooeeInvTime=0;
}


template<class Impl>  
void CayleyFermion5D<Impl>::DminusDag(const FermionField &psi, FermionField &chi)
{
  int Ls=this->Ls;

  this->DW(psi,this->tmp(),DaggerYes);

  for(int s=0;s<Ls;s++){
    axpby_ssp(chi,Coeff_t(1.0),psi,-cs[s],this->tmp(),s,s);// chi = (1-c[s] D_W) psi
  }
}
template<class Impl>  
void CayleyFermion5D<Impl>::M5D   (const FermionField &psi, FermionField &chi)
{
  int Ls=this->Ls;
  std::vector<Coeff_t> diag (Ls,1.0);
  std::vector<Coeff_t> upper(Ls,-1.0); upper[Ls-1]=mass;
  std::vector<Coeff_t> lower(Ls,-1.0); lower[0]   =mass;
  M5D(psi,chi,chi,lower,diag,upper);
}
template<class Impl>
void CayleyFermion5D<Impl>::Meooe5D    (const FermionField &psi, FermionField &Din)
{
  int Ls=this->Ls;
  std::vector<Coeff_t> diag = bs;
  std::vector<Coeff_t> upper= cs;
  std::vector<Coeff_t> lower= cs; 
  upper[Ls-1]=-mass*upper[Ls-1];
  lower[0]   =-mass*lower[0];
  M5D(psi,psi,Din,lower,diag,upper);
}
// FIXME Redunant with the above routine; check this and eliminate
template<class Impl> void CayleyFermion5D<Impl>::Meo5D     (const FermionField &psi, FermionField &chi)
{
  int Ls=this->Ls;
  std::vector<Coeff_t> diag = beo;
  std::vector<Coeff_t> upper(Ls);
  std::vector<Coeff_t> lower(Ls);
  for(int i=0;i<Ls;i++) {
    upper[i]=-ceo[i];
    lower[i]=-ceo[i];
  }
  upper[Ls-1]=-mass*upper[Ls-1];
  lower[0]   =-mass*lower[0];
  M5D(psi,psi,chi,lower,diag,upper);
}
template<class Impl>
void CayleyFermion5D<Impl>::Mooee       (const FermionField &psi, FermionField &chi)
{
  int Ls=this->Ls;
  std::vector<Coeff_t> diag = bee;
  std::vector<Coeff_t> upper(Ls);
  std::vector<Coeff_t> lower(Ls);
  for(int i=0;i<Ls;i++) {
    upper[i]=-cee[i];
    lower[i]=-cee[i];
  }
  upper[Ls-1]=-mass*upper[Ls-1];
  lower[0]   =-mass*lower[0];
  M5D(psi,psi,chi,lower,diag,upper);
}

template<class Impl>
void CayleyFermion5D<Impl>::MooeeDag    (const FermionField &psi, FermionField &chi)
{
  int Ls=this->Ls;
  std::vector<Coeff_t> diag = bee;
  std::vector<Coeff_t> upper(Ls);
  std::vector<Coeff_t> lower(Ls);

  for (int s=0;s<Ls;s++){
    // Assemble the 5d matrix
    if ( s==0 ) {
      upper[s] = -cee[s+1] ;
      lower[s] = mass*cee[Ls-1];
    } else if ( s==(Ls-1)) { 
      upper[s] = mass*cee[0];
      lower[s] = -cee[s-1];
    } else {
      upper[s]=-cee[s+1];
      lower[s]=-cee[s-1];
    }
  }

  M5Ddag(psi,psi,chi,lower,diag,upper);
}

template<class Impl>
void CayleyFermion5D<Impl>::M5Ddag (const FermionField &psi, FermionField &chi)
{
  int Ls=this->Ls;
  std::vector<Coeff_t> diag(Ls,1.0);
  std::vector<Coeff_t> upper(Ls,-1.0);
  std::vector<Coeff_t> lower(Ls,-1.0);
  upper[Ls-1]=-mass*upper[Ls-1];
  lower[0]   =-mass*lower[0];
  M5Ddag(psi,chi,chi,lower,diag,upper);
}

template<class Impl>
void CayleyFermion5D<Impl>::MeooeDag5D    (const FermionField &psi, FermionField &Din)
{
  int Ls=this->Ls;
  std::vector<Coeff_t> diag =bs;
  std::vector<Coeff_t> upper=cs;
  std::vector<Coeff_t> lower=cs;
  upper[Ls-1]=-mass*upper[Ls-1];
  lower[0]   =-mass*lower[0];
  M5Ddag(psi,psi,Din,lower,diag,upper);
}

template<class Impl>
RealD CayleyFermion5D<Impl>::M    (const FermionField &psi, FermionField &chi)
{
  int Ls=this->Ls;
  
  FermionField Din(psi._grid);
  
  // Assemble Din
  Meooe5D(psi,Din);
  
  this->DW(Din,chi,DaggerNo);
  // ((b D_W + D_w hop terms +1) on s-diag
  axpby(chi,1.0,1.0,chi,psi); 
  
  M5D(psi,chi);
  return(norm2(chi));
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
  
  M5Ddag(psi,chi);
  // ((b D_W + D_w hop terms +1) on s-diag
  axpby (chi,1.0,1.0,chi,psi); 
  return norm2(chi);
}

// half checkerboard operations
template<class Impl>
void CayleyFermion5D<Impl>::Meooe       (const FermionField &psi, FermionField &chi)
{
  int Ls=this->Ls;

  Meooe5D(psi,this->tmp()); 

  if ( psi.checkerboard == Odd ) {
    this->DhopEO(this->tmp(),chi,DaggerNo);
  } else {
    this->DhopOE(this->tmp(),chi,DaggerNo);
  }
}

template<class Impl>
void CayleyFermion5D<Impl>::MeooeDag    (const FermionField &psi, FermionField &chi)
{
  // Apply 4d dslash
  if ( psi.checkerboard == Odd ) {
    this->DhopEO(psi,this->tmp(),DaggerYes);
  } else {
    this->DhopOE(psi,this->tmp(),DaggerYes);
  }
  MeooeDag5D(this->tmp(),chi); 
}

template<class Impl>
void  CayleyFermion5D<Impl>::Mdir (const FermionField &psi, FermionField &chi,int dir,int disp){
  Meo5D(psi,this->tmp());
  // Apply 4d dslash fragment
  this->DhopDir(this->tmp(),chi,dir,disp);
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
  std::vector<Coeff_t> gamma(this->Ls);
  for(int s=0;s<this->Ls;s++) gamma[s] = zdata->gamma[s];
  SetCoefficientsInternal(1.0,gamma,b,c);
}
//Zolo
template<class Impl>
void CayleyFermion5D<Impl>::SetCoefficientsZolotarev(RealD zolo_hi,Approx::zolotarev_data *zdata,RealD b,RealD c)
{
  std::vector<Coeff_t> gamma(this->Ls);
  for(int s=0;s<this->Ls;s++) gamma[s] = zdata->gamma[s];
  SetCoefficientsInternal(zolo_hi,gamma,b,c);
}
//Zolo
template<class Impl>
void CayleyFermion5D<Impl>::SetCoefficientsInternal(RealD zolo_hi,std::vector<Coeff_t> & gamma,RealD b,RealD c)
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
    omega[i] = gamma[i]*zolo_hi; //NB reciprocal relative to Chroma NEF code
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
    Coeff_t delta_d=mass*cee[Ls-1];
    for(int j=0;j<Ls-1;j++) delta_d *= cee[j]/bee[j];
    dee[Ls-1] += delta_d;
  }  

  int inv=1;
  this->MooeeInternalCompute(0,inv,MatpInv,MatmInv);
  this->MooeeInternalCompute(1,inv,MatpInvDag,MatmInvDag);

}


template<class Impl>
void CayleyFermion5D<Impl>::MooeeInternalCompute(int dag, int inv,
						 Vector<iSinglet<Simd> > & Matp,
						 Vector<iSinglet<Simd> > & Matm)
{
  int Ls=this->Ls;

  GridBase *grid = this->FermionRedBlackGrid();
  int LLs = grid->_rdimensions[0];

  if ( LLs == Ls ) return; // Not vectorised in 5th direction

  Eigen::MatrixXcd Pplus  = Eigen::MatrixXcd::Zero(Ls,Ls);
  Eigen::MatrixXcd Pminus = Eigen::MatrixXcd::Zero(Ls,Ls);
  
  for(int s=0;s<Ls;s++){
    Pplus(s,s) = bee[s];
    Pminus(s,s)= bee[s];
  }
  
  for(int s=0;s<Ls-1;s++){
    Pminus(s,s+1) = -cee[s];
  }
  
  for(int s=0;s<Ls-1;s++){
    Pplus(s+1,s) = -cee[s+1];
  }
  Pplus (0,Ls-1) = mass*cee[0];
  Pminus(Ls-1,0) = mass*cee[Ls-1];
  
  Eigen::MatrixXcd PplusMat ;
  Eigen::MatrixXcd PminusMat;
  
  if ( inv ) {
    PplusMat =Pplus.inverse();
    PminusMat=Pminus.inverse();
  } else { 
    PplusMat =Pplus;
    PminusMat=Pminus;
  }
  
  if(dag){
    PplusMat.adjointInPlace();
    PminusMat.adjointInPlace();
  }
  
  typedef typename SiteHalfSpinor::scalar_type scalar_type;
  const int Nsimd=Simd::Nsimd();
  Matp.resize(Ls*LLs);
  Matm.resize(Ls*LLs);

  for(int s2=0;s2<Ls;s2++){
  for(int s1=0;s1<LLs;s1++){
    int istride = LLs;
    int ostride = 1;
    Simd Vp;
    Simd Vm;
    scalar_type *sp = (scalar_type *)&Vp;
    scalar_type *sm = (scalar_type *)&Vm;
    for(int l=0;l<Nsimd;l++){
      if ( switcheroo<Coeff_t>::iscomplex() ) {
	sp[l] = PplusMat (l*istride+s1*ostride,s2);
	sm[l] = PminusMat(l*istride+s1*ostride,s2);
      } else { 
      // if real
	scalar_type tmp;
	tmp = PplusMat (l*istride+s1*ostride,s2);
	sp[l] = scalar_type(tmp.real(),tmp.real());
	tmp = PminusMat(l*istride+s1*ostride,s2);
	sm[l] = scalar_type(tmp.real(),tmp.real());
      }
    }
    Matp[LLs*s2+s1] = Vp;
    Matm[LLs*s2+s1] = Vm;
  }}
}


  FermOpTemplateInstantiate(CayleyFermion5D);
  GparityFermOpTemplateInstantiate(CayleyFermion5D);

}}


