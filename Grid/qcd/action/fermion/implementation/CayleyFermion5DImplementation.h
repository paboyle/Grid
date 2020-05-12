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

#include <Grid/Grid_Eigen_Dense.h>
#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/CayleyFermion5D.h>

NAMESPACE_BEGIN(Grid);

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

///////////////////////////////////////////////////////////////
// Physical surface field utilities
///////////////////////////////////////////////////////////////
template<class Impl>  
void CayleyFermion5D<Impl>::ExportPhysicalFermionSolution(const FermionField &solution5d,FermionField &exported4d)
{
  int Ls = this->Ls;
  FermionField tmp(this->FermionGrid());
  tmp = solution5d;
  conformable(solution5d.Grid(),this->FermionGrid());
  conformable(exported4d.Grid(),this->GaugeGrid());
  axpby_ssp_pminus(tmp, 0., solution5d, 1., solution5d, 0, 0);
  axpby_ssp_pplus (tmp, 1., tmp       , 1., solution5d, 0, Ls-1);
  ExtractSlice(exported4d, tmp, 0, 0);
}
template<class Impl>  
void CayleyFermion5D<Impl>::P(const FermionField &psi, FermionField &chi)
{
  int Ls= this->Ls;
  chi=Zero();
  for(int s=0;s<Ls;s++){
    axpby_ssp_pminus(chi,1.0,chi,1.0,psi,s,s);
    axpby_ssp_pplus (chi,1.0,chi,1.0,psi,s,(s+1)%Ls);
  }
}
template<class Impl>  
void CayleyFermion5D<Impl>::Pdag(const FermionField &psi, FermionField &chi)
{
  int Ls= this->Ls;
  chi=Zero();
  for(int s=0;s<Ls;s++){
    axpby_ssp_pminus(chi,1.0,chi,1.0,psi,s,s);
    axpby_ssp_pplus (chi,1.0,chi,1.0,psi,s,(s-1+Ls)%Ls);
  }
}
template<class Impl>  
void CayleyFermion5D<Impl>::ExportPhysicalFermionSource(const FermionField &solution5d,FermionField &exported4d)
{
  int Ls = this->Ls;
  FermionField tmp(this->FermionGrid());
  tmp = solution5d;
  conformable(solution5d.Grid(),this->FermionGrid());
  conformable(exported4d.Grid(),this->GaugeGrid());
  axpby_ssp_pplus (tmp, 0., solution5d, 1., solution5d, 0, 0);
  axpby_ssp_pminus(tmp, 1., tmp       , 1., solution5d, 0, Ls-1);
  ExtractSlice(exported4d, tmp, 0, 0);
}
template<class Impl>
void CayleyFermion5D<Impl>::ImportUnphysicalFermion(const FermionField &input4d,FermionField &imported5d)
{
  int Ls = this->Ls;
  FermionField tmp(this->FermionGrid());
  conformable(imported5d.Grid(),this->FermionGrid());
  conformable(input4d.Grid()   ,this->GaugeGrid());
  tmp = Zero();
  InsertSlice(input4d, tmp, 0   , 0);
  InsertSlice(input4d, tmp, Ls-1, 0);
  axpby_ssp_pplus (tmp, 0., tmp, 1., tmp, 0, 0);
  axpby_ssp_pminus(tmp, 0., tmp, 1., tmp, Ls-1, Ls-1);
  imported5d=tmp;
}

template<class Impl>  
void CayleyFermion5D<Impl>::ImportPhysicalFermionSource(const FermionField &input4d,FermionField &imported5d)
{
  int Ls = this->Ls;
  FermionField tmp(this->FermionGrid());
  conformable(imported5d.Grid(),this->FermionGrid());
  conformable(input4d.Grid()   ,this->GaugeGrid());
  tmp = Zero();
  InsertSlice(input4d, tmp, 0   , 0);
  InsertSlice(input4d, tmp, Ls-1, 0);
  axpby_ssp_pplus (tmp, 0., tmp, 1., tmp, 0, 0);
  axpby_ssp_pminus(tmp, 0., tmp, 1., tmp, Ls-1, Ls-1);
  Dminus(tmp,imported5d);
}
template<class Impl>  
void CayleyFermion5D<Impl>::Dminus(const FermionField &psi, FermionField &chi)
{
  int Ls=this->Ls;

  FermionField tmp_f(this->FermionGrid());
  this->DW(psi,tmp_f,DaggerNo);

  for(int s=0;s<Ls;s++){
    axpby_ssp(chi,Coeff_t(1.0),psi,-cs[s],tmp_f,s,s);// chi = (1-c[s] D_W) psi
  }
}
template<class Impl>  
void CayleyFermion5D<Impl>::DminusDag(const FermionField &psi, FermionField &chi)
{
  int Ls=this->Ls;

  FermionField tmp_f(this->FermionGrid());
  this->DW(psi,tmp_f,DaggerYes);

  for(int s=0;s<Ls;s++){
    axpby_ssp(chi,Coeff_t(1.0),psi,conjugate(-cs[s]),tmp_f,s,s);// chi = (1-c[s] D_W) psi
  }
}

template<class Impl> void CayleyFermion5D<Impl>::CayleyReport(void)
{
  this->Report();
  Coordinate latt = GridDefaultLatt();          
  RealD volume = this->Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt[mu];
  RealD NP     = this->_FourDimGrid->_Nprocessors;
  if ( M5Dcalls > 0 ) {
    std::cout << GridLogMessage << "#### M5D calls report " << std::endl;
    std::cout << GridLogMessage << "CayleyFermion5D Number of M5D Calls     : " << M5Dcalls   << std::endl;
    std::cout << GridLogMessage << "CayleyFermion5D ComputeTime/Calls       : " << M5Dtime / M5Dcalls << " us" << std::endl;

    // Flops = 10.0*(Nc*Ns) *Ls*vol
    RealD mflops = 10.0*(Nc*Ns)*volume*M5Dcalls/M5Dtime/2; // 2 for red black counting
    std::cout << GridLogMessage << "Average mflops/s per call                : " << mflops << std::endl;
    std::cout << GridLogMessage << "Average mflops/s per call per rank       : " << mflops/NP << std::endl;

    // Bytes = sizeof(Real) * (Nc*Ns*Nreim) * Ls * vol * (read+write) (/2 for red black counting)
    // read = 2 ( psi[ss+s+1] and psi[ss+s-1] count as 1 )
    // write = 1
    RealD Gbytes = sizeof(Real) * (Nc*Ns*2) * volume * 3 /2. * 1.e-9;
    std::cout << GridLogMessage << "Average bandwidth (GB/s)                 : " << Gbytes/M5Dtime*M5Dcalls*1.e6 << std::endl;
  }

  if ( MooeeInvCalls > 0 ) {

    std::cout << GridLogMessage << "#### MooeeInv calls report " << std::endl;
    std::cout << GridLogMessage << "CayleyFermion5D Number of MooeeInv Calls     : " << MooeeInvCalls   << std::endl;
    std::cout << GridLogMessage << "CayleyFermion5D ComputeTime/Calls            : " << MooeeInvTime / MooeeInvCalls << " us" << std::endl;
#ifdef GRID_NVCC
    RealD mflops = ( -16.*Nc*Ns+this->Ls*(1.+18.*Nc*Ns) )*volume*MooeeInvCalls/MooeeInvTime/2; // 2 for red black counting
    std::cout << GridLogMessage << "Average mflops/s per call                : " << mflops << std::endl;
    std::cout << GridLogMessage << "Average mflops/s per call per rank       : " << mflops/NP << std::endl;
#else
    // Flops = MADD * Ls *Ls *4dvol * spin/colour/complex
    RealD mflops = 2.0*24*this->Ls*volume*MooeeInvCalls/MooeeInvTime/2; // 2 for red black counting
    std::cout << GridLogMessage << "Average mflops/s per call                : " << mflops << std::endl;
    std::cout << GridLogMessage << "Average mflops/s per call per rank       : " << mflops/NP << std::endl;
#endif
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
void CayleyFermion5D<Impl>::M5D   (const FermionField &psi, FermionField &chi)
{
  int Ls=this->Ls;
  Vector<Coeff_t> diag (Ls,1.0);
  Vector<Coeff_t> upper(Ls,-1.0); upper[Ls-1]=mass;
  Vector<Coeff_t> lower(Ls,-1.0); lower[0]   =mass;
  M5D(psi,chi,chi,lower,diag,upper);
}
template<class Impl>
void CayleyFermion5D<Impl>::Meooe5D    (const FermionField &psi, FermionField &Din)
{
  int Ls=this->Ls;
  Vector<Coeff_t> diag = bs;
  Vector<Coeff_t> upper= cs;
  Vector<Coeff_t> lower= cs; 
  upper[Ls-1]=-mass*upper[Ls-1];
  lower[0]   =-mass*lower[0];
  M5D(psi,psi,Din,lower,diag,upper);
}
// FIXME Redunant with the above routine; check this and eliminate
template<class Impl> void CayleyFermion5D<Impl>::Meo5D     (const FermionField &psi, FermionField &chi)
{
  int Ls=this->Ls;
  Vector<Coeff_t> diag = beo;
  Vector<Coeff_t> upper(Ls);
  Vector<Coeff_t> lower(Ls);
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
  Vector<Coeff_t> diag = bee;
  Vector<Coeff_t> upper(Ls);
  Vector<Coeff_t> lower(Ls);
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
  Vector<Coeff_t> diag = bee;
  Vector<Coeff_t> upper(Ls);
  Vector<Coeff_t> lower(Ls);

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
  // Conjugate the terms 
  for (int s=0;s<Ls;s++){
    diag[s] =conjugate(diag[s]);
    upper[s]=conjugate(upper[s]);
    lower[s]=conjugate(lower[s]);
  }
  M5Ddag(psi,psi,chi,lower,diag,upper);
}

template<class Impl>
void CayleyFermion5D<Impl>::M5Ddag (const FermionField &psi, FermionField &chi)
{
  int Ls=this->Ls;
  Vector<Coeff_t> diag(Ls,1.0);
  Vector<Coeff_t> upper(Ls,-1.0);
  Vector<Coeff_t> lower(Ls,-1.0);
  upper[Ls-1]=-mass*upper[Ls-1];
  lower[0]   =-mass*lower[0];
  M5Ddag(psi,chi,chi,lower,diag,upper);
}

template<class Impl>
void CayleyFermion5D<Impl>::MeooeDag5D    (const FermionField &psi, FermionField &Din)
{
  int Ls=this->Ls;
  Vector<Coeff_t> diag =bs;
  Vector<Coeff_t> upper=cs;
  Vector<Coeff_t> lower=cs; 

  for (int s=0;s<Ls;s++){
    if ( s== 0 ) {
      upper[s] = cs[s+1];
      lower[s] =-mass*cs[Ls-1];
    } else if ( s==(Ls-1) ) { 
      upper[s] =-mass*cs[0];
      lower[s] = cs[s-1];
    } else { 
      upper[s] = cs[s+1];
      lower[s] = cs[s-1];
    }
    upper[s] = conjugate(upper[s]);
    lower[s] = conjugate(lower[s]);
    diag[s]  = conjugate(diag[s]);
  }
  M5Ddag(psi,psi,Din,lower,diag,upper);
}

template<class Impl>
void CayleyFermion5D<Impl>::M    (const FermionField &psi, FermionField &chi)
{
  FermionField Din(psi.Grid());
  
  // Assemble Din
  Meooe5D(psi,Din);
  
  this->DW(Din,chi,DaggerNo);
  // ((b D_W + D_w hop terms +1) on s-diag
  axpby(chi,1.0,1.0,chi,psi); 
  
  M5D(psi,chi);
}

template<class Impl>
void CayleyFermion5D<Impl>::Mdag (const FermionField &psi, FermionField &chi)
{
  // Under adjoint
  //D1+        D1- P-    ->   D1+^dag   P+ D2-^dag
  //D2- P+     D2+            P-D1-^dag D2+dag
  
  FermionField Din(psi.Grid());
  // Apply Dw
  this->DW(psi,Din,DaggerYes); 
  
  MeooeDag5D(Din,chi);
  
  M5Ddag(psi,chi);
  // ((b D_W + D_w hop terms +1) on s-diag
  axpby (chi,1.0,1.0,chi,psi); 
}

// half checkerboard operations
template<class Impl>
void CayleyFermion5D<Impl>::Meooe       (const FermionField &psi, FermionField &chi)
{
  Meooe5D(psi,this->tmp()); 

  if ( psi.Checkerboard() == Odd ) {
    this->DhopEO(this->tmp(),chi,DaggerNo);
  } else {
    this->DhopOE(this->tmp(),chi,DaggerNo);
  }
}

template<class Impl>
void CayleyFermion5D<Impl>::MeooeDag    (const FermionField &psi, FermionField &chi)
{
  // Apply 4d dslash
  if ( psi.Checkerboard() == Odd ) {
    this->DhopEO(psi,this->tmp(),DaggerYes);
  } else {
    this->DhopOE(psi,this->tmp(),DaggerYes);
  }
  MeooeDag5D(this->tmp(),chi); 
}

template<class Impl>
void  CayleyFermion5D<Impl>::Mdir (const FermionField &psi, FermionField &chi,int dir,int disp)
{
  FermionField tmp(psi.Grid());
  Meo5D(psi,tmp);
  this->DhopDir(tmp,chi,dir,disp);
}
template<class Impl>
void  CayleyFermion5D<Impl>::MdirAll(const FermionField &psi, std::vector<FermionField> &out)
{
  FermionField tmp(psi.Grid());
  Meo5D(psi,tmp);
  this->DhopDirAll(tmp,out);
}

// force terms; five routines; default to Dhop on diagonal
template<class Impl>
void CayleyFermion5D<Impl>::MDeriv  (GaugeField &mat,const FermionField &U,const FermionField &V,int dag)
{
  FermionField Din(V.Grid());
  
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
  FermionField Din(V.Grid());
  
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
  FermionField Din(V.Grid());
  
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
  Vector<Coeff_t> gamma(this->Ls);
  for(int s=0;s<this->Ls;s++) gamma[s] = zdata->gamma[s];
  SetCoefficientsInternal(1.0,gamma,b,c);
}
//Zolo
template<class Impl>
void CayleyFermion5D<Impl>::SetCoefficientsZolotarev(RealD zolo_hi,Approx::zolotarev_data *zdata,RealD b,RealD c)
{
  Vector<Coeff_t> gamma(this->Ls);
  for(int s=0;s<this->Ls;s++) gamma[s] = zdata->gamma[s];
  SetCoefficientsInternal(zolo_hi,gamma,b,c);
}
//Zolo
template<class Impl>
void CayleyFermion5D<Impl>::SetCoefficientsInternal(RealD zolo_hi,Vector<Coeff_t> & gamma,RealD b,RealD c)
{
  int Ls=this->Ls;

  ///////////////////////////////////////////////////////////
  // The Cayley coeffs (unprec)
  ///////////////////////////////////////////////////////////
  assert(gamma.size()==Ls);

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
  _b = b;
  _c = c;
  _gamma  = gamma; // Save the parameters so we can change mass later.
  _zolo_hi= zolo_hi;
  for(int i=0; i < Ls; i++){
    as[i] = 1.0;
    omega[i] = _gamma[i]*_zolo_hi; //NB reciprocal relative to Chroma NEF code
    assert(omega[i]!=Coeff_t(0.0));
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
    assert(bee[i]!=Coeff_t(0.0));
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

      assert(bee[i]!=Coeff_t(0.0));
      assert(bee[0]!=Coeff_t(0.0));
      
      lee[i] =-cee[i+1]/bee[i]; // sub-diag entry on the ith column
      
      leem[i]=mass*cee[Ls-1]/bee[0];
      for(int j=0;j<i;j++) {
	assert(bee[j+1]!=Coeff_t(0.0));
	leem[i]*= aee[j]/bee[j+1];
      }
      
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
    for(int j=0;j<Ls-1;j++) {
      assert(bee[j] != Coeff_t(0.0));
      delta_d *= cee[j]/bee[j];
    }
    dee[Ls-1] += delta_d;
  }  

  //  int inv=1;
  //  this->MooeeInternalCompute(0,inv,MatpInv,MatmInv);
  //  this->MooeeInternalCompute(1,inv,MatpInvDag,MatmInvDag);
}


template <class Impl>
void CayleyFermion5D<Impl>::ContractJ5q(FermionField &q_in,ComplexField &J5q)
{
  conformable(this->GaugeGrid(), J5q.Grid());
  conformable(q_in.Grid(), this->FermionGrid());
  Gamma G5(Gamma::Algebra::Gamma5);
  // 4d field
  int Ls = this->Ls;
  FermionField psi(this->GaugeGrid());
  FermionField p_plus (this->GaugeGrid());
  FermionField p_minus(this->GaugeGrid());
  FermionField p(this->GaugeGrid());

  ExtractSlice(p_plus , q_in, Ls/2-1 , 0);
  ExtractSlice(p_minus, q_in, Ls/2   , 0);
  p_plus = p_plus + G5*p_plus;
  p_minus= p_minus - G5*p_minus;
  p=0.5*(p_plus+p_minus);
  J5q = localInnerProduct(p,p);
}

template <class Impl>
void CayleyFermion5D<Impl>::ContractJ5q(PropagatorField &q_in,ComplexField &J5q)
{
  conformable(this->GaugeGrid(), J5q.Grid());
  conformable(q_in.Grid(), this->FermionGrid());
  Gamma G5(Gamma::Algebra::Gamma5);
  // 4d field
  int Ls = this->Ls;
  PropagatorField psi(this->GaugeGrid());
  PropagatorField p_plus (this->GaugeGrid());
  PropagatorField p_minus(this->GaugeGrid());
  PropagatorField p(this->GaugeGrid());

  ExtractSlice(p_plus , q_in, Ls/2-1 , 0);
  ExtractSlice(p_minus, q_in, Ls/2   , 0);
  p_plus = p_plus + G5*p_plus;
  p_minus= p_minus - G5*p_minus;
  p=0.5*(p_plus+p_minus);
  J5q = localInnerProduct(p,p);
}

#define Pp(Q) (0.5*(Q+g5*Q))
#define Pm(Q) (0.5*(Q-g5*Q))
#define Q_4d(Q) (Pm((Q)[0]) + Pp((Q)[Ls-1]))
#define TopRowWithSource(Q) (phys_src + (1.0-mass)*Q_4d(Q))

template <class Impl> 
void CayleyFermion5D<Impl>::ContractConservedCurrent( PropagatorField &q_in_1,
						      PropagatorField &q_in_2,
						      PropagatorField &q_out,
						      PropagatorField &phys_src,
						      Current curr_type,
						      unsigned int mu)
{
#ifndef GRID_NVCC
  Gamma::Algebra Gmu [] = {
    Gamma::Algebra::GammaX,
    Gamma::Algebra::GammaY,
    Gamma::Algebra::GammaZ,
    Gamma::Algebra::GammaT,
    Gamma::Algebra::Gamma5
  };

  auto UGrid= this->GaugeGrid();
  auto FGrid= this->FermionGrid();
  RealD sgn=1.0;
  if ( curr_type == Current::Axial ) sgn = -1.0;

  int Ls = this->Ls;

  std::vector<PropagatorField> L_Q(Ls,UGrid); 
  std::vector<PropagatorField> R_Q(Ls,UGrid); 
  for(int s=0;s<Ls;s++){
    ExtractSlice(L_Q[s], q_in_1, s , 0);
    ExtractSlice(R_Q[s], q_in_2, s , 0);
  }

  Gamma g5(Gamma::Algebra::Gamma5);
  PropagatorField C(UGrid); 
  PropagatorField p5d(UGrid); 
  PropagatorField us_p5d(UGrid); 
  PropagatorField gp5d(UGrid); 
  PropagatorField gus_p5d(UGrid); 

  PropagatorField L_TmLsGq0(UGrid); 
  PropagatorField L_TmLsTmp(UGrid);
  PropagatorField R_TmLsGq0(UGrid); 
  PropagatorField R_TmLsTmp(UGrid);
  {
    PropagatorField TermA(UGrid);
    PropagatorField TermB(UGrid);
    PropagatorField TermC(UGrid);
    PropagatorField TermD(UGrid);
    TermA = (Pp(Q_4d(L_Q)));
    TermB = (Pm(Q_4d(L_Q)));
    TermC = (Pm(TopRowWithSource(L_Q)));
    TermD = (Pp(TopRowWithSource(L_Q)));

    L_TmLsGq0 = (TermD - TermA + TermB);
    L_TmLsTmp = (TermC - TermB + TermA);

    TermA = (Pp(Q_4d(R_Q)));
    TermB = (Pm(Q_4d(R_Q)));
    TermC = (Pm(TopRowWithSource(R_Q)));
    TermD = (Pp(TopRowWithSource(R_Q)));

    R_TmLsGq0 = (TermD - TermA + TermB);
    R_TmLsTmp = (TermC - TermB + TermA);
  }

  std::vector<PropagatorField> R_TmLsGq(Ls,UGrid);
  std::vector<PropagatorField> L_TmLsGq(Ls,UGrid);
  for(int s=0;s<Ls;s++){
    R_TmLsGq[s] = (Pm((R_Q)[(s)]) + Pp((R_Q)[((s)-1+Ls)%Ls]));
    L_TmLsGq[s] = (Pm((L_Q)[(s)]) + Pp((L_Q)[((s)-1+Ls)%Ls]));
  }

  Gamma gmu=Gamma(Gmu[mu]);

  q_out = Zero();
  PropagatorField tmp(UGrid); 
  for(int s=0;s<Ls;s++){

    int sp = (s+1)%Ls;
    int sr = Ls-1-s;
    int srp= (sr+1)%Ls;

    // Mobius parameters
    auto b=this->bs[s];
    auto c=this->cs[s];
    auto bpc = 1.0/(b+c);  // -0.5 factor in gauge links
    if (s == 0) {
      p5d    =(b*Pm(L_TmLsGq[Ls-1])+ c*Pp(L_TmLsGq[Ls-1]) + b*Pp(L_TmLsTmp)   + c*Pm(L_TmLsTmp     ));
      tmp    =(b*Pm(R_TmLsGq0)     + c*Pp(R_TmLsGq0     ) + b*Pp(R_TmLsGq[1]) + c*Pm(R_TmLsGq[1]));
    } else if (s == Ls-1) {
      p5d    =(b*Pm(L_TmLsGq0)     + c*Pp(L_TmLsGq0     ) + b*Pp(L_TmLsGq[1]) + c*Pm(L_TmLsGq[1]));
      tmp    =(b*Pm(R_TmLsGq[Ls-1])+ c*Pp(R_TmLsGq[Ls-1]) + b*Pp(R_TmLsTmp)   + c*Pm(R_TmLsTmp   ));
    } else {
      p5d    =(b*Pm(L_TmLsGq[sr]) + c*Pp(L_TmLsGq[sr])+ b*Pp(L_TmLsGq[srp])+ c*Pm(L_TmLsGq[srp]));
      tmp    =(b*Pm(R_TmLsGq[s])  + c*Pp(R_TmLsGq[s]) + b*Pp(R_TmLsGq[sp ])+ c*Pm(R_TmLsGq[sp]));
    }
    tmp    = Cshift(tmp,mu,1);
    Impl::multLinkField(us_p5d,this->Umu,tmp,mu);
    
    gp5d=g5*p5d*g5;
    gus_p5d=gmu*us_p5d;

    C = bpc*(adj(gp5d)*us_p5d);
    C-= bpc*(adj(gp5d)*gus_p5d);

    if (s == 0) {
      p5d    =(b*Pm(R_TmLsGq0)     + c*Pp(R_TmLsGq0  )    + b*Pp(R_TmLsGq[1]) + c*Pm(R_TmLsGq[1]));
      tmp    =(b*Pm(L_TmLsGq[Ls-1])+ c*Pp(L_TmLsGq[Ls-1]) + b*Pp(L_TmLsTmp)   + c*Pm(L_TmLsTmp  ));
    } else if (s == Ls-1) {
      p5d    =(b*Pm(R_TmLsGq[Ls-1])+ c*Pp(R_TmLsGq[Ls-1]) + b*Pp(R_TmLsTmp)   + c*Pm(R_TmLsTmp  ));
      tmp    =(b*Pm(L_TmLsGq0)     + c*Pp(L_TmLsGq0  )    + b*Pp(L_TmLsGq[1]) + c*Pm(L_TmLsGq[1]));
    } else {
      p5d    =(b*Pm(R_TmLsGq[s])  + c*Pp(R_TmLsGq[s])  + b*Pp(R_TmLsGq[sp ])+ c*Pm(R_TmLsGq[sp]));
      tmp    =(b*Pm(L_TmLsGq[sr]) + c*Pp(L_TmLsGq[sr]) + b*Pp(L_TmLsGq[srp])+ c*Pm(L_TmLsGq[srp]));
    }
    tmp    = Cshift(tmp,mu,1);
    Impl::multLinkField(us_p5d,this->Umu,tmp,mu);

    gp5d=gmu*p5d;
    gus_p5d=g5*us_p5d*g5;

    C-= bpc*(adj(gus_p5d)*gp5d);
    C-= bpc*(adj(gus_p5d)*p5d);

    if (s < Ls/2) q_out += sgn*C;
    else          q_out +=     C;
    
  }
#endif
}

template <class Impl>
void CayleyFermion5D<Impl>::SeqConservedCurrent(PropagatorField &q_in, 
                                                PropagatorField &q_out,
                                                PropagatorField &phys_src,
                                                Current curr_type, 
                                                unsigned int mu,
                                                unsigned int tmin, 
                                                unsigned int tmax,
						ComplexField &ph)// Complex phase factor
{
  assert(mu>=0);
  assert(mu<Nd);


#if 0
  int tshift = (mu == Nd-1) ? 1 : 0;
  ////////////////////////////////////////////////
  // SHAMIR CASE 
  ////////////////////////////////////////////////
  int Ls = this->Ls;
  auto UGrid= this->GaugeGrid();
  auto FGrid= this->FermionGrid();
  Gamma::Algebra Gmu [] = {
    Gamma::Algebra::GammaX,
    Gamma::Algebra::GammaY,
    Gamma::Algebra::GammaZ,
    Gamma::Algebra::GammaT
  };
  Gamma gmu=Gamma(Gmu[mu]);

  PropagatorField L_Q(UGrid); 
  PropagatorField R_Q(UGrid); 

  PropagatorField tmp(UGrid);
  PropagatorField Utmp(UGrid);
  LatticeInteger zz (UGrid);   zz=0.0;
  LatticeInteger lcoor(UGrid); LatticeCoordinate(lcoor,Nd-1);
  for (int s=0;s<Ls;s++) {

    RealD G_s = (curr_type == Current::Axial  ) ? ((s < Ls/2) ? -1 : 1) : 1;

    ExtractSlice(R_Q, q_in, s , 0);

    tmp    = Cshift(R_Q,mu,1);
    Impl::multLinkField(Utmp,this->Umu,tmp,mu);
    tmp    = G_s*( Utmp*ph - gmu*Utmp*ph ); // Forward hop
    tmp    = where((lcoor>=tmin),tmp,zz); // Mask the time
    tmp    = where((lcoor<=tmax),tmp,zz);
    L_Q = tmp;

    tmp    = R_Q*ph;
    tmp    = Cshift(tmp,mu,-1);
    Impl::multLinkField(Utmp,this->Umu,tmp,mu+Nd);// Adjoint link
    tmp    = -G_s*( Utmp + gmu*Utmp ); 
    tmp    = where((lcoor>=tmin+tshift),tmp,zz); // Mask the time 
    tmp    = where((lcoor<=tmax+tshift),tmp,zz); // Position of current complicated
    L_Q= L_Q+tmp;

    InsertSlice(L_Q, q_out, s , 0);
  }
#endif

#ifndef GRID_NVCC
  int tshift = (mu == Nd-1) ? 1 : 0;
  ////////////////////////////////////////////////
  // GENERAL CAYLEY CASE
  ////////////////////////////////////////////////
  Gamma::Algebra Gmu [] = {
    Gamma::Algebra::GammaX,
    Gamma::Algebra::GammaY,
    Gamma::Algebra::GammaZ,
    Gamma::Algebra::GammaT,
    Gamma::Algebra::Gamma5
  };
  Gamma gmu=Gamma(Gmu[mu]);
  Gamma g5(Gamma::Algebra::Gamma5);

  int Ls = this->Ls;
  auto UGrid= this->GaugeGrid();
  auto FGrid= this->FermionGrid();

  std::vector<PropagatorField> R_Q(Ls,UGrid); 
  PropagatorField L_Q(UGrid); 
  PropagatorField tmp(UGrid);
  PropagatorField Utmp(UGrid);

  LatticeInteger zz (UGrid);   zz=0.0;
  LatticeInteger lcoor(UGrid); LatticeCoordinate(lcoor,Nd-1);

  for(int s=0;s<Ls;s++){
    ExtractSlice(R_Q[s], q_in, s , 0);
  }

  PropagatorField R_TmLsGq0(UGrid); 
  PropagatorField R_TmLsTmp(UGrid);
  {
    PropagatorField TermA(UGrid);
    PropagatorField TermB(UGrid);
    PropagatorField TermC(UGrid);
    PropagatorField TermD(UGrid);

    TermA = (Pp(Q_4d(R_Q)));
    TermB = (Pm(Q_4d(R_Q)));
    TermC = (Pm(TopRowWithSource(R_Q)));
    TermD = (Pp(TopRowWithSource(R_Q)));

    R_TmLsGq0 = (TermD - TermA + TermB);
    R_TmLsTmp = (TermC - TermB + TermA);
  }

  std::vector<PropagatorField> R_TmLsGq(Ls,UGrid);
  for(int s=0;s<Ls;s++){
    R_TmLsGq[s] = (Pm((R_Q)[(s)]) + Pp((R_Q)[((s)-1+Ls)%Ls]));
  }

  std::vector<RealD> G_s(Ls,1.0);
  if ( curr_type == Current::Axial ) {
    for(int s=0;s<Ls/2;s++){
      G_s[s] = -1.0;
    }
  }

  for(int s=0;s<Ls;s++){

    int sp = (s+1)%Ls;
    int sr = Ls-1-s;
    int srp= (sr+1)%Ls;

    // Mobius parameters
    auto b=this->bs[s];
    auto c=this->cs[s];
    //    auto bpc = G_s[s]*1.0/(b+c);  // -0.5 factor in gauge links

    if (s == 0) {
      tmp    =(b*Pm(R_TmLsGq0)     + c*Pp(R_TmLsGq0     ) + b*Pp(R_TmLsGq[1]) + c*Pm(R_TmLsGq[1]));
    } else if (s == Ls-1) {
      tmp    =(b*Pm(R_TmLsGq[Ls-1])+ c*Pp(R_TmLsGq[Ls-1]) + b*Pp(R_TmLsTmp)   + c*Pm(R_TmLsTmp   ));
    } else {
      tmp    =(b*Pm(R_TmLsGq[s])  + c*Pp(R_TmLsGq[s]) + b*Pp(R_TmLsGq[sp ])+ c*Pm(R_TmLsGq[sp]));
    }

    tmp    = Cshift(tmp,mu,1);
    Impl::multLinkField(Utmp,this->Umu,tmp,mu);
    tmp    = G_s[s]*( Utmp*ph - gmu*Utmp*ph ); // Forward hop
    tmp    = where((lcoor>=tmin),tmp,zz); // Mask the time 
    L_Q    = where((lcoor<=tmax),tmp,zz); // Position of current complicated

    if (s == 0) {
      tmp    =(b*Pm(R_TmLsGq0)     + c*Pp(R_TmLsGq0  )    + b*Pp(R_TmLsGq[1]) + c*Pm(R_TmLsGq[1]));
    } else if (s == Ls-1) {
      tmp    =(b*Pm(R_TmLsGq[Ls-1])+ c*Pp(R_TmLsGq[Ls-1]) + b*Pp(R_TmLsTmp)   + c*Pm(R_TmLsTmp  ));
    } else {
      tmp    =(b*Pm(R_TmLsGq[s])   + c*Pp(R_TmLsGq[s])    + b*Pp(R_TmLsGq[sp])+ c*Pm(R_TmLsGq[sp]));
    }
    tmp    = tmp *ph;
    tmp    = Cshift(tmp,mu,-1);
    Impl::multLinkField(Utmp,this->Umu,tmp,mu+Nd); // Adjoint link
    tmp = -G_s[s]*( Utmp + gmu*Utmp );
    tmp    = where((lcoor>=tmin+tshift),tmp,zz); // Mask the time 
    L_Q   += where((lcoor<=tmax+tshift),tmp,zz); // Position of current complicated

    InsertSlice(L_Q, q_out, s , 0);
  }
#endif
}
#undef Pp
#undef Pm
#undef Q_4d
#undef TopRowWithSource



#if 0
template<class Impl>
void CayleyFermion5D<Impl>::MooeeInternalCompute(int dag, int inv,
						 Vector<iSinglet<Simd> > & Matp,
						 Vector<iSinglet<Simd> > & Matm)
{
  int Ls=this->Ls;

  GridBase *grid = this->FermionRedBlackGrid();
  int LLs = grid->_rdimensions[0];

  if ( LLs == Ls ) {
    return; // Not vectorised in 5th direction
  }

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
#endif

NAMESPACE_END(Grid);


