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
#include <Grid.h>


namespace Grid {
namespace QCD {
  /*
   * Dense matrix versions of routines
   */
template<class Impl>
void CayleyFermion5D<Impl>::MooeeInvDag (const FermionField &psi, FermionField &chi)
{
  this->MooeeInternal(psi,chi,DaggerYes,InverseYes);
}
  
template<class Impl>
void CayleyFermion5D<Impl>::MooeeInv(const FermionField &psi, FermionField &chi)
{
  this->MooeeInternal(psi,chi,DaggerNo,InverseYes);
}
template<class Impl>  
void CayleyFermion5D<Impl>::M5D(const FermionField &psi,
				const FermionField &phi, 
				FermionField &chi,
				std::vector<Coeff_t> &lower,
				std::vector<Coeff_t> &diag,
				std::vector<Coeff_t> &upper)
{
  GridBase *grid=psi._grid;
  int Ls   = this->Ls;
  int LLs  = grid->_rdimensions[0];
  const int nsimd= Simd::Nsimd();

  Vector<iSinglet<Simd> > u(LLs);
  Vector<iSinglet<Simd> > l(LLs);
  Vector<iSinglet<Simd> > d(LLs);

  assert(Ls/LLs==nsimd);
  assert(phi.checkerboard == psi.checkerboard);

  chi.checkerboard=psi.checkerboard;

  // just directly address via type pun
  typedef typename Simd::scalar_type scalar_type;
  scalar_type * u_p = (scalar_type *)&u[0];
  scalar_type * l_p = (scalar_type *)&l[0];
  scalar_type * d_p = (scalar_type *)&d[0];

  for(int o=0;o<LLs;o++){ // outer
  for(int i=0;i<nsimd;i++){ //inner
    int s  = o+i*LLs;
    int ss = o*nsimd+i;
    u_p[ss] = upper[s];
    l_p[ss] = lower[s];
    d_p[ss] = diag[s];
  }}


  M5Dcalls++;
  M5Dtime-=usecond();

  assert(Nc==3);

PARALLEL_FOR_LOOP
  for(int ss=0;ss<grid->oSites();ss+=LLs){ // adds LLs
#if 0
      alignas(64) SiteHalfSpinor hp;
      alignas(64) SiteHalfSpinor hm;
      alignas(64) SiteSpinor fp;
      alignas(64) SiteSpinor fm;

      for(int v=0;v<LLs;v++){

	int vp=(v+1)%LLs;
	int vm=(v+LLs-1)%LLs;

	spProj5m(hp,psi[ss+vp]);
	spProj5p(hm,psi[ss+vm]);

	if ( vp<=v ) rotate(hp,hp,1);
	if ( vm>=v ) rotate(hm,hm,nsimd-1);
	
	hp=0.5*hp;
        hm=0.5*hm;

	spRecon5m(fp,hp);
	spRecon5p(fm,hm);

	chi[ss+v] = d[v]*phi[ss+v];
	chi[ss+v] = chi[ss+v]     +u[v]*fp;
	chi[ss+v] = chi[ss+v]     +l[v]*fm;

      }
#else
      for(int v=0;v<LLs;v++){

	vprefetch(psi[ss+v+LLs]);
	//	vprefetch(phi[ss+v+LLs]);

	int vp= (v==LLs-1) ? 0     : v+1;
	int vm= (v==0    ) ? LLs-1 : v-1;
	
	Simd hp_00 = psi[ss+vp]()(2)(0); 
	Simd hp_01 = psi[ss+vp]()(2)(1); 
	Simd hp_02 = psi[ss+vp]()(2)(2); 
	Simd hp_10 = psi[ss+vp]()(3)(0); 
	Simd hp_11 = psi[ss+vp]()(3)(1); 
	Simd hp_12 = psi[ss+vp]()(3)(2); 
	
	Simd hm_00 = psi[ss+vm]()(0)(0); 
	Simd hm_01 = psi[ss+vm]()(0)(1); 
	Simd hm_02 = psi[ss+vm]()(0)(2); 
	Simd hm_10 = psi[ss+vm]()(1)(0); 
	Simd hm_11 = psi[ss+vm]()(1)(1); 
	Simd hm_12 = psi[ss+vm]()(1)(2); 

	//	if ( ss==0) std::cout << " hp_00 " <<hp_00<<std::endl;
	//	if ( ss==0) std::cout << " hm_00 " <<hm_00<<std::endl;

	if ( vp<=v ) {
	  hp_00.v = Optimization::Rotate::tRotate<2>(hp_00.v);
	  hp_01.v = Optimization::Rotate::tRotate<2>(hp_01.v);
	  hp_02.v = Optimization::Rotate::tRotate<2>(hp_02.v);
	  hp_10.v = Optimization::Rotate::tRotate<2>(hp_10.v);
	  hp_11.v = Optimization::Rotate::tRotate<2>(hp_11.v);
	  hp_12.v = Optimization::Rotate::tRotate<2>(hp_12.v);
	}
	if ( vm>=v ) {
	  hm_00.v = Optimization::Rotate::tRotate<2*Simd::Nsimd()-2>(hm_00.v);
	  hm_01.v = Optimization::Rotate::tRotate<2*Simd::Nsimd()-2>(hm_01.v);
	  hm_02.v = Optimization::Rotate::tRotate<2*Simd::Nsimd()-2>(hm_02.v);
	  hm_10.v = Optimization::Rotate::tRotate<2*Simd::Nsimd()-2>(hm_10.v);
	  hm_11.v = Optimization::Rotate::tRotate<2*Simd::Nsimd()-2>(hm_11.v);
	  hm_12.v = Optimization::Rotate::tRotate<2*Simd::Nsimd()-2>(hm_12.v);
	}

	/*
	if ( ss==0) std::cout << " dphi_00 " <<d[v]()()() * phi[ss+v]()(0)(0) <<std::endl;
	if ( ss==0) std::cout << " dphi_10 " <<d[v]()()() * phi[ss+v]()(1)(0) <<std::endl;
	if ( ss==0) std::cout << " dphi_20 " <<d[v]()()() * phi[ss+v]()(2)(0) <<std::endl;
	if ( ss==0) std::cout << " dphi_30 " <<d[v]()()() * phi[ss+v]()(3)(0) <<std::endl;
	*/	
	Simd p_00  = d[v]()()() * phi[ss+v]()(0)(0)  + l[v]()()()*hm_00; 
	Simd p_01  = d[v]()()() * phi[ss+v]()(0)(1)  + l[v]()()()*hm_01; 
	Simd p_02  = d[v]()()() * phi[ss+v]()(0)(2)  + l[v]()()()*hm_02; 
	Simd p_10  = d[v]()()() * phi[ss+v]()(1)(0)  + l[v]()()()*hm_10; 
	Simd p_11  = d[v]()()() * phi[ss+v]()(1)(1)  + l[v]()()()*hm_11; 
	Simd p_12  = d[v]()()() * phi[ss+v]()(1)(2)  + l[v]()()()*hm_12; 
	Simd p_20  = d[v]()()() * phi[ss+v]()(2)(0)  + u[v]()()()*hp_00; 
	Simd p_21  = d[v]()()() * phi[ss+v]()(2)(1)  + u[v]()()()*hp_01; 
	Simd p_22  = d[v]()()() * phi[ss+v]()(2)(2)  + u[v]()()()*hp_02;  
	Simd p_30  = d[v]()()() * phi[ss+v]()(3)(0)  + u[v]()()()*hp_10; 
	Simd p_31  = d[v]()()() * phi[ss+v]()(3)(1)  + u[v]()()()*hp_11; 
	Simd p_32  = d[v]()()() * phi[ss+v]()(3)(2)  + u[v]()()()*hp_12; 

	
	//	if ( ss==0){
	/*
	std::cout << ss<<" "<< v<< " good "<< chi[ss+v]()(0)(0) << " bad "<<p_00<<" diff "<<chi[ss+v]()(0)(0)-p_00<<std::endl;
	std::cout << ss<<" "<< v<< " good "<< chi[ss+v]()(0)(1) << " bad "<<p_01<<" diff "<<chi[ss+v]()(0)(1)-p_01<<std::endl;
	std::cout << ss<<" "<< v<< " good "<< chi[ss+v]()(0)(2) << " bad "<<p_02<<" diff "<<chi[ss+v]()(0)(2)-p_02<<std::endl;
	std::cout << ss<<" "<< v<< " good "<< chi[ss+v]()(1)(0) << " bad "<<p_10<<" diff "<<chi[ss+v]()(1)(0)-p_10<<std::endl;
	std::cout << ss<<" "<< v<< " good "<< chi[ss+v]()(1)(1) << " bad "<<p_11<<" diff "<<chi[ss+v]()(1)(1)-p_11<<std::endl;
	std::cout << ss<<" "<< v<< " good "<< chi[ss+v]()(1)(2) << " bad "<<p_12<<" diff "<<chi[ss+v]()(1)(2)-p_12<<std::endl;
	std::cout << ss<<" "<< v<< " good "<< chi[ss+v]()(2)(0) << " bad "<<p_20<<" diff "<<chi[ss+v]()(2)(0)-p_20<<std::endl;
	std::cout << ss<<" "<< v<< " good "<< chi[ss+v]()(2)(1) << " bad "<<p_21<<" diff "<<chi[ss+v]()(2)(1)-p_21<<std::endl;
	std::cout << ss<<" "<< v<< " good "<< chi[ss+v]()(2)(2) << " bad "<<p_22<<" diff "<<chi[ss+v]()(2)(2)-p_22<<std::endl;
	std::cout << ss<<" "<< v<< " good "<< chi[ss+v]()(3)(0) << " bad "<<p_30<<" diff "<<chi[ss+v]()(3)(0)-p_30<<std::endl;
	std::cout << ss<<" "<< v<< " good "<< chi[ss+v]()(3)(1) << " bad "<<p_31<<" diff "<<chi[ss+v]()(3)(1)-p_31<<std::endl;
	std::cout << ss<<" "<< v<< " good "<< chi[ss+v]()(3)(2) << " bad "<<p_32<<" diff "<<chi[ss+v]()(3)(2)-p_32<<std::endl;
	}
	*/
	vstream(chi[ss+v]()(0)(0),p_00);
	vstream(chi[ss+v]()(0)(1),p_01);
	vstream(chi[ss+v]()(0)(2),p_02);
	vstream(chi[ss+v]()(1)(0),p_10);
	vstream(chi[ss+v]()(1)(1),p_11);
	vstream(chi[ss+v]()(1)(2),p_12);
	vstream(chi[ss+v]()(2)(0),p_20);
	vstream(chi[ss+v]()(2)(1),p_21);
	vstream(chi[ss+v]()(2)(2),p_22);
	vstream(chi[ss+v]()(3)(0),p_30);
	vstream(chi[ss+v]()(3)(1),p_31);
	vstream(chi[ss+v]()(3)(2),p_32);

      }
#endif
  }
  M5Dtime+=usecond();
}

template<class Impl>  
void CayleyFermion5D<Impl>::M5Ddag(const FermionField &psi,
				   const FermionField &phi, 
				   FermionField &chi,
				   std::vector<Coeff_t> &lower,
				   std::vector<Coeff_t> &diag,
				   std::vector<Coeff_t> &upper)
{
  GridBase *grid=psi._grid;
  int Ls   = this->Ls;
  int LLs  = grid->_rdimensions[0];
  int nsimd= Simd::Nsimd();

  Vector<iSinglet<Simd> > u(LLs);
  Vector<iSinglet<Simd> > l(LLs);
  Vector<iSinglet<Simd> > d(LLs);

  assert(Ls/LLs==nsimd);
  assert(phi.checkerboard == psi.checkerboard);

  chi.checkerboard=psi.checkerboard;

  // just directly address via type pun
  typedef typename Simd::scalar_type scalar_type;
  scalar_type * u_p = (scalar_type *)&u[0];
  scalar_type * l_p = (scalar_type *)&l[0];
  scalar_type * d_p = (scalar_type *)&d[0];

  for(int o=0;o<LLs;o++){ // outer
  for(int i=0;i<nsimd;i++){ //inner
    int s  = o+i*LLs;
    int ss = o*nsimd+i;
    u_p[ss] = upper[s];
    l_p[ss] = lower[s];
    d_p[ss] = diag[s];
  }}

  M5Dcalls++;
  M5Dtime-=usecond();
PARALLEL_FOR_LOOP
  for(int ss=0;ss<grid->oSites();ss+=LLs){ // adds LLs

    alignas(64) SiteHalfSpinor hp;
    alignas(64) SiteHalfSpinor hm;
    alignas(64) SiteSpinor fp;
    alignas(64) SiteSpinor fm;

    for(int v=0;v<LLs;v++){

      int vp=(v+1)%LLs;
      int vm=(v+LLs-1)%LLs;

      spProj5p(hp,psi[ss+vp]);
      spProj5m(hm,psi[ss+vm]);

      if ( vp<=v ) rotate(hp,hp,1);
      if ( vm>=v ) rotate(hm,hm,nsimd-1);
      
      hp=hp*0.5;
      hm=hm*0.5;
      spRecon5p(fp,hp);
      spRecon5m(fm,hm);

      chi[ss+v] = d[v]*phi[ss+v]+u[v]*fp;
      chi[ss+v] = chi[ss+v]     +l[v]*fm;

    }
  }
  M5Dtime+=usecond();
}
template<class Impl>
void CayleyFermion5D<Impl>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv)
{
  int Ls=this->Ls;
  int LLs = psi._grid->_rdimensions[0];
  int vol = psi._grid->oSites()/LLs;

  chi.checkerboard=psi.checkerboard;
  
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
  Vector<iSinglet<Simd> > Matp(Ls*LLs);
  Vector<iSinglet<Simd> > Matm(Ls*LLs);

  for(int s2=0;s2<Ls;s2++){
  for(int s1=0;s1<LLs;s1++){
    int istride = LLs;
    int ostride = 1;
      Simd Vp;
      Simd Vm;
      scalar_type *sp = (scalar_type *)&Vp;
      scalar_type *sm = (scalar_type *)&Vm;
      for(int l=0;l<Nsimd;l++){
	sp[l] = PplusMat (l*istride+s1*ostride ,s2);
	sm[l] = PminusMat(l*istride+s1*ostride,s2);
      }
      Matp[LLs*s2+s1] = Vp;
      Matm[LLs*s2+s1] = Vm;
    }
  }
  
  MooeeInvCalls++;
  MooeeInvTime-=usecond();
  // Dynamic allocate on stack to get per thread without serialised heap acces
#pragma omp parallel  
  {

    Vector<SiteHalfSpinor> SitePplus(LLs);
    Vector<SiteHalfSpinor> SitePminus(LLs);
    Vector<SiteHalfSpinor> SiteChiP(LLs);
    Vector<SiteHalfSpinor> SiteChiM(LLs);
    Vector<SiteSpinor>     SiteChi(LLs);

    SiteHalfSpinor BcastP;
    SiteHalfSpinor BcastM;

#pragma omp for 
  for(auto site=0;site<vol;site++){

    for(int s=0;s<LLs;s++){
      int lex = s+LLs*site;
      spProj5p(SitePplus[s] ,psi[lex]);
      spProj5m(SitePminus[s],psi[lex]);
      SiteChiP[s]=zero;
      SiteChiM[s]=zero;
    }
      
    int s=0;
    for(int  l=0; l<Simd::Nsimd();l++){ // simd lane
      for(int s2=0;s2<LLs;s2++){ // Column loop of right hand side
	vbroadcast(BcastP,SitePplus [s2],l);
	vbroadcast(BcastM,SitePminus[s2],l);
	for(int s1=0;s1<LLs;s1++){ // Column loop of reduction variables
	  SiteChiP[s1]=SiteChiP[s1]+Matp[LLs*s+s1]*BcastP;
	  SiteChiM[s1]=SiteChiM[s1]+Matm[LLs*s+s1]*BcastM;
	}
      s++;
    }}

    for(int s=0;s<LLs;s++){
      int lex = s+LLs*site;
      spRecon5p(SiteChi[s],SiteChiP[s]);
      accumRecon5m(SiteChi[s],SiteChiM[s]);
      chi[lex] = SiteChi[s]*0.5;
    }
  }
  }
  MooeeInvTime+=usecond();
}

INSTANTIATE_DPERP(DomainWallVec5dImplD);
INSTANTIATE_DPERP(DomainWallVec5dImplF);
INSTANTIATE_DPERP(ZDomainWallVec5dImplD);
INSTANTIATE_DPERP(ZDomainWallVec5dImplF);

template void CayleyFermion5D<DomainWallVec5dImplF>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<DomainWallVec5dImplD>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<ZDomainWallVec5dImplF>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<ZDomainWallVec5dImplD>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);

}}
