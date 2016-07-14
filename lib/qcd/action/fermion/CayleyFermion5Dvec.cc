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
				std::vector<RealD> &lower,
				std::vector<RealD> &diag,
				std::vector<RealD> &upper)
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

PARALLEL_FOR_LOOP
  for(int ss=0;ss<grid->oSites();ss+=LLs){ // adds LLs

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

      hp=hp*0.5;
      hm=hm*0.5;
      spRecon5m(fp,hp);
      spRecon5p(fm,hm);

      chi[ss+v] = d[v]*phi[ss+v]+u[v]*fp;
      chi[ss+v] = chi[ss+v]     +l[v]*fm;

    }
  }
}

template<class Impl>  
void CayleyFermion5D<Impl>::M5Ddag(const FermionField &psi,
				   const FermionField &phi, 
				   FermionField &chi,
				   std::vector<RealD> &lower,
				   std::vector<RealD> &diag,
				   std::vector<RealD> &upper)
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
}

template<class Impl>
void CayleyFermion5D<Impl>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv)
{
  int Ls=this->Ls;
  int LLs = psi._grid->_rdimensions[0];
  int vol = psi._grid->oSites()/LLs;

  chi.checkerboard=psi.checkerboard;
  
  Eigen::MatrixXd Pplus  = Eigen::MatrixXd::Zero(Ls,Ls);
  Eigen::MatrixXd Pminus = Eigen::MatrixXd::Zero(Ls,Ls);
  
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
  
  Eigen::MatrixXd PplusMat ;
  Eigen::MatrixXd PminusMat;
  
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
  
  // Dynamic allocate on stack to get per thread without serialised heap acces
PARALLEL_FOR_LOOP
  for(auto site=0;site<vol;site++){
    
    //    SiteHalfSpinor *SitePplus =(SiteHalfSpinor *) alloca(LLs*sizeof(SiteHalfSpinor));
    //    SiteHalfSpinor *SitePminus=(SiteHalfSpinor *) alloca(LLs*sizeof(SiteHalfSpinor));
    //    SiteSpinor     *SiteChi   =(SiteSpinor *)     alloca(LLs*sizeof(SiteSpinor));

    Vector<SiteHalfSpinor> SitePplus(LLs);
    Vector<SiteHalfSpinor> SitePminus(LLs);
    Vector<SiteHalfSpinor> SiteChiP(LLs);
    Vector<SiteHalfSpinor> SiteChiM(LLs);
    Vector<SiteSpinor>     SiteChi(LLs);

    SiteHalfSpinor BcastP;
    SiteHalfSpinor BcastM;

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

  FermOp5dVecTemplateInstantiate(CayleyFermion5D);

}}
