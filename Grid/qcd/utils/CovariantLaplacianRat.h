/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/scalar/CovariantLaplacianRat.h

Copyright (C) 2021

Author: Chulwoo Jung <chulwoo@bnl.gov>

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
			   /*  END LEGAL */
#pragma once 
#define MIXED_CG
//enable/disable push_back
#undef USE_CHRONO 

//#include <roctracer/roctx.h>

NAMESPACE_BEGIN(Grid);

struct LaplacianRatParams {

  RealD offset;
  int order;
  std::vector<RealD> a0;
  std::vector<RealD> a1;
  std::vector<RealD> b0;
  std::vector<RealD> b1;
  RealD b2; //for debugging
  int   MaxIter;
  RealD tolerance;
  int   precision;
  
  // constructor 
  LaplacianRatParams(int ord = 1,
                  int maxit     = 1000,
                  RealD tol     = 1.0e-8, 
                  int precision = 64)
    : offset(1.), order(ord),b2(1.),
      MaxIter(maxit),
      tolerance(tol),
      precision(precision){ 
      a0.resize(ord,0.);
      a1.resize(ord,0.);
      b0.resize(ord,0.);
      b1.resize(ord,0.);
      };
};



////////////////////////////////////////////////////////////
// Laplacian operator L on adjoint fields
//
// phi: adjoint field
// L: D_mu^dag D_mu
//
// L phi(x) = Sum_mu [ U_mu(x)phi(x+mu)U_mu(x)^dag + 
//                     U_mu(x-mu)^dag phi(x-mu)U_mu(x-mu)
//                     -2phi(x)]
//
// Operator designed to be encapsulated by
// an HermitianLinearOperator<.. , ..>
////////////////////////////////////////////////////////////

template <class Impl, class ImplF>
class LaplacianAdjointRat: public Metric<typename Impl::Field> {
  OperatorFunction<typename Impl::Field> &Solver;
  LaplacianRatParams Gparam;
  LaplacianRatParams Mparam;
  GridBase *grid;
  GridBase *grid_f;
  CovariantAdjointLaplacianStencil<Impl,typename Impl::LinkField> LapStencil;
  CovariantAdjointLaplacianStencil<ImplF,typename ImplF::LinkField> LapStencilF;
public:
  INHERIT_GIMPL_TYPES(Impl);
//   typedef typename GImpl::LinkField GaugeLinkField; \
//  typedef typename GImpl::Field GaugeField;         
  typedef typename ImplF::Field GaugeFieldF;
  typedef typename ImplF::LinkField GaugeLinkFieldF; \
  GaugeField Usav;
  GaugeFieldF UsavF;
  std::vector< std::vector<GaugeLinkField> > prev_solnsM;
  std::vector< std::vector<GaugeLinkField> > prev_solnsMinv;
  std::vector< std::vector<GaugeLinkField> > prev_solnsMDeriv;
  std::vector< std::vector<GaugeLinkField> > prev_solnsMinvDeriv;

	  LaplacianAdjointRat(GridBase* _grid, GridBase* _grid_f, OperatorFunction<GaugeField>& S, LaplacianRatParams& gpar, LaplacianRatParams& mpar)
    : grid(_grid),grid_f(_grid_f), LapStencil(_grid), LapStencilF(_grid_f), U(Nd, _grid), Solver(S), Gparam(gpar), Mparam(mpar),Usav(_grid), UsavF(_grid_f),
      prev_solnsM(4),prev_solnsMinv(4),prev_solnsMDeriv(4),prev_solnsMinvDeriv(4) {
//    std::cout<<GridLogMessage << "Generating degree "<<param.degree<<" for x^(1/2)"<<std::endl;
    this->triv=0;
        

  };
  LaplacianAdjointRat(){this->triv=0; printf("triv=%d\n",this->Trivial());}
  void Mdir(const GaugeField&, GaugeField&, int, int){ assert(0);}
  void MdirAll(const GaugeField&, std::vector<GaugeField> &){ assert(0);}
  void Mdiag(const GaugeField&, GaugeField&){ assert(0);}

  void ImportGauge(const GaugeField& _U) {
    RealD total=0.;
    for (int mu = 0; mu < Nd; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(_U, mu);
      total += norm2(U[mu]);
    }
    Usav = _U;
    precisionChange(UsavF,Usav);
    std::cout <<GridLogDebug << "ImportGauge:norm2(_U) = "<<" "<<total<<std::endl;
  }

  void MDerivLink(const GaugeLinkField& left, const GaugeLinkField& right,
              GaugeField& der) {
    RealD factor = -1. / (double(4 * Nd));

    for (int mu = 0; mu < Nd; mu++) {
      GaugeLinkField der_mu(der.Grid());
      der_mu = Zero();
//      for (int nu = 0; nu < Nd; nu++) {
//        GaugeLinkField left_nu = PeekIndex<LorentzIndex>(left, nu);
//        GaugeLinkField right_nu = PeekIndex<LorentzIndex>(right, nu);
        der_mu += U[mu] * Cshift(left, mu, 1) * adj(U[mu]) * right;
        der_mu += U[mu] * Cshift(right, mu, 1) * adj(U[mu]) * left;
//      }
      PokeIndex<LorentzIndex>(der, -factor * der_mu, mu);
    }
    std::cout << GridLogDebug <<"MDerivLink:  norm2(der) = "<<norm2(der)<<std::endl;
  }

  void MDerivInt(LaplacianRatParams &par, const GaugeField& left, const GaugeField& right,
              GaugeField& der ,  std::vector< std::vector<GaugeLinkField> >& prev_solns ) {

// get rid of this please
    std::cout<<GridLogMessage << "LaplaceStart " <<std::endl;
    RealD fac =  - 1. / (double(4 * Nd)) ;
    RealD coef=0.5;
    LapStencil.GaugeImport(Usav);
    LapStencilF.GaugeImport(UsavF);


    for (int nu=0;nu<Nd;nu++){
        GaugeLinkField right_nu = PeekIndex<LorentzIndex>(right, nu);
        GaugeLinkField left_nu = PeekIndex<LorentzIndex>(left, nu);
        GaugeLinkField LMinvMom(left.Grid());
    
        GaugeLinkField GMom(left.Grid());
        GaugeLinkField LMinvGMom(left.Grid());
    
        GaugeLinkField AGMom(left.Grid());
        GaugeLinkField MinvAGMom(left.Grid());
        GaugeLinkField LMinvAGMom(left.Grid());
    
        GaugeLinkField AMinvMom(left.Grid());
        GaugeLinkField LMinvAMom(left.Grid());
        GaugeLinkField temp(left.Grid());
        GaugeLinkField temp2(left.Grid());
    
        std::vector<GaugeLinkField> MinvMom(par.order,left.Grid());
    
        GaugeLinkField MinvGMom(left.Grid());
        GaugeLinkField Gtemp(left.Grid());
        GaugeLinkField Gtemp2(left.Grid());
    
    
        ConjugateGradient<GaugeLinkField> CG(par.tolerance,10000,false);
    //    ConjugateGradient<GaugeFieldF> CG_f(par.tolerance,10000,false);
        LaplacianParams LapPar(0.0001, 1.0, 10000, 1e-8, 12, 64);
    
        ChronoForecast< QuadLinearOperator<CovariantAdjointLaplacianStencil<Impl,GaugeLinkField>,GaugeLinkField> , GaugeLinkField> Forecast;
    
        GMom = par.offset * right_nu;
    
        for(int i =0;i<par.order;i++){
        QuadLinearOperator<CovariantAdjointLaplacianStencil<Impl,typename Impl::LinkField>,GaugeLinkField> QuadOp(LapStencil,par.b0[i],fac*par.b1[i],fac*fac*par.b2);
#if USE_CHRONO
        MinvMom[i] = Forecast(QuadOp, right_nu, prev_solns[nu]);
#endif
#ifndef MIXED_CG
        CG(QuadOp,right_nu,MinvMom[i]);
#else
        QuadLinearOperator<CovariantAdjointLaplacianStencil<ImplF,typename ImplF::LinkField>,GaugeLinkFieldF> QuadOpF(LapStencilF,par.b0[i],fac*par.b1[i],fac*fac*par.b2);
    //    QuadLinearOperator<LaplacianAdjointField<ImplF>,GaugeLinkFieldF> QuadOpF(LapStencilF,par.b0[i],par.b1[i],par.b2);
        MixedPrecisionConjugateGradient<GaugeLinkField,GaugeLinkFieldF> MixedCG(par.tolerance,10000,10000,grid_f,QuadOpF,QuadOp);
        MixedCG.InnerTolerance=par.tolerance;
        MixedCG(right_nu,MinvMom[i]);
    #endif
    #if USE_CHRONO
        prev_solns[nu].push_back(MinvMom[i]);
    #endif
        
        GMom += par.a0[i]*MinvMom[i]; 
        LapStencil.M(MinvMom[i],Gtemp2);
        GMom += par.a1[i]*fac*Gtemp2; 
        }
        for(int i =0;i<par.order;i++){
        QuadLinearOperator<CovariantAdjointLaplacianStencil<Impl,typename Impl::LinkField>,GaugeLinkField> QuadOp(LapStencil,par.b0[i],fac*par.b1[i],fac*fac*par.b2);
    
        MinvGMom = Forecast(QuadOp, GMom, prev_solns[nu]);
    #ifndef MIXED_CG
        CG(QuadOp,GMom,MinvGMom);
        LapStencil.M(MinvGMom, Gtemp2); LMinvGMom=fac*Gtemp2;
        CG(QuadOp,right_nu,MinvMom[i]);
    #else
        QuadLinearOperator<CovariantAdjointLaplacianStencil<ImplF,typename ImplF::LinkField>,GaugeLinkFieldF> QuadOpF(LapStencilF,par.b0[i],fac*par.b1[i],fac*fac*par.b2);
    //    QuadLinearOperator<LaplacianAdjointField<ImplF>,GaugeLinkFieldF> QuadOpF(LapStencilF,par.b0[i],par.b1[i],par.b2);
        MixedPrecisionConjugateGradient<GaugeLinkField,GaugeLinkFieldF> MixedCG(par.tolerance,10000,10000,grid_f,QuadOpF,QuadOp);
        MixedCG.InnerTolerance=par.tolerance;
        MixedCG(GMom,MinvGMom);
        LapStencil.M(MinvGMom, Gtemp2); LMinvGMom=fac*Gtemp2;
    //    Laplacian.M(MinvGMom, LMinvGMom);
        MixedCG(right_nu,MinvMom[i]);
    #endif
#if USE_CHRONO
        prev_solns[nu].push_back(MinvGMom);
#endif
    
        LapStencil.M(MinvMom[i], Gtemp2); LMinvMom=fac*Gtemp2;
        AMinvMom = par.a1[i]*LMinvMom;
        AMinvMom += par.a0[i]*MinvMom[i];
    
        LapStencil.M(AMinvMom, Gtemp2); LMinvAMom=fac*Gtemp2;
        LapStencil.M(MinvGMom, Gtemp2); temp=fac*Gtemp2;
        MinvAGMom = par.a1[i]*temp;
        MinvAGMom += par.a0[i]*MinvGMom;
        LapStencil.M(MinvAGMom, Gtemp2); LMinvAGMom=fac*Gtemp2;
    
    
        GaugeField tempDer(left.Grid());
        std::cout<<GridLogMessage << "force contraction "<< i <<std::endl;
    //    roctxRangePushA("RMHMC force contraction");
        MDerivLink(GMom,MinvMom[i],tempDer); der += coef*2*par.a1[i]*tempDer;
        MDerivLink(left_nu,MinvGMom,tempDer); der += coef*2*par.a1[i]*tempDer;
        MDerivLink(LMinvAGMom,MinvMom[i],tempDer); der += coef*-2.*par.b2*tempDer;
        MDerivLink(LMinvAMom,MinvGMom,tempDer); der += coef*-2.*par.b2*tempDer;
        MDerivLink(MinvAGMom,LMinvMom,tempDer); der += coef*-2.*par.b2*tempDer;
        MDerivLink(AMinvMom,LMinvGMom,tempDer); der += coef*-2.*par.b2*tempDer;
        MDerivLink(MinvAGMom,MinvMom[i],tempDer); der += coef*-2.*par.b1[i]*tempDer;
        MDerivLink(AMinvMom,MinvGMom,tempDer); der += coef*-2.*par.b1[i]*tempDer;
        std::cout<<GridLogMessage << "coef =  force contraction "<< i << "done "<< coef <<std::endl;
    //    roctxRangePop();
    
        }
    }
    std::cout<<GridLogMessage << "LaplaceEnd " <<std::endl;
//  exit(-42);
  }

  void MDeriv(const GaugeField& in, GaugeField& der) {
    MDeriv(in,in, der);
  }

  void MDeriv(const GaugeField& left, const GaugeField& right,
              GaugeField& der) {

    der=Zero();
    MDerivInt(Mparam, left, right, der,prev_solnsMDeriv );
    std::cout <<GridLogDebug << "MDeriv:norm2(der) = "<<norm2(der)<<std::endl;
  }

  void MinvDeriv(const GaugeField& in, GaugeField& der) {
    std::vector< std::vector<GaugeLinkField> > prev_solns(4);
    der=Zero();
    MDerivInt(Gparam, in, in, der,prev_solnsMinvDeriv);
    std::cout <<GridLogDebug << "MinvDeriv:norm2(der) = "<<norm2(der)<<std::endl;
  }


  void MSquareRootInt(LaplacianRatParams &par, GaugeField& P, std::vector< std::vector<GaugeLinkField> > & prev_solns ){

    std::cout<<GridLogMessage << "LaplaceStart " <<std::endl;
    RealD fac = -1. / (double(4 * Nd));
    LapStencil.GaugeImport(Usav);
    LapStencilF.GaugeImport(UsavF);
    for(int nu=0; nu<Nd;nu++){
        GaugeLinkField P_nu = PeekIndex<LorentzIndex>(P, nu);
        GaugeLinkField Gp(P.Grid());
        Gp = par.offset * P_nu;
        ConjugateGradient<GaugeLinkField> CG(par.tolerance,10000);
    //    ConjugateGradient<GaugeLinkFieldF> CG_f(1.0e-8,10000);
    
        ChronoForecast< QuadLinearOperator<CovariantAdjointLaplacianStencil<Impl,typename Impl::LinkField>,GaugeLinkField> , GaugeLinkField> Forecast;
    
        GaugeLinkField Gtemp(P.Grid());
        GaugeLinkField Gtemp2(P.Grid());
    
    
        for(int i =0;i<par.order;i++){
        QuadLinearOperator<CovariantAdjointLaplacianStencil<Impl,typename Impl::LinkField>,GaugeLinkField> QuadOp(LapStencil,par.b0[i],fac*par.b1[i],fac*fac*par.b2);
    
        Gtemp = Forecast(QuadOp, P_nu, prev_solns[nu]);
    #ifndef MIXED_CG
        CG(QuadOp,P_nu,Gtemp);
    #else
        QuadLinearOperator<CovariantAdjointLaplacianStencil<ImplF,typename ImplF::LinkField>,GaugeLinkFieldF> QuadOpF(LapStencilF,par.b0[i],fac*par.b1[i],fac*fac*par.b2);
    //    QuadLinearOperator<LaplacianAdjointField<ImplF>,GaugeFieldF> QuadOpF(LapStencilF,par.b0[i],par.b1[i],par.b2);
        MixedPrecisionConjugateGradient<GaugeLinkField,GaugeLinkFieldF> MixedCG(par.tolerance,10000,10000,grid_f,QuadOpF,QuadOp);
        MixedCG.InnerTolerance=par.tolerance;
        MixedCG(P_nu,Gtemp);
    #endif
    #if USE_CHRONO
        prev_solns[nu].push_back(Gtemp);
    #endif
    
        Gp += par.a0[i]*Gtemp; 
        LapStencil.M(Gtemp,Gtemp2);
        Gp += par.a1[i]*fac*Gtemp2; 
        }
        PokeIndex<LorentzIndex>(P, Gp, nu);
    }
    std::cout<<GridLogMessage << "LaplaceEnd " <<std::endl;
  }

  void MSquareRoot(GaugeField& P){
    std::vector< std::vector<GaugeLinkField> > prev_solns(4);
    MSquareRootInt(Mparam,P,prev_solns);
    std::cout <<GridLogDebug << "MSquareRoot:norm2(P) = "<<norm2(P)<<std::endl;
  }

  void MInvSquareRoot(GaugeField& P){
    std::vector< std::vector<GaugeLinkField> > prev_solns(4);
    MSquareRootInt(Gparam,P,prev_solns);
    std::cout <<GridLogDebug << "MInvSquareRoot:norm2(P) = "<<norm2(P)<<std::endl;
  }

  void M(const GaugeField& in, GaugeField& out) {
      out = in;
      std::vector< std::vector<GaugeLinkField> > prev_solns(4);
      MSquareRootInt(Mparam,out,prev_solns);
      MSquareRootInt(Mparam,out,prev_solns);
      std::cout <<GridLogDebug << "M:norm2(out) = "<<norm2(out)<<std::endl;
  }

  void Minv(const GaugeField& in, GaugeField& inverted){
      inverted = in;
      std::vector< std::vector<GaugeLinkField> > prev_solns(4);
      MSquareRootInt(Gparam,inverted,prev_solns);
      MSquareRootInt(Gparam,inverted,prev_solns);
      std::cout <<GridLogDebug << "Minv:norm2(inverted) = "<<norm2(inverted)<<std::endl;
  }



private:
  std::vector<GaugeLinkField> U;
};
#undef MIXED_CG

NAMESPACE_END(Grid);
