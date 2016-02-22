    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/approx/Chebyshev.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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
#ifndef GRID_CHEBYSHEV_H
#define GRID_CHEBYSHEV_H

#include<Grid.h>
#include<algorithms/LinearOperator.h>

namespace Grid {

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Simple general polynomial with user supplied coefficients
  ////////////////////////////////////////////////////////////////////////////////////////////
  template<class Field>
  class HermOpOperatorFunction : public OperatorFunction<Field> {
    void operator() (LinearOperatorBase<Field> &Linop, const Field &in, Field &out) {
      Linop.HermOp(in,out);
    };
  };

  template<class Field>
  class Polynomial : public OperatorFunction<Field> {
  private:
    std::vector<RealD> Coeffs;
  public:
    Polynomial(std::vector<RealD> &_Coeffs) : Coeffs(_Coeffs) { };

    // Implement the required interface
    void operator() (LinearOperatorBase<Field> &Linop, const Field &in, Field &out) {

      Field AtoN(in._grid);
      Field Mtmp(in._grid);
      AtoN = in;
      out = AtoN*Coeffs[0];
      //      std::cout <<"Poly in " <<norm2(in)<<std::endl;
      //      std::cout <<"0 " <<norm2(out)<<std::endl;
      for(int n=1;n<Coeffs.size();n++){
	Mtmp = AtoN;
	Linop.HermOp(Mtmp,AtoN);
	out=out+AtoN*Coeffs[n];
	//	std::cout << n<<" " <<norm2(out)<<std::endl;
      }
    };
  };

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Generic Chebyshev approximations
  ////////////////////////////////////////////////////////////////////////////////////////////
  template<class Field>
  class Chebyshev : public OperatorFunction<Field> {
  private:
    std::vector<RealD> Coeffs;
    int order;
    RealD hi;
    RealD lo;

  public:
    void csv(std::ostream &out){
      for (RealD x=lo; x<hi; x+=(hi-lo)/1000) {
	RealD f = approx(x);
	out<< x<<" "<<f<<std::endl;
      }
      return;
    }

    // Convenience for plotting the approximation
    void   PlotApprox(std::ostream &out) {
      out<<"Polynomial approx ["<<lo<<","<<hi<<"]"<<std::endl;
      for(RealD x=lo;x<hi;x+=(hi-lo)/50.0){
	out <<x<<"\t"<<approx(x)<<std::endl;
      }
    };

    Chebyshev(){};
    Chebyshev(RealD _lo,RealD _hi,int _order, RealD (* func)(RealD) ) {Init(_lo,_hi,_order,func);};
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // c.f. numerical recipes "chebft"/"chebev". This is sec 5.8 "Chebyshev approximation".
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    void Init(RealD _lo,RealD _hi,int _order, RealD (* func)(RealD))
    {
      lo=_lo;
      hi=_hi;
      order=_order;
      
      if(order < 2) exit(-1);
      Coeffs.resize(order);
      for(int j=0;j<order;j++){
	RealD s=0;
	for(int k=0;k<order;k++){
	  RealD y=std::cos(M_PI*(k+0.5)/order);
	  RealD x=0.5*(y*(hi-lo)+(hi+lo));
	  RealD f=func(x);
	  s=s+f*std::cos( j*M_PI*(k+0.5)/order );
	}
	Coeffs[j] = s * 2.0/order;
      }
    };

    
    void JacksonSmooth(void){
      RealD M=order;
      RealD alpha = M_PI/(M+2);
      RealD lmax = std::cos(alpha);
      RealD sumUsq =0;
      std::vector<RealD> U(M);
      std::vector<RealD> a(M);
      std::vector<RealD> g(M);
      for(int n=0;n<=M;n++){
	U[n] = std::sin((n+1)*std::acos(lmax))/std::sin(std::acos(lmax));
	sumUsq += U[n]*U[n];
      }      
      sumUsq = std::sqrt(sumUsq);

      for(int i=1;i<=M;i++){
	a[i] = U[i]/sumUsq;
      }
      g[0] = 1.0;
      for(int m=1;m<=M;m++){
	g[m] = 0;
	for(int i=0;i<=M-m;i++){
	  g[m]+= a[i]*a[m+i];
	}
      }
      for(int m=1;m<=M;m++){
	Coeffs[m]*=g[m];
      }
    }
    RealD approx(RealD x) // Convenience for plotting the approximation
    {
      RealD Tn;
      RealD Tnm;
      RealD Tnp;
      
      RealD y=( x-0.5*(hi+lo))/(0.5*(hi-lo));
      
      RealD T0=1;
      RealD T1=y;
      
      RealD sum;
      sum = 0.5*Coeffs[0]*T0;
      sum+= Coeffs[1]*T1;
      
      Tn =T1;
      Tnm=T0;
      for(int i=2;i<order;i++){
	Tnp=2*y*Tn-Tnm;
	Tnm=Tn;
	Tn =Tnp;
	sum+= Tn*Coeffs[i];
      }
      return sum;
    };

    // Implement the required interface
    void operator() (LinearOperatorBase<Field> &Linop, const Field &in, Field &out) {

      GridBase *grid=in._grid;

      int vol=grid->gSites();

      Field T0(grid); T0 = in;  
      Field T1(grid); 
      Field T2(grid);
      Field y(grid);
      
      Field *Tnm = &T0;
      Field *Tn  = &T1;
      Field *Tnp = &T2;

      // Tn=T1 = (xscale M + mscale)in
      RealD xscale = 2.0/(hi-lo);
      RealD mscale = -(hi+lo)/(hi-lo);
      Linop.HermOp(T0,y);
      T1=y*xscale+in*mscale;

      // sum = .5 c[0] T0 + c[1] T1
      out = (0.5*Coeffs[0])*T0 + Coeffs[1]*T1;
      for(int n=2;n<order;n++){
	
	Linop.HermOp(*Tn,y);

	y=xscale*y+mscale*(*Tn);

	*Tnp=2.0*y-(*Tnm);

	out=out+Coeffs[n]* (*Tnp);

	// Cycle pointers to avoid copies
	Field *swizzle = Tnm;
	Tnm    =Tn;
	Tn     =Tnp;
	Tnp    =swizzle;
	  
      }
    }
  };


  template<class Field>
  class ChebyshevLanczos : public Chebyshev<Field> {
  private:
    std::vector<RealD> Coeffs;
    int order;
    RealD alpha;
    RealD beta;
    RealD mu;

  public:
    ChebyshevLanczos(RealD _alpha,RealD _beta,RealD _mu,int _order) :
    alpha(_alpha),
      beta(_beta),
          mu(_mu)
    {
      order=_order;
      Coeffs.resize(order);
      for(int i=0;i<_order;i++){
	Coeffs[i] = 0.0;
      }
      Coeffs[order-1]=1.0;
    };

    void csv(std::ostream &out){
      for (RealD x=-1.2*alpha; x<1.2*alpha; x+=(2.0*alpha)/10000) {
	RealD f = approx(x);
	out<< x<<" "<<f<<std::endl;
      }
      return;
    }

    RealD approx(RealD xx) // Convenience for plotting the approximation
    {
      RealD Tn;
      RealD Tnm;
      RealD Tnp;
      Real aa = alpha * alpha;
      Real bb = beta  *  beta;
      
      RealD x = ( 2.0 * (xx-mu)*(xx-mu) - (aa+bb) ) / (aa-bb);

      RealD y= x;
      
      RealD T0=1;
      RealD T1=y;
      
      RealD sum;
      sum = 0.5*Coeffs[0]*T0;
      sum+= Coeffs[1]*T1;
      
      Tn =T1;
      Tnm=T0;
      for(int i=2;i<order;i++){
	Tnp=2*y*Tn-Tnm;
	Tnm=Tn;
	Tn =Tnp;
	sum+= Tn*Coeffs[i];
      }
      return sum;
    };

    // shift_Multiply in Rudy's code
    void AminusMuSq(LinearOperatorBase<Field> &Linop, const Field &in, Field &out) 
    {
      GridBase *grid=in._grid;
      Field tmp(grid);

      RealD aa= alpha*alpha;
      RealD bb= beta * beta;

      Linop.HermOp(in,out);
      out = out - mu*in;

      Linop.HermOp(out,tmp);
      tmp = tmp - mu * out;

      out = (2.0/ (aa-bb) ) * tmp -  ((aa+bb)/(aa-bb))*in;
    };
    // Implement the required interface
    void operator() (LinearOperatorBase<Field> &Linop, const Field &in, Field &out) {

      GridBase *grid=in._grid;

      int vol=grid->gSites();

      Field T0(grid); T0 = in;  
      Field T1(grid); 
      Field T2(grid);
      Field  y(grid);
      
      Field *Tnm = &T0;
      Field *Tn  = &T1;
      Field *Tnp = &T2;

      // Tn=T1 = (xscale M )*in
      AminusMuSq(Linop,T0,T1);

      // sum = .5 c[0] T0 + c[1] T1
      out = (0.5*Coeffs[0])*T0 + Coeffs[1]*T1;
      for(int n=2;n<order;n++){
	
	AminusMuSq(Linop,*Tn,y);

	*Tnp=2.0*y-(*Tnm);

	out=out+Coeffs[n]* (*Tnp);

	// Cycle pointers to avoid copies
	Field *swizzle = Tnm;
	Tnm    =Tn;
	Tn     =Tnp;
	Tnp    =swizzle;
	  
      }
    }
  };
}
#endif
