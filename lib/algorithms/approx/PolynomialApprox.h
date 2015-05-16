#ifndef GRID_POLYNOMIAL_APPROX_H
#define GRID_POLYNOMIAL_APPROX_H

#include<Grid.h>
#include<algorithms/LinearOperator.h>

namespace Grid {

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Simple general polynomial with user supplied coefficients
  ////////////////////////////////////////////////////////////////////////////////////////////
  template<class Field>
  class Polynomial : public OperatorFunction<Field> {
  private:
    std::vector<double> _oeffs;
  public:
    Polynomial(std::vector<double> &_Coeffs) : Coeffs(_Coeffs) {};

    // Implement the required interface
    void operator() (LinearOperatorBase<Field> &Linop, const Field &in, Field &out) {

      Field AtoN = in;
      out = AtoN*Coeffs[0];

      for(int n=1;n<Coeffs.size();n++){
	Field Mtmp=AtoN;
	Linop.Op(Mtmp,AtoN);
	out=out+AtoN*Coeffs[n];
      }
    };
  };

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Generic Chebyshev approximations
  ////////////////////////////////////////////////////////////////////////////////////////////
  template<class Field>
  class Chebyshev : public OperatorFunction<Field> {
  private:
    std::vector<double> Coeffs;
    int order;
    double hi;
    double lo;

  public:
    Chebyshev(double _lo,double _hi,int _order, double (* func)(double) ){
      lo=_lo;
      hi=_hi;
      order=_order;
      
      if(order < 2) exit(-1);
      Coeffs.resize(order);
      for(int j=0;j<order;j++){
	double s=0;
	for(int k=0;k<order;k++){
	  double y=cos(M_PI*(k+0.5)/order);
	  double x=0.5*(y*(hi-lo)+(hi+lo));
	  double f=func(x);
	  s=s+f*cos( j*M_PI*(k+0.5)/order );
	}
	Coeffs[j] = s * 2.0/order;
      }
    };

    double Evaluate(double x) // Convenience for plotting the approximation
    {
      double Tn;
      double Tnm;
      double Tnp;
      
      double y=( x-0.5*(hi+lo))/(0.5*(hi-lo));
      
      double T0=1;
      double T1=y;
      
      double sum;
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

    // Convenience for plotting the approximation
    void   PlotApprox(std::ostream &out) {
      out<<"Polynomial approx ["<<lo<<","<<hi<<"]"<<std::endl;
      for(double x=lo;x<hi;x+=(hi-lo)/50.0){
	out <<x<<"\t"<<Evaluate(x)<<std::endl;
      }
    };

    // Implement the required interface; could require Lattice base class
    void operator() (LinearOperatorBase<Field> &Linop, const Field &in, Field &out) {

      Field T0 = in;
      Field T1 = T0; // Field T1(T0._grid); more efficient but hardwires Lattice class
      Field T2 = T1;
      
      // use a pointer trick to eliminate copies
      Field *Tnm = &T0;
      Field *Tn  = &T1;
      Field *Tnp = &T2;
      Field y   = in;
  
      double xscale = 2.0/(hi-lo);
      double mscale = -(hi+lo)/(hi-lo);

      Field *T0=Tnm;
      Field *T1=Tn;
      

      // Tn=T1 = (xscale M + mscale)in
      Linop.Op(T0,y);

      T1=y*xscale+in*mscale;

      // sum = .5 c[0] T0 + c[1] T1
      out = (0.5*coeffs[0])*T0 + coeffs[1]*T1;

      for(int n=2;n<order;n++){
	
	Linop.Op(*Tn,y);

	y=xscale*y+mscale*(*Tn);

	*Tnp=2.0*y-(*Tnm);
	
	out=out+coeffs[n]* (*Tnp);

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
