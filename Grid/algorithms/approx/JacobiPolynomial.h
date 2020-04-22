#ifndef GRID_JACOBIPOLYNOMIAL_H
#define GRID_JACOBIPOLYNOMIAL_H

#include <Grid/algorithms/LinearOperator.h>

NAMESPACE_BEGIN(Grid);

template<class Field>
class JacobiPolynomial : public OperatorFunction<Field> {
 private:
  using OperatorFunction<Field>::operator();

  int order;
  RealD hi;
  RealD lo;
  RealD alpha;
  RealD beta;

 public:
  void csv(std::ostream &out){
    csv(out,lo,hi);
  }
  void csv(std::ostream &out,RealD llo,RealD hhi){
    RealD diff = hhi-llo;
    RealD delta = diff*1.0e-5;
    for (RealD x=llo-delta; x<=hhi; x+=delta) {
      RealD f = approx(x);
      out<< x<<" "<<f <<std::endl;
    }
    return;
  }

  JacobiPolynomial(){};
  JacobiPolynomial(RealD _lo,RealD _hi,int _order,RealD _alpha, RealD _beta)
  {
      lo=_lo;
      hi=_hi;
      alpha=_alpha;
      beta=_beta;
      order=_order;
  };

  RealD approx(RealD x) // Convenience for plotting the approximation                                                       
  {
    RealD Tn;
    RealD Tnm;
    RealD Tnp;

    RealD y=( x-0.5*(hi+lo))/(0.5*(hi-lo));

    RealD T0=1.0;
    RealD T1=(alpha-beta)*0.5+(alpha+beta+2.0)*0.5*y;

    Tn =T1;
    Tnm=T0;
    for(int n=2;n<=order;n++){
      RealD cnp = 2.0*n*(n+alpha+beta)*(2.0*n-2.0+alpha+beta);
      RealD cny = (2.0*n-2.0+alpha+beta)*(2.0*n-1.0+alpha+beta)*(2.0*n+alpha+beta);
      RealD cn1 = (2.0*n+alpha+beta-1.0)*(alpha*alpha-beta*beta);
      RealD cnm = - 2.0*(n+alpha-1.0)*(n+beta-1.0)*(2.0*n+alpha+beta);
      Tnp= ( cny * y *Tn + cn1 * Tn + cnm * Tnm )/ cnp;
      Tnm=Tn;
      Tn =Tnp;
    }
    return Tnp;
  };

  // Implement the required interface                                                                                       
  void operator() (LinearOperatorBase<Field> &Linop, const Field &in, Field &out) {
    GridBase *grid=in.Grid();

    int vol=grid->gSites();

    Field T0(grid);
    Field T1(grid);
    Field T2(grid);
    Field y(grid);


    Field *Tnm = &T0;
    Field *Tn  = &T1;
    Field *Tnp = &T2;

    //    RealD T0=1.0;                                                                                                     
    T0=in;

    //    RealD y=( x-0.5*(hi+lo))/(0.5*(hi-lo));                                                                           
    //           = x * 2/(hi-lo) - (hi+lo)/(hi-lo)                                                                          
    Linop.HermOp(T0,y);
    RealD xscale = 2.0/(hi-lo);
    RealD mscale = -(hi+lo)/(hi-lo);
    Linop.HermOp(T0,y);
    y=y*xscale+in*mscale;

    // RealD T1=(alpha-beta)*0.5+(alpha+beta+2.0)*0.5*y;
    RealD halfAmB  = (alpha-beta)*0.5;
    RealD halfApBp2= (alpha+beta+2.0)*0.5;
    T1 = halfAmB * in + halfApBp2*y;

    for(int n=2;n<=order;n++){

      Linop.HermOp(*Tn,y);
      y=xscale*y+mscale*(*Tn);

      RealD cnp = 2.0*n*(n+alpha+beta)*(2.0*n-2.0+alpha+beta);
      RealD cny = (2.0*n-2.0+alpha+beta)*(2.0*n-1.0+alpha+beta)*(2.0*n+alpha+beta);
      RealD cn1 = (2.0*n+alpha+beta-1.0)*(alpha*alpha-beta*beta);
      RealD cnm = - 2.0*(n+alpha-1.0)*(n+beta-1.0)*(2.0*n+alpha+beta);

      //      Tnp= ( cny * y *Tn + cn1 * Tn + cnm * Tnm )/ cnp;                                                             
      cny=cny/cnp;
      cn1=cn1/cnp;
      cn1=cn1/cnp;
      cnm=cnm/cnp;

      *Tnp=cny*y + cn1 *(*Tn) + cnm * (*Tnm);

      // Cycle pointers to avoid copies                                                                                     
      Field *swizzle = Tnm;
      Tnm    =Tn;
      Tn     =Tnp;
      Tnp    =swizzle;
    }
    out=*Tnp;

  }
};
NAMESPACE_END(Grid);
#endif
