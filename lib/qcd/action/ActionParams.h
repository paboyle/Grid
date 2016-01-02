#ifndef GRID_QCD_ACTION_PARAMS_H
#define GRID_QCD_ACTION_PARAMS_H

namespace Grid {
namespace QCD {

    // These can move into a params header and be given MacroMagic serialisation
    struct GparityWilsonImplParams {
      std::vector<int> twists; 
      GparityWilsonImplParams () : twists(Nd,0) {};

    };

    struct WilsonImplParams { };

    struct OneFlavourRationalParams { 
      RealD  lo;
      RealD  hi;
      int MaxIter;   // Vector?
      RealD tolerance; // Vector? 
      int    degree=10;
      int precision=64;

      OneFlavourRationalParams (RealD _lo,RealD _hi,int _maxit,RealD tol=1.0e-8,int _degree = 10,int _precision=64) :
        lo(_lo), hi(_hi), MaxIter(_maxit), tolerance(tol), degree(_degree), precision(_precision)
      {};
    };

}}

#endif
