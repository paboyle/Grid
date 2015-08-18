#ifndef GRID_QCD_ACTION_PARAMS_H
#define GRID_QCD_ACTION_PARAMS_H

namespace Grid {
namespace QCD {

    // These can move into a params header and be given MacroMagic serialisation
    struct GparityWilsonImplParams {
      std::vector<int> twists; 
    };

    struct WilsonImplParams { };

    struct OneFlavourRationalParams { 
      RealD  lo;
      RealD  hi;
      int precision=64;
      int    degree=10;
      RealD tolerance; // Vector? 
      RealD MaxIter;   // Vector?
      OneFlavourRationalParams (RealD lo,RealD hi,int precision=64,int degree = 10);
    };

}}

#endif
