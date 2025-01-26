#pragma once
namespace Grid {

class Band
{
  RealD lo, hi;
public:
  Band(RealD _lo,RealD _hi)
  {
    lo=_lo;
    hi=_hi;
  }
  RealD operator() (RealD x){
    if ( x>lo && x<hi ){
      return 1.0;
    } else {
      return 0.0;
    }
  }
};

class PowerSpectrum
{ 
 public: 

  template<typename T>  static RealD normalise(T& v) 
  {
    RealD nn = norm2(v);
    nn = sqrt(nn);
    v = v * (1.0/nn);
    return nn;
  }

  std::vector<RealD> ranges;
  std::vector<int> order;
  
  PowerSpectrum(  std::vector<RealD> &bins, std::vector<int> &_order ) : ranges(bins), order(_order)  { };

  template<class Field>
  RealD operator()(LinearOperatorBase<Field> &HermOp, const Field &src) 
  { 
    GridBase *grid = src.Grid(); 
    int N=ranges.size();
    RealD hi = ranges[N-1];

    RealD lo_band = 0.0;
    RealD hi_band;
    RealD nn=norm2(src);
    RealD ss=0.0;

    Field tmp = src;

    for(int b=0;b<N;b++){
      hi_band = ranges[b];
      Band Notch(lo_band,hi_band);
      
      Chebyshev<Field> polynomial;
      polynomial.Init(0.0,hi,order[b],Notch);
      polynomial.JacksonSmooth();

      polynomial(HermOp,src,tmp) ;

      RealD p=norm2(tmp);
      ss=ss+p;
      std::cout << GridLogMessage << " PowerSpectrum Band["<<lo_band<<","<<hi_band<<"] power "<<norm2(tmp)/nn<<std::endl;
      
      lo_band=hi_band;
    }
    std::cout << GridLogMessage << " PowerSpectrum total power "<<ss/nn<<std::endl;
    std::cout << GridLogMessage << " PowerSpectrum total power (unnormalised) "<<nn<<std::endl;

    return 0;
  };
};
  
}
