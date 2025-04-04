#pragma once
namespace Grid {
template<class Field> class PowerMethod  
{ 
 public: 

  template<typename T>  static RealD normalise(T& v) 
  {
    RealD nn = norm2(v);
    nn = sqrt(nn);
    v = v * (1.0/nn);
    return nn;
  }

  RealD operator()(LinearOperatorBase<Field> &HermOp, const Field &src) 
  { 
    GridBase *grid = src.Grid(); 
    
    // quickly get an idea of the largest eigenvalue to more properly normalize the residuum 
    RealD evalMaxApprox = 0.0; 
    auto src_n = src; 
    auto tmp = src; 
    const int _MAX_ITER_EST_ = 200; 

    for (int i=0;i<_MAX_ITER_EST_;i++) { 
      
      normalise(src_n); 
      HermOp.HermOp(src_n,tmp); 
      RealD vnum = real(innerProduct(src_n,tmp)); // HermOp. 
      RealD vden = norm2(src_n); 
      RealD na = vnum/vden; 

      std::cout << GridLogMessage << "PowerMethod: Current approximation of largest eigenvalue " << na << std::endl;
      
      //      if ( (fabs(evalMaxApprox/na - 1.0) < 0.0001) || (i==_MAX_ITER_EST_-1) ) { 
	// 	evalMaxApprox = na; 
	// 	return evalMaxApprox; 
      //      } 
      evalMaxApprox = na; 
      src_n = tmp;
    }
    std::cout << GridLogMessage << " Approximation of largest eigenvalue: " << evalMaxApprox << std::endl;
    return evalMaxApprox;
  }
};
}
