#include <Grid.h>
#include "simd/Grid_vector_types.h"
#include <parallelIO/GridNerscIO.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

class funcPlus {
public:
  funcPlus() {};
  template<class vec> void operator()(vec &rr,vec &i1,vec &i2) const { rr = i1+i2;}
  std::string name(void) const { return std::string("Plus"); }
};
class funcMinus {
public:
  funcMinus() {};
  template<class vec> void operator()(vec &rr,vec &i1,vec &i2) const { rr = i1-i2;}
  std::string name(void) const { return std::string("Minus"); }
};
class funcTimes {
public:
  funcTimes() {};
  template<class vec> void operator()(vec &rr,vec &i1,vec &i2) const { rr = i1*i2;}
  std::string name(void) const { return std::string("Times"); }
};
class funcConj {
public:
  funcConj() {};
  template<class vec> void operator()(vec &rr,vec &i1,vec &i2) const { rr = conj(i1);}
  std::string name(void) const { return std::string("Conj"); }
};
class funcAdj {
public:
  funcAdj() {};
  template<class vec> void operator()(vec &rr,vec &i1,vec &i2) const { rr = adj(i1);}
  std::string name(void) const { return std::string("Adj"); }
};

class funcTimesI {
public:
  funcTimesI() {};
  template<class vec> void operator()(vec &rr,vec &i1,vec &i2) const { rr = timesI(i1);}
  std::string name(void) const { return std::string("timesI"); }
};

class funcTimesMinusI {
public:
  funcTimesMinusI() {};
  template<class vec> void operator()(vec &rr,vec &i1,vec &i2) const { rr = timesMinusI(i1);}
  std::string name(void) const { return std::string("timesMinusI"); }
};

template<class scal, class vec,class functor > 
void Tester(const functor &func)
{
  GridSerialRNG          sRNG;
  sRNG.SeedRandomDevice();
  
  int Nsimd = vec::Nsimd();

  std::vector<scal> input1(Nsimd);
  std::vector<scal> input2(Nsimd);
  std::vector<scal> result(Nsimd);
  std::vector<scal> reference(Nsimd);

  std::vector<vec,alignedAllocator<vec> > buf(3);
  vec & v_input1 = buf[0];
  vec & v_input2 = buf[1];
  vec & v_result = buf[2];


  for(int i=0;i<Nsimd;i++){
    random(sRNG,input1[i]);
    random(sRNG,input2[i]);
    random(sRNG,result[i]);
  }

  merge<vec,scal>(v_input1,input1);
  merge<vec,scal>(v_input2,input2);
  merge<vec,scal>(v_result,result);

  func(v_result,v_input1,v_input2);

  for(int i=0;i<Nsimd;i++) {
    func(reference[i],input1[i],input2[i]);
  }

  extract<vec,scal>(v_result,result);
  std::cout << " " << func.name()<<std::endl;

  int ok=0;
  for(int i=0;i<Nsimd;i++){
    if ( abs(reference[i]-result[i])>0){
      std::cout<< "*****" << std::endl;
      std::cout<< "["<<i<<"] "<< abs(reference[i]-result[i]) << " " <<reference[i]<< " " << result[i]<<std::endl;
      ok++;
    }
  }
  if ( ok==0 ) std::cout << " OK!" <<std::endl;

}



int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(4,MyComplexF::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
    
  GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
  std::vector<int> seeds({1,2,3,4});

  // Insist that operations on random scalars gives
  // identical results to on vectors.

  std::cout << "==================================="<<  std::endl;
  std::cout << "Testing MyComplexF "<<std::endl;
  std::cout << "==================================="<<  std::endl;

  Tester<ComplexF,MyComplexF>(funcTimesI());
  Tester<ComplexF,MyComplexF>(funcTimesMinusI());
  Tester<ComplexF,MyComplexF>(funcPlus());
  Tester<ComplexF,MyComplexF>(funcMinus());
  Tester<ComplexF,MyComplexF>(funcTimes());
  Tester<ComplexF,MyComplexF>(funcConj());
  Tester<ComplexF,MyComplexF>(funcAdj());

  std::cout << "==================================="<<  std::endl;
  std::cout << "Testing MyComplexD "<<std::endl;
  std::cout << "==================================="<<  std::endl;


  Tester<ComplexD,MyComplexD>(funcTimesI());
  Tester<ComplexD,MyComplexD>(funcTimesMinusI());
  Tester<ComplexD,MyComplexD>(funcPlus());
  Tester<ComplexD,MyComplexD>(funcMinus());
  Tester<ComplexD,MyComplexD>(funcTimes());
  Tester<ComplexD,MyComplexD>(funcConj());
  Tester<ComplexD,MyComplexD>(funcAdj());

  std::cout << "==================================="<<  std::endl;
  std::cout << "Testing MyRealF "<<std::endl;
  std::cout << "==================================="<<  std::endl;


  Tester<RealF,MyRealF>(funcPlus());
  Tester<RealF,MyRealF>(funcMinus());
  Tester<RealF,MyRealF>(funcTimes());
  Tester<RealF,MyRealF>(funcAdj());

  std::cout << "==================================="<<  std::endl;
  std::cout << "Testing MyRealD "<<std::endl;
  std::cout << "==================================="<<  std::endl;

  Tester<RealD,MyRealD>(funcPlus());
  Tester<RealD,MyRealD>(funcMinus());
  Tester<RealD,MyRealD>(funcTimes());
  Tester<RealD,MyRealD>(funcAdj());

  Grid_finalize();
}
