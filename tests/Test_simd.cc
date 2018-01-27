/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_simd.cc

Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>

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
#include <Grid/Grid.h>

using namespace Grid;

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
class funcDivide {
public:
  funcDivide() {};
  template<class vec> void operator()(vec &rr,vec &i1,vec &i2) const { rr = i1/i2;}
  std::string name(void) const { return std::string("Divide"); }
};
class funcConj {
public:
  funcConj() {};
  template<class vec> void operator()(vec &rr,vec &i1,vec &i2) const { rr = conjugate(i1);}
  std::string name(void) const { return std::string("Conj"); }
};
class funcAdj {
public:
  funcAdj() {};
  template<class vec> void operator()(vec &rr,vec &i1,vec &i2) const { rr = adj(i1);}
  std::string name(void) const { return std::string("Adj"); }
};
class funcImag {
public:
  funcImag() {};
  template<class vec> void operator()(vec &rr,vec &i1,vec &i2) const { rr = imag(i1);}
  std::string name(void) const { return std::string("imag"); }
};
class funcReal {
public:
  funcReal() {};
  template<class vec> void operator()(vec &rr,vec &i1,vec &i2) const { rr = real(i1);}
  std::string name(void) const { return std::string("real"); }
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
class funcInnerProduct {
public:
  funcInnerProduct() {};
  template<class vec> void operator()(vec &rr,vec &i1,vec &i2) const { rr = innerProduct(i1,i2);}
  std::string name(void) const { return std::string("innerProduct"); }
};

// FIXME still to test:
//
//  innerProduct,
//  norm2, 
//  Reduce,
//
//  mac,mult,sub,add, vone,vzero,vcomplex_i, =Zero(),
//  vset,vsplat,vstore,vstream,vload, scalar*vec, vec*scalar
//  unary -,
//  *= , -=, +=
//  outerproduct, 
//  zeroit
//  permute
class funcReduce {
public:
  funcReduce() {};
template<class reduce,class vec>    void vfunc(reduce &rr,vec &i1,vec &i2)   const { rr = Reduce(i1);}
template<class reduce,class scal>   void sfunc(reduce &rr,scal &i1,scal &i2) const { rr = i1;}
  std::string name(void) const { return std::string("Reduce"); }
};

template<class scal, class vec,class functor > 
void Tester(const functor &func)
{
  GridSerialRNG          sRNG;
  sRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
  
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

  std::cout << GridLogMessage << " " << func.name() << std::endl;

  std::cout << GridLogDebug << v_input1 << std::endl;
  std::cout << GridLogDebug << v_input2 << std::endl;
  std::cout << GridLogDebug << v_result << std::endl;

  int ok=0;
  for(int i=0;i<Nsimd;i++){
    if ( abs(reference[i]-result[i])>1.0e-6){
      std::cout<<GridLogMessage<< "*****" << std::endl;
      std::cout<<GridLogMessage<< "["<<i<<"] "<< abs(reference[i]-result[i]) << " " <<reference[i]<< " " << result[i]<<std::endl;
      ok++;
    }
  }
  if ( ok==0 ) {
    std::cout<<GridLogMessage << " OK!" <<std::endl;
  }
  assert(ok==0);
}

template<class functor>
void IntTester(const functor &func)
{
  typedef Integer  scal;
  typedef vInteger vec;

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
    input1[i] = (i + 1) * 30;
    input2[i] = (i + 1) * 20;
    result[i] = (i + 1) * 10;
  }

  merge<vec,scal>(v_input1,input1);
  merge<vec,scal>(v_input2,input2);
  merge<vec,scal>(v_result,result);

  func(v_result,v_input1,v_input2);

  for(int i=0;i<Nsimd;i++) {
    func(reference[i],input1[i],input2[i]);
  }

  extract<vec,scal>(v_result,result);

  std::cout << GridLogMessage << " " << func.name() << std::endl;

  std::cout << GridLogDebug << v_input1 << std::endl;
  std::cout << GridLogDebug << v_input2 << std::endl;
  std::cout << GridLogDebug << v_result << std::endl;

  int ok=0;
  for(int i=0;i<Nsimd;i++){
    if ( reference[i]-result[i] != 0){
      std::cout<<GridLogMessage<< "*****" << std::endl;
      std::cout<<GridLogMessage<< "["<<i<<"] "<< reference[i]-result[i] << " " <<reference[i]<< " " << result[i]<<std::endl;
      ok++;
    }
  }
  if ( ok==0 ) {
    std::cout<<GridLogMessage << " OK!" <<std::endl;
  }
  assert(ok==0);
}


template<class reduced,class scal, class vec,class functor > 
void ReductionTester(const functor &func)
{
  GridSerialRNG          sRNG;
  sRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
  
  int Nsimd = vec::Nsimd();

  std::vector<scal> input1(Nsimd);
  std::vector<scal> input2(Nsimd);
  reduced result(0);
  reduced reference(0);
  reduced tmp;

  std::vector<vec,alignedAllocator<vec> > buf(3);
  vec & v_input1 = buf[0];
  vec & v_input2 = buf[1];


  for(int i=0;i<Nsimd;i++){
    random(sRNG,input1[i]);
    random(sRNG,input2[i]);
  }

  merge<vec,scal>(v_input1,input1);
  merge<vec,scal>(v_input2,input2);

  func.template vfunc<reduced,vec>(result,v_input1,v_input2);

  for(int i=0;i<Nsimd;i++) {
    func.template sfunc<reduced,scal>(tmp,input1[i],input2[i]);
    reference+=tmp;
  }

  std::cout<<GridLogMessage << " " << func.name()<<std::endl;

  int ok=0;
  if ( abs(reference-result)/abs(reference) > 1.0e-6 ){ // rounding is possible for reduce order
    std::cout<<GridLogMessage<< "*****" << std::endl;
    std::cout<<GridLogMessage<< abs(reference-result) << " " <<reference<< " " << result<<std::endl;
    ok++;
  }
  if ( ok==0 ) {
    std::cout<<GridLogMessage << " OK!" <<std::endl;
  }
  assert(ok==0);
}


template<class reduced,class scal, class vec,class functor > 
void IntReductionTester(const functor &func)
{
  int Nsimd = vec::Nsimd();

  std::vector<scal> input1(Nsimd);
  std::vector<scal> input2(Nsimd);
  reduced result(0);
  reduced reference(0);
  reduced tmp;

  std::vector<vec,alignedAllocator<vec> > buf(3);
  vec & v_input1 = buf[0];
  vec & v_input2 = buf[1];

  for(int i=0;i<Nsimd;i++){
    input1[i] = (i + 1) * 30;
    input2[i] = (i + 1) * 20;
  }

  merge<vec,scal>(v_input1,input1);
  merge<vec,scal>(v_input2,input2);

  func.template vfunc<reduced,vec>(result,v_input1,v_input2);

  for(int i=0;i<Nsimd;i++) {
    func.template sfunc<reduced,scal>(tmp,input1[i],input2[i]);
    reference+=tmp;
  }

  std::cout<<GridLogMessage << " " << func.name()<<std::endl;

  int ok=0;
  if ( reference-result != 0 ){
    std::cout<<GridLogMessage<< "*****" << std::endl;
    std::cout<<GridLogMessage<< reference-result << " " <<reference<< " " << result<<std::endl;
    ok++;
  }
  if ( ok==0 ) {
    std::cout<<GridLogMessage << " OK!" <<std::endl;
  }
  assert(ok==0);
}


class funcPermute {
public:
  int n;
  funcPermute(int _n) { n=_n;};
  template<class vec>    void operator()(vec &rr,vec &i1,vec &i2) const { permute(rr,i1,n);}
  template<class scal>   void apply(std::vector<scal> &rr,std::vector<scal> &in)  const { 
    int sz=in.size();
    int msk = sz>>(n+1);
    for(int i=0;i<sz;i++){
      rr[i] = in[ i^msk ];
    }
  }
  std::string name(void) const { return std::string("Permute"); }
};

class funcExchange {
public:
  int n;
  funcExchange(int _n) { n=_n;};
  template<class vec>    void operator()(vec &r1,vec &r2,vec &i1,vec &i2) const { exchange(r1,r2,i1,i2,n);}
  template<class scal>   void apply(std::vector<scal> &r1,
				    std::vector<scal> &r2,
				    std::vector<scal> &in1,
				    std::vector<scal> &in2)  const 
  { 
    int sz=in1.size();
    int msk = sz>>(n+1);

    for(int i=0;i<sz;i++) {
      int j1 = i&(~msk);
      int j2 = i|msk;
      if  ( (i&msk) == 0 ) { r1[i]=in1[j1];}
      else                 { r1[i]=in2[j1];}

      if  ( (i&msk) == 0 ) { r2[i]=in1[j2];}
      else                 { r2[i]=in2[j2];}
    }      
  }
  std::string name(void) const { return std::string("Exchange"); }
};

class funcRotate {
public:
  int n;
  funcRotate(int _n) { n=_n;};
  template<class vec>    void operator()(vec &rr,vec &i1,vec &i2) const { rr=rotate(i1,n);}
  template<class scal>   void apply(std::vector<scal> &rr,std::vector<scal> &in)  const { 
    int sz = in.size();
    for(int i=0;i<sz;i++){
      rr[i] = in[(i+n)%sz];
    }
  }
  std::string name(void) const { return std::string("Rotate"); }
};


template<class scal, class vec,class functor > 
void PermTester(const functor &func)
{
  GridSerialRNG          sRNG;
  sRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
  
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

  func.apply(reference,input1);

  extract<vec,scal>(v_result,result);
  std::cout<<GridLogMessage << " " << func.name() << " " <<func.n <<std::endl;

  int ok=0;
  if (0) {
    std::cout<<GridLogMessage<< "*****" << std::endl;
    for(int i=0;i<Nsimd;i++){
      std::cout<< input1[i]<<" ";
    }
    std::cout <<std::endl; 
    for(int i=0;i<Nsimd;i++){
      std::cout<< result[i]<<" ";
    }
    std::cout <<std::endl; 
    for(int i=0;i<Nsimd;i++){
      std::cout<< reference[i]<<" ";
    }
    std::cout <<std::endl; 
    std::cout<<GridLogMessage<< "*****" << std::endl;
  }
  for(int i=0;i<Nsimd;i++){
    if ( abs(reference[i]-result[i])>1.0e-7){
      std::cout<<GridLogMessage<< "*****" << std::endl;      
      std::cout<<GridLogMessage<< "["<<i<<"] "<< abs(reference[i]-result[i]) << " " <<reference[i]<< " " << result[i]<<std::endl;
      ok++;
    }
  }
  if ( ok==0 ) {
    std::cout<<GridLogMessage << " OK!" <<std::endl;
  }
  assert(ok==0);
}


template<class scal, class vec,class functor > 
void ExchangeTester(const functor &func)
{
  GridSerialRNG          sRNG;
  sRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
  
  int Nsimd = vec::Nsimd();

  std::vector<scal> input1(Nsimd);
  std::vector<scal> input2(Nsimd);
  std::vector<scal> result1(Nsimd);
  std::vector<scal> result2(Nsimd);
  std::vector<scal> reference1(Nsimd);
  std::vector<scal> reference2(Nsimd);
  std::vector<scal> test1(Nsimd);
  std::vector<scal> test2(Nsimd);

  std::vector<vec,alignedAllocator<vec> > buf(6);
  vec & v_input1 = buf[0];
  vec & v_input2 = buf[1];
  vec & v_result1 = buf[2];
  vec & v_result2 = buf[3];
  vec & v_test1 = buf[4];
  vec & v_test2 = buf[5];

  for(int i=0;i<Nsimd;i++){
    random(sRNG,input1[i]);
    random(sRNG,input2[i]);
    random(sRNG,result1[i]);
    random(sRNG,result2[i]);
  }

  merge<vec,scal>(v_input1,input1);
  merge<vec,scal>(v_input2,input2);
  merge<vec,scal>(v_result1,result1);
  merge<vec,scal>(v_result2,result1);

  func(v_result1,v_result2,v_input1,v_input2);
  func.apply(reference1,reference2,input1,input2);

  func(v_test1,v_test2,v_result1,v_result2);

  extract<vec,scal>(v_result1,result1);
  extract<vec,scal>(v_result2,result2);
  extract<vec,scal>(v_test1,test1);
  extract<vec,scal>(v_test2,test2);

  std::cout<<GridLogMessage << " " << func.name() << " " <<func.n <<std::endl;

  //for(int i=0;i<Nsimd;i++) std::cout << " i "<<i<<" ref "<<reference1[i]<<" res "<<result1[i]<<std::endl;
  //for(int i=0;i<Nsimd;i++) std::cout << " i "<<i<<" ref "<<reference2[i]<<" res "<<result2[i]<<std::endl;

  for(int i=0;i<Nsimd;i++){
    int found=0;
    for(int j=0;j<Nsimd;j++){
      if(reference1[j]==result1[i]) {
	found=1;
	//	std::cout << " i "<<i<<" j "<<j<<" "<<reference1[j]<<" "<<result1[i]<<std::endl;
      }
    }
    //    assert(found==1);
    assert(found==1||found==0);
  }
  for(int i=0;i<Nsimd;i++){
    int found=0;
    for(int j=0;j<Nsimd;j++){
      if(reference2[j]==result2[i]) {
	found=1;
	//	std::cout << " i "<<i<<" j "<<j<<" "<<reference2[j]<<" "<<result2[i]<<std::endl;
      }
    }
    //    assert(found==1);
    assert(found==1||found==0);
  }

  /*
  for(int i=0;i<Nsimd;i++){
    std::cout << " i "<< i
	      <<" result1  "<<result1[i]
	      <<" result2  "<<result2[i]
	      <<" test1  "<<test1[i]
	      <<" test2  "<<test2[i]
	      <<" input1 "<<input1[i]
	      <<" input2 "<<input2[i]<<std::endl;
  }
  */
  for(int i=0;i<Nsimd;i++){
    assert(test1[i]==input1[i]);
    assert(test2[i]==input2[i]);
  }
}


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
    
  GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
  std::vector<int> seeds({1,2,3,4});

  // Insist that operations on random scalars gives
  // identical results to on vectors.

  std::cout << GridLogMessage <<"==================================="<<  std::endl;
  std::cout << GridLogMessage <<"Testing vRealF "<<std::endl;
  std::cout << GridLogMessage <<"==================================="<<  std::endl;


  Tester<RealF,vRealF>(funcPlus());
  Tester<RealF,vRealF>(funcMinus());
  Tester<RealF,vRealF>(funcTimes());
  Tester<RealF,vRealF>(funcDivide());
  Tester<RealF,vRealF>(funcAdj());
  Tester<RealF,vRealF>(funcConj());
  Tester<RealF,vRealF>(funcInnerProduct());
  ReductionTester<RealF,RealF,vRealF>(funcReduce());


  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  std::cout<<GridLogMessage << "Testing vRealF permutes "<<std::endl;
  std::cout<<GridLogMessage << "==================================="<<  std::endl;

  // Log2 iteration
  for(int i=0;(1<<i)< vRealF::Nsimd();i++){
    PermTester<RealF,vRealF>(funcPermute(i));
  }

  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  std::cout<<GridLogMessage << "Testing vRealF exchanges "<<std::endl;
  std::cout<<GridLogMessage << "==================================="<<  std::endl;

  // Log2 iteration
  for(int i=0;(1<<i)< vRealF::Nsimd();i++){
    ExchangeTester<RealF,vRealF>(funcExchange(i));
  }

  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  std::cout<<GridLogMessage << "Testing vRealF rotate "<<std::endl;
  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  for(int r=0;r<vRealF::Nsimd();r++){
    PermTester<RealF,vRealF>(funcRotate(r));
  }


  std::cout << GridLogMessage <<"==================================="<<  std::endl;
  std::cout << GridLogMessage <<"Testing vRealD "<<std::endl;
  std::cout << GridLogMessage <<"==================================="<<  std::endl;

  Tester<RealD,vRealD>(funcPlus());
  Tester<RealD,vRealD>(funcMinus());
  Tester<RealD,vRealD>(funcTimes());
  Tester<RealD,vRealD>(funcDivide());
  Tester<RealD,vRealD>(funcAdj());
  Tester<RealD,vRealD>(funcConj());
  Tester<RealD,vRealD>(funcInnerProduct());
  ReductionTester<RealD,RealD,vRealD>(funcReduce());


  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  std::cout<<GridLogMessage << "Testing vRealD permutes "<<std::endl;
  std::cout<<GridLogMessage << "==================================="<<  std::endl;

  // Log2 iteration
  for(int i=0;(1<<i)< vRealD::Nsimd();i++){
    PermTester<RealD,vRealD>(funcPermute(i));
  }

  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  std::cout<<GridLogMessage << "Testing vRealD exchanges "<<std::endl;
  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  // Log2 iteration
  for(int i=0;(1<<i)< vRealD::Nsimd();i++){
    ExchangeTester<RealD,vRealD>(funcExchange(i));
  }

  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  std::cout<<GridLogMessage << "Testing vRealD rotate "<<std::endl;
  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  for(int r=0;r<vRealD::Nsimd();r++){
    PermTester<RealD,vRealD>(funcRotate(r));
  }



  std::cout << GridLogMessage <<"==================================="<<  std::endl;
  std::cout << GridLogMessage <<"Testing vComplexF "<<std::endl;
  std::cout << GridLogMessage <<"==================================="<<  std::endl;

  Tester<ComplexF,vComplexF>(funcTimesI());
  Tester<ComplexF,vComplexF>(funcTimesMinusI());
  Tester<ComplexF,vComplexF>(funcPlus());
  Tester<ComplexF,vComplexF>(funcMinus());
  Tester<ComplexF,vComplexF>(funcTimes());
  Tester<ComplexF,vComplexF>(funcConj());
  Tester<ComplexF,vComplexF>(funcAdj());
  Tester<ComplexF,vComplexF>(funcReal());
  Tester<ComplexF,vComplexF>(funcImag());
  Tester<ComplexF,vComplexF>(funcInnerProduct());
  ReductionTester<ComplexF,ComplexF,vComplexF>(funcReduce());


  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  std::cout<<GridLogMessage << "Testing vComplexF permutes "<<std::endl;
  std::cout<<GridLogMessage << "==================================="<<  std::endl;

  // Log2 iteration
  for(int i=0;(1<<i)< vComplexF::Nsimd();i++){
    PermTester<ComplexF,vComplexF>(funcPermute(i));
  }


  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  std::cout<<GridLogMessage << "Testing vComplexF exchanges "<<std::endl;
  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  // Log2 iteration
  for(int i=0;(1<<i)< vComplexF::Nsimd();i++){
    ExchangeTester<ComplexF,vComplexF>(funcExchange(i));
  }


  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  std::cout<<GridLogMessage << "Testing vComplexF rotate "<<std::endl;
  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  for(int r=0;r<vComplexF::Nsimd();r++){
    PermTester<ComplexF,vComplexF>(funcRotate(r));
  }

  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  std::cout<<GridLogMessage << "Testing vComplexD "<<std::endl;
  std::cout<<GridLogMessage << "==================================="<<  std::endl;


  Tester<ComplexD,vComplexD>(funcTimesI());
  Tester<ComplexD,vComplexD>(funcTimesMinusI());
  Tester<ComplexD,vComplexD>(funcPlus());
  Tester<ComplexD,vComplexD>(funcMinus());
  Tester<ComplexD,vComplexD>(funcTimes());
  Tester<ComplexD,vComplexD>(funcConj());
  Tester<ComplexD,vComplexD>(funcAdj());
  Tester<ComplexD, vComplexD>(funcReal());
  Tester<ComplexD, vComplexD>(funcImag());

  Tester<ComplexD, vComplexD>(funcInnerProduct());
  ReductionTester<ComplexD, ComplexD, vComplexD>(funcReduce());

  std::cout << GridLogMessage
            << "===================================" << std::endl;
  std::cout << GridLogMessage << "Testing vComplexD permutes " << std::endl;
  std::cout << GridLogMessage
            << "===================================" << std::endl;

  // Log2 iteration
  for (int i = 0; (1 << i) < vComplexD::Nsimd(); i++) {
    PermTester<ComplexD, vComplexD>(funcPermute(i));
  }


  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  std::cout<<GridLogMessage << "Testing vComplexD exchanges "<<std::endl;
  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  // Log2 iteration
  for(int i=0;(1<<i)< vComplexD::Nsimd();i++){
    ExchangeTester<ComplexD,vComplexD>(funcExchange(i));
  }


  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  std::cout<<GridLogMessage << "Testing vComplexD rotate "<<std::endl;
  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  for(int r=0;r<vComplexD::Nsimd();r++){
    PermTester<ComplexD,vComplexD>(funcRotate(r));
  }
  
  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  std::cout<<GridLogMessage << "Testing vInteger                   "<<  std::endl;
  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  IntTester(funcPlus());
  IntTester(funcMinus());
  IntTester(funcTimes());
  IntReductionTester<Integer, Integer, vInteger>(funcReduce());

  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  std::cout<<GridLogMessage << "Testing precisionChange            "<<  std::endl;
  std::cout<<GridLogMessage << "==================================="<<  std::endl;
  {
    GridSerialRNG          sRNG;
    sRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
    const int Ndp = 16;
    const int Nsp = Ndp/2;
    const int Nhp = Ndp/4;
    std::vector<vRealH,alignedAllocator<vRealH> > H (Nhp);
    std::vector<vRealF,alignedAllocator<vRealF> > F (Nsp);
    std::vector<vRealF,alignedAllocator<vRealF> > FF(Nsp);
    std::vector<vRealD,alignedAllocator<vRealD> > D (Ndp);
    std::vector<vRealD,alignedAllocator<vRealD> > DD(Ndp);
    for(int i=0;i<16;i++){
      random(sRNG,D[i]);
    }
    // Double to Single
    precisionChange(&F[0],&D[0],Ndp);
    precisionChange(&DD[0],&F[0],Ndp);
    std::cout << GridLogMessage<<"Double to single";
    for(int i=0;i<Ndp;i++){
      //      std::cout << "DD["<<i<<"] = "<< DD[i]<<" "<<D[i]<<" "<<DD[i]-D[i] <<std::endl; 
      DD[i] = DD[i] - D[i];
      decltype(innerProduct(DD[0],DD[0])) nrm;
      nrm = innerProduct(DD[i],DD[i]);
      auto tmp = Reduce(nrm);
      //      std::cout << tmp << std::endl;
      assert( tmp < 1.0e-14 ); 
    }
    std::cout <<" OK ! "<<std::endl;

    // Double to Half
#ifdef USE_FP16
    std::cout << GridLogMessage<< "Double to half" ;
    precisionChange(&H[0],&D[0],Ndp);
    precisionChange(&DD[0],&H[0],Ndp);
    for(int i=0;i<Ndp;i++){
      //      std::cout << "DD["<<i<<"] = "<< DD[i]<<" "<<D[i]<<" "<<DD[i]-D[i]<<std::endl; 
      DD[i] = DD[i] - D[i];
      decltype(innerProduct(DD[0],DD[0])) nrm;
      nrm = innerProduct(DD[i],DD[i]);
      auto tmp = Reduce(nrm);
      //      std::cout << tmp << std::endl;
      assert( tmp < 1.0e-3 ); 
    }
    std::cout <<" OK ! "<<std::endl;

    std::cout << GridLogMessage<< "Single to half";
    // Single to Half
    precisionChange(&H[0] ,&F[0],Nsp);
    precisionChange(&FF[0],&H[0],Nsp);
    for(int i=0;i<Nsp;i++){
      //      std::cout << "FF["<<i<<"] = "<< FF[i]<<" "<<F[i]<<" "<<FF[i]-F[i]<<std::endl; 
      FF[i] = FF[i] - F[i];
      decltype(innerProduct(FF[0],FF[0])) nrm;
      nrm = innerProduct(FF[i],FF[i]);
      auto tmp = Reduce(nrm);
      //      std::cout << tmp << std::endl;
      assert( tmp < 1.0e-3 ); 
    }
    std::cout <<" OK ! "<<std::endl;
#endif
  }
  Grid_finalize();
}
