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
#include <Grid.h>

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
  template<class vec> void operator()(vec &rr,vec &i1,vec &i2) const { rr = conjugate(i1);}
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
//  mac,mult,sub,add, vone,vzero,vcomplex_i, =zero,
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
  std::cout<<GridLogMessage << " " << func.name()<<std::endl;

  int ok=0;
  for(int i=0;i<Nsimd;i++){
    if ( abs(reference[i]-result[i])>0){
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


template<class reduced,class scal, class vec,class functor > 
void ReductionTester(const functor &func)
{
  GridSerialRNG          sRNG;
  sRNG.SeedRandomDevice();
  
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
  Tester<RealF,vRealF>(funcAdj());
  Tester<RealF,vRealF>(funcConj());
  Tester<RealF,vRealF>(funcInnerProduct());
  ReductionTester<RealF,RealF,vRealF>(funcReduce());

  std::cout << GridLogMessage <<"==================================="<<  std::endl;
  std::cout << GridLogMessage <<"Testing vRealD "<<std::endl;
  std::cout << GridLogMessage <<"==================================="<<  std::endl;

  Tester<RealD,vRealD>(funcPlus());
  Tester<RealD,vRealD>(funcMinus());
  Tester<RealD,vRealD>(funcTimes());
  Tester<RealD,vRealD>(funcAdj());
  Tester<RealD,vRealD>(funcConj());
  Tester<RealD,vRealD>(funcInnerProduct());
  ReductionTester<RealD,RealD,vRealD>(funcReduce());

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
  Tester<ComplexF,vComplexF>(funcInnerProduct());
  ReductionTester<ComplexF,ComplexF,vComplexF>(funcReduce());

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
  Tester<ComplexD,vComplexD>(funcInnerProduct());
  ReductionTester<ComplexD,ComplexD,vComplexD>(funcReduce());

  Grid_finalize();
}
