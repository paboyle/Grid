#include <vector>
#include <complex>
#include <type_traits>
#include <iostream>

typedef std::complex<double> ComplexD;

template <class T> class TypeMapper {
public:
  enum { NestLevel = T::NestLevel };
};

template<> class TypeMapper<ComplexD> {
public:
  enum { NestLevel = 0 };
};

template<class obj> class Container {
 public:
  std::vector<obj> data;
  Container(int size) : data (size){};
};

template<class obj> class Recursive {
public:
  enum { NestLevel = TypeMapper<obj>::NestLevel + 1};
  obj internal;
};

template<int N,class obj,typename std::enable_if<N==obj::NestLevel >::type * = nullptr > auto function(const obj &arg)-> obj
{
  std::cout<<"Leaf "<<obj::NestLevel<<std::endl;
  return arg;
}
template<int N,class obj,typename std::enable_if<N!=obj::NestLevel >::type * = nullptr > auto function(const obj &arg)-> obj
{
  std::cout<<"Node "<<obj::NestLevel<<std::endl;
  obj ret;
  ret.internal=function<N>(arg.internal);
  return ret;
}

template<int N,class obj> auto function(const Container<obj> & arg)-> Container<decltype(function<N>(arg.data[0]))>
{
  Container<decltype(function<N>(arg.data[0]))> ret(arg.data.size());
  for(int ss=0;ss<arg.data.size();ss++){
    ret.data[ss] = function<N>(arg.data[ss]);
  }
  return ret;
}


int main(int argc,char **argv)
{
  Container<Recursive<Recursive<ComplexD> > > array(10);
  Container<Recursive<Recursive<ComplexD> > > ret(10);
  
  ret = function<1>(array);
}
