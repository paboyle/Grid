    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./gcc-bug-report/broken.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
  std::cout<<GridLogMessage<<"Leaf "<<obj::NestLevel<<std::endl;
  return arg;
}
template<int N,class obj,typename std::enable_if<N!=obj::NestLevel >::type * = nullptr > auto function(const obj &arg)-> obj
{
  std::cout<<GridLogMessage<<"Node "<<obj::NestLevel<<std::endl;
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
  ret = function<2>(array);
}
