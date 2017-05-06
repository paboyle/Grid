/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/gauge/WilsonGaugeAction.h

Copyright (C) 2016

Author: Guido Cossu <guido.cossu@ed.ac.uk>

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#ifndef HMC_MODULES_H
#define HMC_MODULES_H

/*
Define loadable, serializable modules
for the HMC execution
*/

namespace Grid {

// Empty class for no parameters
class NoParameters{};


/*
Base class for modules with parameters
*/
template < class P >
class Parametrized{
public:
  typedef P Parameters;

  Parametrized(Parameters Par):Par_(Par){};

  template <class ReaderClass>
  Parametrized(Reader<ReaderClass> & R, std::string section_name = "parameters"){
    read(R, section_name, Par_);
  }

  void set_parameters(Parameters Par){
        Par_ = Par;
  }

  void print_parameters(){
    std::cout << Par_ << std::endl;
  }

protected:
  Parameters Par_;
private:
  std::string section_name;
};


template <>
class Parametrized<NoParameters>{
        public:
  typedef NoParameters Parameters;

  Parametrized(Parameters Par){};

  template <class ReaderClass>
  Parametrized(Reader<ReaderClass> & Reader){};

  void set_parameters(Parameters Par){}

  void print_parameters(){}

};



////////////////////////////////////////
// Lowest level abstract module class
////////////////////////////////////////
template <class Prod>
class HMCModuleBase {
 public:
  typedef Prod Product;

  virtual Prod* getPtr() = 0;

  // add a getReference? 
  
  virtual void print_parameters(){};  // default to nothing
};


/////////////////////////////////////////////
// Registration class
/////////////////////////////////////////////

template <class T, class TheFactory>
class Registrar {
 public:
  Registrar(std::string className) {
    // register the class factory function
    TheFactory::getInstance().registerBuilder(className, 
        [&](typename TheFactory::TheReader Reader)
        { 
          return std::unique_ptr<T>(new T(Reader));
        }
        );
  }
};



}


#endif //HMC_MODULES_H