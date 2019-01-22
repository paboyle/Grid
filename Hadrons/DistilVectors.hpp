/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/A2AVectors.hpp

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>
Author: fionnoh <fionnoh@gmail.com>

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
#ifndef Distil_Vectors_hpp_
#define Distil_Vectors_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Environment.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

template<typename LatticeObj>
class Perambulator : Serializable{
	  // TODO: The next line makes friends across all combinations
	  //     (not much of a problem given all public anyway ...)
	  //      FYI, the bug here was that I forgot that the friend is templated
  template<typename T> friend std::ostream & operator<<(std::ostream &os, const Perambulator<T>& p);
protected:
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS( Perambulator,
                                  std::string,             ID,          // Allows owner to specialise
                                  std::string,             Provenance,  // For info only
                                  std::vector<int>,        dimensions,
                                  std::vector<LatticeObj>, perambulator,
                                  // Following items are redundant, but useful
                                  int,     nd,          // Number of dimensions
                                  size_t,  NumElements); // Number of elements
protected:
  // Constructor common code
  inline void ConstructCommon(const int * Dimensions) {
    assert(nd > 0);
    dimensions.resize(nd);
    NumElements = 1;
    for(int i = 0 ; i < nd ; i++) {
      assert(Dimensions[i] > 0);
      NumElements *= (size_t) Dimensions[i];
      dimensions[i] = Dimensions[i];
    }
    //const LatticeObj perambulatorDefault;
    perambulator.resize(NumElements);//,perambulatorDefault);
  }
public:
  // Constructor with dimensions passed as std::vector<int>
  inline Perambulator(const std::vector<int> & Dimensions)
  : nd {(int) Dimensions.size()} {
    ConstructCommon( &Dimensions[0] ); }
  
  // Constructor with dimensions passed as std::vector<int>
  inline Perambulator(const std::vector<int> & Dimensions, const std::string sID)
  : nd {(int) Dimensions.size()}, ID(sID) {
    ConstructCommon( &Dimensions[0] ); }
  
  // Constructor with dimensions passed as std::vector<int>
  inline Perambulator(const std::vector<int> & Dimensions, const std::string sID, const std::string sProvenance)
  : nd {(int) Dimensions.size()}, ID(sID), Provenance(sProvenance) {
    ConstructCommon( &Dimensions[0] ); }
  
  // Constructor with dimensions passed as individual parameters
  // FYI: The caller is free to ignore the names and use the indices however they see fit
  inline Perambulator(int NumNoise, int NumEvec=1, int NumTime=1, int NumSpin=1, int I_k=1, int I_t=1, int I_s=1) {
    int Dimensions[]={NumNoise,NumEvec,NumTime,NumSpin,I_k,I_t,I_s};
    nd = sizeof(Dimensions)/sizeof(Dimensions[0]);
    while( nd > 1 && Dimensions[nd-1] == 1 )
      nd--;
    ConstructCommon( Dimensions );
  }
  
  inline LatticeObj & operator()(size_t count, const int * Coord) {
    assert( count == nd );
    assert( Coord );
    size_t idx = 0;
    // C memory order (???)
    for( int d = 0 ; d < nd ; d++ ) {
      assert( Coord[d] < dimensions[d] );
      idx *= (size_t) dimensions[d];
      idx += (size_t) Coord[d];
    }
    return perambulator[idx];
  }

  inline LatticeObj & operator()(const std::vector<int> Coord) {
    return operator()(Coord.size(), &Coord[0]);
  }
  
  inline LatticeObj & operator()(int idxNoise, int idxEvec=0, int idxTime=0, int idxSpin=0, int I_k=0, int I_t=0, int I_s=0) {
    int MyIndex[]={idxNoise,idxEvec,idxTime,idxSpin,I_k,I_t,I_s};
    int i = sizeof(MyIndex)/sizeof(MyIndex[0]);
    assert( i >= nd );
    while( i > nd )
      assert(MyIndex[--i] == 0);
    return operator()(i, MyIndex);
  }
};


END_HADRONS_NAMESPACE

#endif // Distil_Vectors_hpp_
