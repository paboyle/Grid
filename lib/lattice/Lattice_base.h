/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/lattice/Lattice_base.h

Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#pragma once 

#define STREAMING_STORES

NAMESPACE_BEGIN(Grid);

// TODO: 
//       mac,real,imag

// Functionality:
//     -=,+=,*=,()
//     add,+,sub,-,mult,mac,*
//     adj,conjugate
//     real,imag
//     transpose,transposeIndex  
//     trace,traceIndex
//     peekIndex
//     innerProduct,outerProduct,
//     localNorm2
//     localInnerProduct

extern int GridCshiftPermuteMap[4][16];

////////////////////////////////////////////////
// Basic expressions used in Expression Template
////////////////////////////////////////////////

class LatticeBase
{
public:
  virtual ~LatticeBase(void) = default;
  GridBase *_grid;
};
    
class LatticeExpressionBase {};

template <typename Op, typename T1>                           
class LatticeUnaryExpression  : public std::pair<Op,std::tuple<T1> > , public LatticeExpressionBase {
public:
  LatticeUnaryExpression(const std::pair<Op,std::tuple<T1> > &arg): std::pair<Op,std::tuple<T1> >(arg) {};
};

template <typename Op, typename T1, typename T2>              
class LatticeBinaryExpression : public std::pair<Op,std::tuple<T1,T2> > , public LatticeExpressionBase {
public:
  LatticeBinaryExpression(const std::pair<Op,std::tuple<T1,T2> > &arg): std::pair<Op,std::tuple<T1,T2> >(arg) {};
};

template <typename Op, typename T1, typename T2, typename T3> 
class LatticeTrinaryExpression :public std::pair<Op,std::tuple<T1,T2,T3> >, public LatticeExpressionBase {
public:
  LatticeTrinaryExpression(const std::pair<Op,std::tuple<T1,T2,T3> > &arg): std::pair<Op,std::tuple<T1,T2,T3> >(arg) {};
};

void inline conformable(GridBase *lhs,GridBase *rhs)
{
  assert(lhs == rhs);
}

template<class vobj>
class Lattice : public LatticeBase
{
private:
  //  vobj     * _odata_p;
  //  uint64_t   _odata_size;
public:
  int checkerboard;
  Vector<vobj> _odata;

  // to pthread need a computable loop where loop induction is not required
  accelerator int begin(void) const { return 0;};
  accelerator int end(void)   const { return _odata.size(); }
  accelerator_inline vobj & operator[](int i) { return _odata[i]; };
  accelerator_inline const vobj & operator[](int i) const { return _odata[i]; };

public:
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
  typedef vobj vector_object;

  void resize(uint64_t size)
  {
    _odata.resize(size);
    //    _odata_p    = &_odata[0];
    //    _odata_size = size;
  }   
  ////////////////////////////////////////////////////////////////////////////////
  // Expression Template closure support
  ////////////////////////////////////////////////////////////////////////////////
  template <typename Op, typename T1>                         inline Lattice<vobj> & operator=(const LatticeUnaryExpression<Op,T1> &expr)
  {
    GridBase *egrid(nullptr);
    GridFromExpression(egrid,expr);
    assert(egrid!=nullptr);
    conformable(_grid,egrid);

    int cb=-1;
    CBFromExpression(cb,expr);
    assert( (cb==Odd) || (cb==Even));
    checkerboard=cb;

#ifdef STREAMING_STORES
    accelerator_loop(ss,(*this),{
      vobj tmp = eval(ss,expr);
      vstream(_odata[ss] ,tmp);
    });
#else
    accelerator_loop(ss,(*this),{
      _odata[ss]=eval(ss,expr);
    });
#endif
    return *this;
  }
  template <typename Op, typename T1,typename T2> inline Lattice<vobj> & operator=(const LatticeBinaryExpression<Op,T1,T2> &expr)
  {
    GridBase *egrid(nullptr);
    GridFromExpression(egrid,expr);
    assert(egrid!=nullptr);
    conformable(_grid,egrid);

    int cb=-1;
    CBFromExpression(cb,expr);
    assert( (cb==Odd) || (cb==Even));
    checkerboard=cb;

#ifdef STREAMING_STORES
    accelerator_loop(ss,(*this),{
      vobj tmp = eval(ss,expr);
      vstream(_odata[ss] ,tmp);
    });
#else
    accelerator_loop(ss,(*this),{
      _odata[ss]=eval(ss,expr);
    });
#endif
    return *this;
  }
  template <typename Op, typename T1,typename T2,typename T3> inline Lattice<vobj> & operator=(const LatticeTrinaryExpression<Op,T1,T2,T3> &expr)
  {
    GridBase *egrid(nullptr);
    GridFromExpression(egrid,expr);
    assert(egrid!=nullptr);
    conformable(_grid,egrid);

    int cb=-1;
    CBFromExpression(cb,expr);
    assert( (cb==Odd) || (cb==Even));
    checkerboard=cb;

#ifdef STREAMING_STORES
    accelerator_loop(ss,(*this),{
      //vobj tmp = eval(ss,expr);
      vstream(_odata[ss] ,eval(ss,expr));
    });
#else
    accelerator_loop(ss,(*this),{
      _odata[ss] = eval(ss,expr);
    });
#endif
    return *this;
  }
  //GridFromExpression is tricky to do
  template<class Op,class T1>
  Lattice(const LatticeUnaryExpression<Op,T1> & expr) {
    _grid = nullptr;
    GridFromExpression(_grid,expr);
    assert(_grid!=nullptr);

    int cb=-1;
    CBFromExpression(cb,expr);
    assert( (cb==Odd) || (cb==Even));
    checkerboard=cb;

    resize(_grid->oSites());
#ifdef STREAMING_STORES
    cpu_loop(ss,(*this),{
      vobj tmp = eval(ss,expr);
      vstream(_odata[ss] ,tmp);
    });
#else
    cpu_loop(ss,(*this),{
      _odata[ss]=eval(ss,expr);
    });
#endif
  }
  template<class Op,class T1, class T2>
  Lattice(const LatticeBinaryExpression<Op,T1,T2> & expr) {
    _grid = nullptr;
    GridFromExpression(_grid,expr);
    assert(_grid!=nullptr);

    int cb=-1;
    CBFromExpression(cb,expr);
    assert( (cb==Odd) || (cb==Even));
    checkerboard=cb;
    
    resize(_grid->oSites());

#ifdef STREAMING_STORES
    cpu_loop(ss,(*this),{
      vobj tmp = eval(ss,expr);
      vstream(_odata[ss] ,tmp);
    });
#else
    cpu_loop(ss,(*this),{
      _odata[ss]=eval(ss,expr);
    });
#endif
  }
  template<class Op,class T1, class T2, class T3>
  Lattice(const LatticeTrinaryExpression<Op,T1,T2,T3> & expr) {
    _grid = nullptr;
    GridFromExpression(_grid,expr);
    assert(_grid!=nullptr);

    int cb=-1;
    CBFromExpression(cb,expr);
    assert( (cb==Odd) || (cb==Even));
    checkerboard=cb;

    resize(_grid->oSites());
    cpu_loop(ss,(*this),{
      vstream(_odata[ss] ,eval(ss,expr));
    });
  }

  //////////////////////////////////////////////////////////////////
  // Constructor requires "grid" passed.
  // what about a default grid?
  //////////////////////////////////////////////////////////////////
  Lattice(GridBase *grid) {
    _grid = grid;
    resize(grid->oSites());
    assert((((uint64_t)&_odata[0])&0xF) ==0);
    checkerboard=0;
  }
  
  Lattice(const Lattice& r){ // copy constructor
    _grid = r._grid;
    checkerboard = r.checkerboard;
    resize(_grid->oSites());// essential
    cpu_loop(ss,(*this),{
      _odata[ss]=r._odata[ss];
    });
  }
  
  
  
  virtual ~Lattice(void) = default;
    
  void reset(GridBase* grid) {
    if (_grid != grid) {
      _grid = grid;
      resize(grid->oSites());
      checkerboard = 0;
    }
  }
  

  template<class sobj> inline Lattice<vobj> & operator = (const sobj & r){
    accelerator_loop(ss,(*this),{
      this->_odata[ss]=r;
    });
    return *this;
  }
  
  template<class robj> inline Lattice<vobj> & operator = (const Lattice<robj> & r){
    this->checkerboard = r.checkerboard;
    conformable(*this,r);
    accelerator_loop(ss,(*this),{
      this->_odata[ss]=r._odata[ss];
    });
    return *this;
  }
  
  // *=,+=,-= operators inherit behvour from correspond */+/- operation
  template<class T> inline Lattice<vobj> &operator *=(const T &r) {
    *this = (*this)*r;
    return *this;
  }
  
  template<class T> inline Lattice<vobj> &operator -=(const T &r) {
    *this = (*this)-r;
    return *this;
  }
  template<class T> inline Lattice<vobj> &operator +=(const T &r) {
    *this = (*this)+r;
    return *this;
  }
}; // class Lattice
  
template<class vobj> std::ostream& operator<< (std::ostream& stream, const Lattice<vobj> &o){
  std::vector<int> gcoor;
  typedef typename vobj::scalar_object sobj;
  sobj ss;
  for(int g=0;g<o._grid->_gsites;g++){
    o._grid->GlobalIndexToGlobalCoor(g,gcoor);
    peekSite(ss,o,gcoor);
    stream<<"[";
    for(int d=0;d<gcoor.size();d++){
      stream<<gcoor[d];
      if(d!=gcoor.size()-1) stream<<",";
    }
    stream<<"]\t";
    stream<<ss<<std::endl;
  }
  return stream;
}
  
NAMESPACE_END(Grid);

