#ifndef GRID_LATTICE_BASE_H
#define GRID_LATTICE_BASE_H

namespace Grid {

// TODO: 
//       mac,real,imag

// Functionality:
//     -=,+=,*=,()
//     add,+,sub,-,mult,mac,*
//     adj,conj
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

class LatticeBase {};
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

template<class vobj>
class Lattice : public LatticeBase
{
public:

    GridBase *_grid;
    int checkerboard;
    std::vector<vobj,alignedAllocator<vobj> > _odata;

public:
    typedef typename vobj::scalar_type scalar_type;
    typedef typename vobj::vector_type vector_type;
    typedef vobj vector_object;
    
  ////////////////////////////////////////////////////////////////////////////////
  // Expression Template closure support
  ////////////////////////////////////////////////////////////////////////////////
  template <typename Op, typename T1>                         inline Lattice<vobj> & operator=(const LatticeUnaryExpression<Op,T1> &expr)
  {
    //PARALLEL_FOR_LOOP
#pragma omp parallel for
    for(int ss=0;ss<_grid->oSites();ss++){
      vobj tmp= eval(ss,expr);
      vstream(_odata[ss] ,tmp);
    }
    return *this;
  }
  template <typename Op, typename T1,typename T2>             inline Lattice<vobj> & operator=(const LatticeBinaryExpression<Op,T1,T2> &expr)
  {
    // PARALLEL_FOR_LOOP
#pragma omp parallel for
    for(int ss=0;ss<_grid->oSites();ss++){
      vobj tmp= eval(ss,expr);
      vstream(_odata[ss] ,tmp);
    }
    return *this;
  }
  template <typename Op, typename T1,typename T2,typename T3> inline Lattice<vobj> & operator=(const LatticeTrinaryExpression<Op,T1,T2,T3> &expr)
  {
    //PARALLEL_FOR_LOOP
#pragma omp parallel for
    for(int ss=0;ss<_grid->oSites();ss++){
      vobj tmp= eval(ss,expr);
      vstream(_odata[ss] ,tmp);
    }
    return *this;
  }
  //GridFromExpression is tricky to do
  template<class Op,class T1>
    Lattice(const LatticeUnaryExpression<Op,T1> & expr):    _grid(nullptr){
    GridFromExpression(_grid,expr);
    assert(_grid!=nullptr);
    _odata.resize(_grid->oSites());
PARALLEL_FOR_LOOP
    for(int ss=0;ss<_grid->oSites();ss++){
      _odata[ss] = eval(ss,expr);
    }
  };
  template<class Op,class T1, class T2>
  Lattice(const LatticeBinaryExpression<Op,T1,T2> & expr):    _grid(nullptr){
    GridFromExpression(_grid,expr);
    assert(_grid!=nullptr);
    _odata.resize(_grid->oSites());
PARALLEL_FOR_LOOP
    for(int ss=0;ss<_grid->oSites();ss++){
      _odata[ss] = eval(ss,expr);
    }
  };
  template<class Op,class T1, class T2, class T3>
  Lattice(const LatticeTrinaryExpression<Op,T1,T2,T3> & expr):    _grid(nullptr){
    GridFromExpression(_grid,expr);
    assert(_grid!=nullptr);
    _odata.resize(_grid->oSites());
PARALLEL_FOR_LOOP
    for(int ss=0;ss<_grid->oSites();ss++){
      _odata[ss] = eval(ss,expr);
    }
  };

    //////////////////////////////////////////////////////////////////
    // Constructor requires "grid" passed.
    // what about a default grid?
    //////////////////////////////////////////////////////////////////
 Lattice(GridBase *grid) : _grid(grid), _odata(_grid->oSites()) {
      //        _odata.reserve(_grid->oSites());
      //        _odata.resize(_grid->oSites());
        assert((((uint64_t)&_odata[0])&0xF) ==0);
        checkerboard=0;
    }

    template<class sobj> inline Lattice<vobj> & operator = (const sobj & r){
PARALLEL_FOR_LOOP
        for(int ss=0;ss<_grid->oSites();ss++){
            this->_odata[ss]=r;
        }
        return *this;
    }
    template<class robj> inline Lattice<vobj> & operator = (const Lattice<robj> & r){
      conformable(*this,r);
      std::cout<<"Lattice operator ="<<std::endl;
PARALLEL_FOR_LOOP
        for(int ss=0;ss<_grid->oSites();ss++){
            this->_odata[ss]=r._odata[ss];
        }
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
    
    inline friend Lattice<vobj> operator / (const Lattice<vobj> &lhs,const Lattice<vobj> &rhs){
        conformable(lhs,rhs);
        Lattice<vobj> ret(lhs._grid);
PARALLEL_FOR_LOOP
        for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = lhs._odata[ss]/rhs._odata[ss];
        }
        return ret;
    };

 }; // class Lattice

  template<class vobj> inline std::ostream& operator<< (std::ostream& stream, const Lattice<vobj> &o){
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

}


#include <lattice/Grid_lattice_conformable.h>

#define GRID_LATTICE_EXPRESSION_TEMPLATES
#ifdef  GRID_LATTICE_EXPRESSION_TEMPLATES
#include <lattice/Grid_lattice_ET.h>
#else 
#include <lattice/Grid_lattice_overload.h>
#endif

#include <lattice/Grid_lattice_arith.h>

#include <lattice/Grid_lattice_trace.h>
#include <lattice/Grid_lattice_transpose.h>
#include <lattice/Grid_lattice_local.h>
#include <lattice/Grid_lattice_reduction.h>
#include <lattice/Grid_lattice_peekpoke.h>
#include <lattice/Grid_lattice_reality.h>
#include <Grid_extract.h>
#include <lattice/Grid_lattice_coordinate.h>
#include <lattice/Grid_lattice_rng.h>
#include <lattice/Grid_lattice_transfer.h>



#endif
