    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_ET.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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
#ifndef GRID_LATTICE_ET_H
#define GRID_LATTICE_ET_H

#include <iostream>
#include <vector>
#include <tuple>
#include <typeinfo>

namespace Grid {

  ////////////////////////////////////////////////////
  // Predicated where support
  ////////////////////////////////////////////////////
  template<class iobj,class vobj,class robj>
    inline vobj predicatedWhere(const iobj &predicate,const vobj &iftrue,const robj &iffalse) {

    typename std::remove_const<vobj>::type ret;

    typedef typename vobj::scalar_object scalar_object;
    typedef typename vobj::scalar_type scalar_type;
    typedef typename vobj::vector_type vector_type;

    const int Nsimd = vobj::vector_type::Nsimd();
    const int words = sizeof(vobj)/sizeof(vector_type);

    std::vector<Integer> mask(Nsimd);
    std::vector<scalar_object> truevals (Nsimd);
    std::vector<scalar_object> falsevals(Nsimd);

    extract(iftrue   ,truevals);
    extract(iffalse  ,falsevals);
    extract<vInteger,Integer>(TensorRemove(predicate),mask);

    for(int s=0;s<Nsimd;s++){
      if (mask[s]) falsevals[s]=truevals[s];
    }

    merge(ret,falsevals);
    return ret;
  }

////////////////////////////////////////////
// recursive evaluation of expressions; Could
// switch to generic approach with variadics, a la
// Antonin's Lat Sim but the repack to variadic with popped
// from tuple is hideous; C++14 introduces std::make_index_sequence for this
////////////////////////////////////////////


//leaf eval of lattice ; should enable if protect using traits

template <typename T> using is_lattice      = std::is_base_of<LatticeBase,T >;

template <typename T> using is_lattice_expr = std::is_base_of<LatticeExpressionBase,T >;

template<class sobj>
inline sobj eval(const unsigned int ss, const sobj &arg)
{
  return arg;
}
template<class lobj>
inline const lobj &eval(const unsigned int ss, const Lattice<lobj> &arg)
{
    return arg._odata[ss];
}

// handle nodes in syntax tree
template <typename Op, typename T1>
auto inline eval(const unsigned int ss, const LatticeUnaryExpression<Op,T1 > &expr) // eval one operand
  -> decltype(expr.first.func(eval(ss,std::get<0>(expr.second))))
{
  return expr.first.func(eval(ss,std::get<0>(expr.second)));
}

template <typename Op, typename T1, typename T2>
auto inline eval(const unsigned int ss, const LatticeBinaryExpression<Op,T1,T2> &expr) // eval two operands
  -> decltype(expr.first.func(eval(ss,std::get<0>(expr.second)),eval(ss,std::get<1>(expr.second))))
{
  return expr.first.func(eval(ss,std::get<0>(expr.second)),eval(ss,std::get<1>(expr.second)));
}

template <typename Op, typename T1, typename T2, typename T3>
auto inline eval(const unsigned int ss, const LatticeTrinaryExpression<Op,T1,T2,T3 > &expr) // eval three operands
  -> decltype(expr.first.func(eval(ss,std::get<0>(expr.second)),eval(ss,std::get<1>(expr.second)),eval(ss,std::get<2>(expr.second))))
{
  return expr.first.func(eval(ss,std::get<0>(expr.second)),eval(ss,std::get<1>(expr.second)),eval(ss,std::get<2>(expr.second)) );
}

//////////////////////////////////////////////////////////////////////////
// Obtain the grid from an expression, ensuring conformable. This must follow a tree recursion
//////////////////////////////////////////////////////////////////////////
template<class T1, typename std::enable_if<is_lattice<T1>::value, T1>::type * =nullptr >
inline void GridFromExpression(GridBase * &grid,const T1& lat)   // Lattice leaf
{
  if ( grid ) {
    conformable(grid,lat._grid);
  } 
  grid=lat._grid;
}
template<class T1,typename std::enable_if<!is_lattice<T1>::value, T1>::type * = nullptr >
inline void GridFromExpression(GridBase * &grid,const T1& notlat)   // non-lattice leaf
{
}
template <typename Op, typename T1>
inline void GridFromExpression(GridBase * &grid,const LatticeUnaryExpression<Op,T1 > &expr)
{
  GridFromExpression(grid,std::get<0>(expr.second));// recurse 
}

template <typename Op, typename T1, typename T2>
inline void GridFromExpression(GridBase * &grid,const LatticeBinaryExpression<Op,T1,T2> &expr) 
{
  GridFromExpression(grid,std::get<0>(expr.second));// recurse
  GridFromExpression(grid,std::get<1>(expr.second));
}
template <typename Op, typename T1, typename T2, typename T3>
inline void GridFromExpression( GridBase * &grid,const LatticeTrinaryExpression<Op,T1,T2,T3 > &expr) 
{
  GridFromExpression(grid,std::get<0>(expr.second));// recurse
  GridFromExpression(grid,std::get<1>(expr.second));
  GridFromExpression(grid,std::get<2>(expr.second));
}


//////////////////////////////////////////////////////////////////////////
// Obtain the CB from an expression, ensuring conformable. This must follow a tree recursion
//////////////////////////////////////////////////////////////////////////
template<class T1, typename std::enable_if<is_lattice<T1>::value, T1>::type * =nullptr >
inline void CBFromExpression(int &cb,const T1& lat)   // Lattice leaf
{
  if ( (cb==Odd) || (cb==Even) ) {
    assert(cb==lat.checkerboard);
  } 
  cb=lat.checkerboard;
  //  std::cout<<GridLogMessage<<"Lattice leaf cb "<<cb<<std::endl;
}
template<class T1,typename std::enable_if<!is_lattice<T1>::value, T1>::type * = nullptr >
inline void CBFromExpression(int &cb,const T1& notlat)   // non-lattice leaf
{
  //  std::cout<<GridLogMessage<<"Non lattice leaf cb"<<cb<<std::endl;
}
template <typename Op, typename T1>
inline void CBFromExpression(int &cb,const LatticeUnaryExpression<Op,T1 > &expr)
{
  CBFromExpression(cb,std::get<0>(expr.second));// recurse 
  //  std::cout<<GridLogMessage<<"Unary node cb "<<cb<<std::endl;
}

template <typename Op, typename T1, typename T2>
inline void CBFromExpression(int &cb,const LatticeBinaryExpression<Op,T1,T2> &expr) 
{
  CBFromExpression(cb,std::get<0>(expr.second));// recurse
  CBFromExpression(cb,std::get<1>(expr.second));
  //  std::cout<<GridLogMessage<<"Binary node cb "<<cb<<std::endl;
}
template <typename Op, typename T1, typename T2, typename T3>
inline void CBFromExpression( int &cb,const LatticeTrinaryExpression<Op,T1,T2,T3 > &expr) 
{
  CBFromExpression(cb,std::get<0>(expr.second));// recurse
  CBFromExpression(cb,std::get<1>(expr.second));
  CBFromExpression(cb,std::get<2>(expr.second));
  //  std::cout<<GridLogMessage<<"Trinary node cb "<<cb<<std::endl;
}

////////////////////////////////////////////
// Unary operators and funcs
////////////////////////////////////////////
#define GridUnopClass(name,ret)\
template <class arg> struct name\
{\
  static auto inline func(const arg a)-> decltype(ret) { return ret; } \
};

GridUnopClass(UnarySub,-a);
GridUnopClass(UnaryNot,Not(a));
GridUnopClass(UnaryAdj,adj(a));
GridUnopClass(UnaryConj,conjugate(a));
GridUnopClass(UnaryTrace,trace(a));
GridUnopClass(UnaryTranspose,transpose(a));
GridUnopClass(UnaryTa,Ta(a));
GridUnopClass(UnaryProjectOnGroup,ProjectOnGroup(a));
GridUnopClass(UnaryReal,real(a));
GridUnopClass(UnaryImag,imag(a));
GridUnopClass(UnaryToReal,toReal(a));
GridUnopClass(UnaryToComplex,toComplex(a));
GridUnopClass(UnaryAbs,abs(a));
GridUnopClass(UnarySqrt,sqrt(a));
GridUnopClass(UnaryRsqrt,rsqrt(a));
GridUnopClass(UnarySin,sin(a));
GridUnopClass(UnaryCos,cos(a));
GridUnopClass(UnaryLog,log(a));
GridUnopClass(UnaryExp,exp(a));

////////////////////////////////////////////
// Binary operators
////////////////////////////////////////////
#define GridBinOpClass(name,combination)\
template <class left,class right>\
struct name\
{\
  static auto inline func(const left &lhs,const right &rhs)-> decltype(combination) const \
    {\
      return combination;\
    }\
}
GridBinOpClass(BinaryAdd,lhs+rhs);
GridBinOpClass(BinarySub,lhs-rhs);
GridBinOpClass(BinaryMul,lhs*rhs);

GridBinOpClass(BinaryAnd   ,lhs&rhs);
GridBinOpClass(BinaryOr    ,lhs|rhs);
GridBinOpClass(BinaryAndAnd,lhs&&rhs);
GridBinOpClass(BinaryOrOr  ,lhs||rhs);

////////////////////////////////////////////////////
// Trinary conditional op
////////////////////////////////////////////////////
#define GridTrinOpClass(name,combination)\
template <class predicate,class left, class right>	\
struct name\
{\
  static auto inline func(const predicate &pred,const left &lhs,const right &rhs)-> decltype(combination) const \
    {\
      return combination;\
    }\
}

GridTrinOpClass(TrinaryWhere,(predicatedWhere<predicate, \
			       typename std::remove_reference<left>::type, \
			       typename std::remove_reference<right>::type> (pred,lhs,rhs)));

////////////////////////////////////////////
// Operator syntactical glue
////////////////////////////////////////////
 
#define GRID_UNOP(name)   name<decltype(eval(0, arg))>
#define GRID_BINOP(name)  name<decltype(eval(0, lhs)), decltype(eval(0, rhs))>
#define GRID_TRINOP(name) name<decltype(eval(0, pred)), decltype(eval(0, lhs)), decltype(eval(0, rhs))>

#define GRID_DEF_UNOP(op, name)\
template <typename T1,\
  typename std::enable_if<is_lattice<T1>::value||is_lattice_expr<T1>::value, T1>::type* = nullptr> inline auto op(const T1 &arg) \
  -> decltype(LatticeUnaryExpression<GRID_UNOP(name),const T1&>(std::make_pair(GRID_UNOP(name)(),std::forward_as_tuple(arg)))) \
{ return LatticeUnaryExpression<GRID_UNOP(name), const T1 &>(std::make_pair(GRID_UNOP(name)(),std::forward_as_tuple(arg))); }

#define GRID_BINOP_LEFT(op, name)\
template <typename T1,typename T2,\
          typename std::enable_if<is_lattice<T1>::value||is_lattice_expr<T1>::value, T1>::type* = nullptr>\
inline auto op(const T1 &lhs,const T2&rhs) \
  -> decltype(LatticeBinaryExpression<GRID_BINOP(name),const T1&,const T2 &>(std::make_pair(GRID_BINOP(name)(),\
											    std::forward_as_tuple(lhs, rhs)))) \
{\
 return LatticeBinaryExpression<GRID_BINOP(name), const T1 &, const T2 &>(std::make_pair(GRID_BINOP(name)(),\
									  std::forward_as_tuple(lhs, rhs))); \
}

#define GRID_BINOP_RIGHT(op, name)\
 template <typename T1,typename T2,\
           typename std::enable_if<!is_lattice<T1>::value && !is_lattice_expr<T1>::value, T1>::type* = nullptr,\
           typename std::enable_if< is_lattice<T2>::value ||  is_lattice_expr<T2>::value, T2>::type* = nullptr> \
inline auto op(const T1 &lhs,const T2&rhs)			\
  -> decltype(LatticeBinaryExpression<GRID_BINOP(name),const T1&,const T2 &>(std::make_pair(GRID_BINOP(name)(),\
											    std::forward_as_tuple(lhs, rhs)))) \
{\
 return LatticeBinaryExpression<GRID_BINOP(name), const T1 &, const T2 &>(std::make_pair(GRID_BINOP(name)(),\
								          std::forward_as_tuple(lhs, rhs))); \
}

#define GRID_DEF_BINOP(op, name)\
 GRID_BINOP_LEFT(op,name);\
 GRID_BINOP_RIGHT(op,name);


#define GRID_DEF_TRINOP(op, name)\
template <typename T1,typename T2,typename T3> inline auto op(const T1 &pred,const T2&lhs,const T3 &rhs) \
  -> decltype(LatticeTrinaryExpression<GRID_TRINOP(name),const T1&,const T2 &,const T3&>(std::make_pair(GRID_TRINOP(name)(),\
										   std::forward_as_tuple(pred,lhs,rhs)))) \
{\
  return LatticeTrinaryExpression<GRID_TRINOP(name), const T1 &, const T2 &,const T3&>(std::make_pair(GRID_TRINOP(name)(), \
										 std::forward_as_tuple(pred,lhs, rhs))); \
}
////////////////////////
//Operator definitions
////////////////////////

GRID_DEF_UNOP(operator -,UnarySub);
GRID_DEF_UNOP(Not,UnaryNot);
GRID_DEF_UNOP(operator !,UnaryNot);
GRID_DEF_UNOP(adj,UnaryAdj);
GRID_DEF_UNOP(conjugate,UnaryConj);
GRID_DEF_UNOP(trace,UnaryTrace);
GRID_DEF_UNOP(transpose,UnaryTranspose);
GRID_DEF_UNOP(Ta,UnaryTa);
GRID_DEF_UNOP(ProjectOnGroup,UnaryProjectOnGroup);
GRID_DEF_UNOP(real,UnaryReal);
GRID_DEF_UNOP(imag,UnaryImag);
GRID_DEF_UNOP(toReal,UnaryToReal);
GRID_DEF_UNOP(toComplex,UnaryToComplex);
GRID_DEF_UNOP(abs  ,UnaryAbs); //abs overloaded in cmath C++98; DON'T do the abs-fabs-dabs-labs thing
GRID_DEF_UNOP(sqrt ,UnarySqrt);
GRID_DEF_UNOP(rsqrt,UnaryRsqrt);
GRID_DEF_UNOP(sin  ,UnarySin);
GRID_DEF_UNOP(cos  ,UnaryCos);
GRID_DEF_UNOP(log  ,UnaryLog);
GRID_DEF_UNOP(exp  ,UnaryExp);

GRID_DEF_BINOP(operator+,BinaryAdd);
GRID_DEF_BINOP(operator-,BinarySub);
GRID_DEF_BINOP(operator*,BinaryMul);

GRID_DEF_BINOP(operator&,BinaryAnd);
GRID_DEF_BINOP(operator|,BinaryOr);
GRID_DEF_BINOP(operator&&,BinaryAndAnd);
GRID_DEF_BINOP(operator||,BinaryOrOr);

GRID_DEF_TRINOP(where,TrinaryWhere);

/////////////////////////////////////////////////////////////
// Closure convenience to force expression to evaluate
/////////////////////////////////////////////////////////////
template<class Op,class T1>
  auto closure(const LatticeUnaryExpression<Op,T1> & expr)
  -> Lattice<decltype(expr.first.func(eval(0,std::get<0>(expr.second))))>
{
  Lattice<decltype(expr.first.func(eval(0,std::get<0>(expr.second))))> ret(expr);
  return ret;
}
template<class Op,class T1, class T2>
  auto closure(const LatticeBinaryExpression<Op,T1,T2> & expr)
  -> Lattice<decltype(expr.first.func(eval(0,std::get<0>(expr.second)),
				      eval(0,std::get<1>(expr.second))))>
{
  Lattice<decltype(expr.first.func(eval(0,std::get<0>(expr.second)),
				   eval(0,std::get<1>(expr.second))))> ret(expr);
  return ret;
}
template<class Op,class T1, class T2, class T3>
  auto closure(const LatticeTrinaryExpression<Op,T1,T2,T3> & expr)
  -> Lattice<decltype(expr.first.func(eval(0,std::get<0>(expr.second)),
				      eval(0,std::get<1>(expr.second)),
				      eval(0,std::get<2>(expr.second))))>
{
  Lattice<decltype(expr.first.func(eval(0,std::get<0>(expr.second)),
				   eval(0,std::get<1>(expr.second)),
				   eval(0,std::get<2>(expr.second))))> ret(expr);
  return ret;
}

#undef GRID_UNOP
#undef GRID_BINOP
#undef GRID_TRINOP

#undef GRID_DEF_UNOP
#undef GRID_DEF_BINOP
#undef GRID_DEF_TRINOP

}

#if 0
using namespace Grid;
 	      
 int main(int argc,char **argv){
   
   Lattice<double> v1(16);
   Lattice<double> v2(16);
   Lattice<double> v3(16);

   BinaryAdd<double,double> tmp;
   LatticeBinaryExpression<BinaryAdd<double,double>,Lattice<double> &,Lattice<double> &> 
     expr(std::make_pair(tmp,
	  std::forward_as_tuple(v1,v2)));
   tmp.func(eval(0,v1),eval(0,v2));

   auto var = v1+v2;
   std::cout<<GridLogMessage<<typeid(var).name()<<std::endl;

   v3=v1+v2;
   v3=v1+v2+v1*v2;
 };

void testit(Lattice<double> &v1,Lattice<double> &v2,Lattice<double> &v3)
{
   v3=v1+v2+v1*v2;
}
#endif

#endif
