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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
			   /*  END LEGAL */
#ifndef GRID_LATTICE_ET_H
#define GRID_LATTICE_ET_H

#include <iostream>
#include <tuple>
#include <typeinfo>
#include <vector>

NAMESPACE_BEGIN(Grid);

////////////////////////////////////////////////////
// Predicated where support
////////////////////////////////////////////////////
template <class iobj, class vobj, class robj>
accelerator_inline vobj predicatedWhere(const iobj &predicate, const vobj &iftrue,
                            const robj &iffalse) {
  typename std::remove_const<vobj>::type ret;

  typedef typename vobj::scalar_object scalar_object;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  const int Nsimd = vobj::vector_type::Nsimd();

  ExtractBuffer<Integer> mask(Nsimd);
  ExtractBuffer<scalar_object> truevals(Nsimd);
  ExtractBuffer<scalar_object> falsevals(Nsimd);

  extract(iftrue, truevals);
  extract(iffalse, falsevals);
  extract<vInteger, Integer>(TensorRemove(predicate), mask);

  for (int s = 0; s < Nsimd; s++) {
    if (mask[s]) falsevals[s] = truevals[s];
  }

  merge(ret, falsevals);
  return ret;
}

/////////////////////////////////////////////////////
//Specialization of getVectorType for lattices
/////////////////////////////////////////////////////
template<typename T>
struct getVectorType<Lattice<T> >{
  typedef typename Lattice<T>::vector_object type;
};

////////////////////////////////////////////
//--  recursive evaluation of expressions; --
// handle leaves of syntax tree
///////////////////////////////////////////////////
template<class sobj> accelerator_inline 
sobj eval(const unsigned int ss, const sobj &arg)
{
  return arg;
}

template <class lobj> accelerator_inline 
const lobj & eval(const unsigned int ss, const LatticeView<lobj> &arg) 
{
  return arg[ss];
}
template <class lobj> accelerator_inline 
const lobj & eval(const unsigned int ss, const Lattice<lobj> &arg) 
{
  auto view = arg.View();
  return view[ss];
}

///////////////////////////////////////////////////
// handle nodes in syntax tree- eval one operand
///////////////////////////////////////////////////
template <typename Op, typename T1> accelerator_inline 
auto eval(const unsigned int ss, const LatticeUnaryExpression<Op, T1> &expr)  
  -> decltype(expr.op.func( eval(ss, expr.arg1)))
{
  return expr.op.func( eval(ss, expr.arg1) );
}
///////////////////////
// eval two operands
///////////////////////
template <typename Op, typename T1, typename T2> accelerator_inline
auto eval(const unsigned int ss, const LatticeBinaryExpression<Op, T1, T2> &expr)  
  -> decltype(expr.op.func( eval(ss,expr.arg1),eval(ss,expr.arg2)))
{
  return expr.op.func( eval(ss,expr.arg1), eval(ss,expr.arg2) );
}
///////////////////////
// eval three operands
///////////////////////
template <typename Op, typename T1, typename T2, typename T3> accelerator_inline
auto eval(const unsigned int ss, const LatticeTrinaryExpression<Op, T1, T2, T3> &expr)  
  -> decltype(expr.op.func(eval(ss, expr.arg1), eval(ss, expr.arg2), eval(ss, expr.arg3)))
{
  return expr.op.func(eval(ss, expr.arg1), eval(ss, expr.arg2), eval(ss, expr.arg3));
}

//////////////////////////////////////////////////////////////////////////
// Obtain the grid from an expression, ensuring conformable. This must follow a
// tree recursion; must retain grid pointer in the LatticeView class which sucks
// Use a different method, and make it void *.
// Perhaps a conformable method.
//////////////////////////////////////////////////////////////////////////
template <class T1,typename std::enable_if<is_lattice<T1>::value, T1>::type * = nullptr>
accelerator_inline void GridFromExpression(GridBase *&grid, const T1 &lat)  // Lattice leaf
{
  lat.Conformable(grid);
}

template <class T1,typename std::enable_if<!is_lattice<T1>::value, T1>::type * = nullptr>
accelerator_inline 
void GridFromExpression(GridBase *&grid,const T1 &notlat)  // non-lattice leaf
{}

template <typename Op, typename T1>
accelerator_inline 
void GridFromExpression(GridBase *&grid,const LatticeUnaryExpression<Op, T1> &expr) 
{
  GridFromExpression(grid, expr.arg1);  // recurse
}

template <typename Op, typename T1, typename T2>
accelerator_inline 
void GridFromExpression(GridBase *&grid, const LatticeBinaryExpression<Op, T1, T2> &expr) 
{
  GridFromExpression(grid, expr.arg1);  // recurse
  GridFromExpression(grid, expr.arg2);
}
template <typename Op, typename T1, typename T2, typename T3>
accelerator_inline 
void GridFromExpression(GridBase *&grid, const LatticeTrinaryExpression<Op, T1, T2, T3> &expr) 
{
  GridFromExpression(grid, expr.arg1);  // recurse
  GridFromExpression(grid, expr.arg2);  // recurse
  GridFromExpression(grid, expr.arg3);  // recurse
}

//////////////////////////////////////////////////////////////////////////
// Obtain the CB from an expression, ensuring conformable. This must follow a
// tree recursion
//////////////////////////////////////////////////////////////////////////
template <class T1,typename std::enable_if<is_lattice<T1>::value, T1>::type * = nullptr>
inline void CBFromExpression(int &cb, const T1 &lat)  // Lattice leaf
{
  if ((cb == Odd) || (cb == Even)) {
    assert(cb == lat.Checkerboard());
  }
  cb = lat.Checkerboard();
}
template <class T1,typename std::enable_if<!is_lattice<T1>::value, T1>::type * = nullptr>
inline void CBFromExpression(int &cb, const T1 &notlat)  // non-lattice leaf
{
}

template <typename Op, typename T1> inline 
void CBFromExpression(int &cb,const LatticeUnaryExpression<Op, T1> &expr) 
{
  CBFromExpression(cb, expr.arg1);  // recurse AST
}

template <typename Op, typename T1, typename T2> inline 
void CBFromExpression(int &cb,const LatticeBinaryExpression<Op, T1, T2> &expr) 
{
  CBFromExpression(cb, expr.arg1);  // recurse AST
  CBFromExpression(cb, expr.arg2);  // recurse AST
}
template <typename Op, typename T1, typename T2, typename T3>
inline void CBFromExpression(int &cb, const LatticeTrinaryExpression<Op, T1, T2, T3> &expr) 
{
  CBFromExpression(cb, expr.arg1);  // recurse AST
  CBFromExpression(cb, expr.arg2);  // recurse AST
  CBFromExpression(cb, expr.arg3);  // recurse AST
}

////////////////////////////////////////////
// Unary operators and funcs
////////////////////////////////////////////
#define GridUnopClass(name, ret)					\
  template <class arg>							\
  struct name {								\
    static auto accelerator_inline func(const arg a) -> decltype(ret) { return ret; } \
  };

GridUnopClass(UnarySub, -a);
GridUnopClass(UnaryNot, Not(a));
GridUnopClass(UnaryAdj, adj(a));
GridUnopClass(UnaryConj, conjugate(a));
GridUnopClass(UnaryTrace, trace(a));
GridUnopClass(UnaryTranspose, transpose(a));
GridUnopClass(UnaryTa, Ta(a));
GridUnopClass(UnaryProjectOnGroup, ProjectOnGroup(a));
GridUnopClass(UnaryReal, real(a));
GridUnopClass(UnaryImag, imag(a));
GridUnopClass(UnaryToReal, toReal(a));
GridUnopClass(UnaryToComplex, toComplex(a));
GridUnopClass(UnaryTimesI, timesI(a));
GridUnopClass(UnaryTimesMinusI, timesMinusI(a));
GridUnopClass(UnaryAbs, abs(a));
GridUnopClass(UnarySqrt, sqrt(a));
GridUnopClass(UnaryRsqrt, rsqrt(a));
GridUnopClass(UnarySin, sin(a));
GridUnopClass(UnaryCos, cos(a));
GridUnopClass(UnaryAsin, asin(a));
GridUnopClass(UnaryAcos, acos(a));
GridUnopClass(UnaryLog, log(a));
GridUnopClass(UnaryExp, exp(a));

////////////////////////////////////////////
// Binary operators
////////////////////////////////////////////
#define GridBinOpClass(name, combination)			\
  template <class left, class right>				\
  struct name {							\
    static auto accelerator_inline				\
    func(const left &lhs, const right &rhs)			\
      -> decltype(combination) const				\
    {								\
      return combination;					\
    }								\
  };

GridBinOpClass(BinaryAdd, lhs + rhs);
GridBinOpClass(BinarySub, lhs - rhs);
GridBinOpClass(BinaryMul, lhs *rhs);
GridBinOpClass(BinaryDiv, lhs /rhs);
GridBinOpClass(BinaryAnd, lhs &rhs);
GridBinOpClass(BinaryOr, lhs | rhs);
GridBinOpClass(BinaryAndAnd, lhs &&rhs);
GridBinOpClass(BinaryOrOr, lhs || rhs);

////////////////////////////////////////////////////
// Trinary conditional op
////////////////////////////////////////////////////
#define GridTrinOpClass(name, combination)				\
  template <class predicate, class left, class right>			\
  struct name {								\
    static auto accelerator_inline					\
    func(const predicate &pred, const left &lhs, const right &rhs)	\
      -> decltype(combination) const					\
    {									\
      return combination;						\
    }									\
  };

GridTrinOpClass(TrinaryWhere,
		(predicatedWhere<predicate, 
		 typename std::remove_reference<left>::type,
		 typename std::remove_reference<right>::type>(pred, lhs,rhs)));

////////////////////////////////////////////
// Operator syntactical glue
////////////////////////////////////////////

#define GRID_UNOP(name)   name<decltype(eval(0, arg))>
#define GRID_BINOP(name)  name<decltype(eval(0, lhs)), decltype(eval(0, rhs))>
#define GRID_TRINOP(name) name<decltype(eval(0, pred)), decltype(eval(0, lhs)), decltype(eval(0, rhs))>

#define GRID_DEF_UNOP(op, name)						\
  template <typename T1, typename std::enable_if<is_lattice<T1>::value||is_lattice_expr<T1>::value,T1>::type * = nullptr> \
  inline auto op(const T1 &arg) ->decltype(LatticeUnaryExpression<GRID_UNOP(name),T1>(GRID_UNOP(name)(), arg)) \
  {									\
    return     LatticeUnaryExpression<GRID_UNOP(name),T1>(GRID_UNOP(name)(), arg); \
  }

#define GRID_BINOP_LEFT(op, name)					\
  template <typename T1, typename T2,					\
            typename std::enable_if<is_lattice<T1>::value||is_lattice_expr<T1>::value,T1>::type * = nullptr> \
  inline auto op(const T1 &lhs, const T2 &rhs)				\
    ->decltype(LatticeBinaryExpression<GRID_BINOP(name),T1,T2>(GRID_BINOP(name)(),lhs,rhs)) \
  {									\
    return     LatticeBinaryExpression<GRID_BINOP(name),T1,T2>(GRID_BINOP(name)(),lhs,rhs);\
  }

#define GRID_BINOP_RIGHT(op, name)					\
  template <typename T1, typename T2,					\
            typename std::enable_if<!is_lattice<T1>::value&&!is_lattice_expr<T1>::value,T1>::type * = nullptr, \
            typename std::enable_if< is_lattice<T2>::value|| is_lattice_expr<T2>::value,T2>::type * = nullptr> \
  inline auto op(const T1 &lhs, const T2 &rhs)				\
    ->decltype(LatticeBinaryExpression<GRID_BINOP(name),T1,T2>(GRID_BINOP(name)(),lhs, rhs)) \
  {									\
    return     LatticeBinaryExpression<GRID_BINOP(name),T1,T2>(GRID_BINOP(name)(),lhs, rhs); \
  }

#define GRID_DEF_BINOP(op, name)		\
  GRID_BINOP_LEFT(op, name);			\
  GRID_BINOP_RIGHT(op, name);

#define GRID_DEF_TRINOP(op, name)					\
  template <typename T1, typename T2, typename T3>			\
  inline auto op(const T1 &pred, const T2 &lhs, const T3 &rhs)		\
    ->decltype(LatticeTrinaryExpression<GRID_TRINOP(name),T1,T2,T3>(GRID_TRINOP(name)(),pred, lhs, rhs)) \
  {									\
    return LatticeTrinaryExpression<GRID_TRINOP(name),T1,T2,T3>(GRID_TRINOP(name)(),pred, lhs, rhs); \
  }

////////////////////////
// Operator definitions
////////////////////////
GRID_DEF_UNOP(operator-, UnarySub);
GRID_DEF_UNOP(Not, UnaryNot);
GRID_DEF_UNOP(operator!, UnaryNot);
GRID_DEF_UNOP(adj, UnaryAdj);
GRID_DEF_UNOP(conjugate, UnaryConj);
GRID_DEF_UNOP(trace, UnaryTrace);
GRID_DEF_UNOP(transpose, UnaryTranspose);
GRID_DEF_UNOP(Ta, UnaryTa);
GRID_DEF_UNOP(ProjectOnGroup, UnaryProjectOnGroup);
GRID_DEF_UNOP(real, UnaryReal);
GRID_DEF_UNOP(imag, UnaryImag);
GRID_DEF_UNOP(toReal, UnaryToReal);
GRID_DEF_UNOP(toComplex, UnaryToComplex);
GRID_DEF_UNOP(timesI, UnaryTimesI);
GRID_DEF_UNOP(timesMinusI, UnaryTimesMinusI);
GRID_DEF_UNOP(abs, UnaryAbs);  // abs overloaded in cmath C++98; DON'T do the
                               // abs-fabs-dabs-labs thing
GRID_DEF_UNOP(sqrt, UnarySqrt);
GRID_DEF_UNOP(rsqrt, UnaryRsqrt);
GRID_DEF_UNOP(sin, UnarySin);
GRID_DEF_UNOP(cos, UnaryCos);
GRID_DEF_UNOP(asin, UnaryAsin);
GRID_DEF_UNOP(acos, UnaryAcos);
GRID_DEF_UNOP(log, UnaryLog);
GRID_DEF_UNOP(exp, UnaryExp);

GRID_DEF_BINOP(operator+, BinaryAdd);
GRID_DEF_BINOP(operator-, BinarySub);
GRID_DEF_BINOP(operator*, BinaryMul);
GRID_DEF_BINOP(operator/, BinaryDiv);

GRID_DEF_BINOP(operator&, BinaryAnd);
GRID_DEF_BINOP(operator|, BinaryOr);
GRID_DEF_BINOP(operator&&, BinaryAndAnd);
GRID_DEF_BINOP(operator||, BinaryOrOr);

GRID_DEF_TRINOP(where, TrinaryWhere);

/////////////////////////////////////////////////////////////
// Closure convenience to force expression to evaluate
/////////////////////////////////////////////////////////////
template <class Op, class T1>
auto closure(const LatticeUnaryExpression<Op, T1> &expr)
  -> Lattice<decltype(expr.op.func(eval(0, expr.arg1)))> 
{
  Lattice<decltype(expr.op.func(eval(0, expr.arg1)))> ret(expr);
  return ret;
}
template <class Op, class T1, class T2>
auto closure(const LatticeBinaryExpression<Op, T1, T2> &expr)
  -> Lattice<decltype(expr.op.func(eval(0, expr.arg1),eval(0, expr.arg2)))> 
{
  Lattice<decltype(expr.op.func(eval(0, expr.arg1),eval(0, expr.arg2)))> ret(expr);
  return ret;
}
template <class Op, class T1, class T2, class T3>
auto closure(const LatticeTrinaryExpression<Op, T1, T2, T3> &expr)
  -> Lattice<decltype(expr.op.func(eval(0, expr.arg1),
				   eval(0, expr.arg2),
				   eval(0, expr.arg3)))> 
{
  Lattice<decltype(expr.op.func(eval(0, expr.arg1),
				eval(0, expr.arg2),
				eval(0, expr.arg3)))>  ret(expr);
  return ret;
}

#undef GRID_UNOP
#undef GRID_BINOP
#undef GRID_TRINOP

#undef GRID_DEF_UNOP
#undef GRID_DEF_BINOP
#undef GRID_DEF_TRINOP

NAMESPACE_END(Grid);

#endif
