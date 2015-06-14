#ifndef GRID_TENSOR_LOGICAL_H
#define GRID_TENSOR_LOGICAL_H

namespace Grid {

#define LOGICAL_BINOP(Op)\
template<class v> strong_inline iScalar<v> operator Op (const iScalar<v>& lhs,const iScalar<v>& rhs) \
{\
  iScalar<v> ret;\
  ret._internal = lhs._internal Op rhs._internal ;\
  return ret;\
}\
template<class l> strong_inline iScalar<l> operator Op (const iScalar<l>& lhs,Integer rhs) \
{\
  typename iScalar<l>::scalar_type t; t=rhs;\
  typename iScalar<l>::tensor_reduced srhs; srhs=t;\
  return lhs Op srhs;\
}\
template<class l> strong_inline iScalar<l> operator Op (Integer lhs,const iScalar<l>& rhs) \
{\
  typename iScalar<l>::scalar_type t;t=lhs;\
  typename iScalar<l>::tensor_reduced slhs;slhs=t;\
  return slhs Op rhs;\
}

LOGICAL_BINOP(|);
LOGICAL_BINOP(&);
LOGICAL_BINOP(||);
LOGICAL_BINOP(&&);

}
#endif
