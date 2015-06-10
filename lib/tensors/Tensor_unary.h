#ifndef GRID_TENSOR_UNARY_H
#define GRID_TENSOR_UNARY_H
namespace Grid {

#define UNARY_REAL(func)\
template<class obj> inline auto func(const iScalar<obj> &z) -> iScalar<obj>\
{\
    iScalar<obj> ret;\
    ret._internal = func( (z._internal));\
    return ret;\
}\
template<class obj,int N> inline auto func(const iVector<obj,N> &z) -> iVector<obj,N>\
{\
    iVector<obj,N> ret;\
    for(int c1=0;c1<N;c1++){\
      ret._internal[c1] = func( (z._internal[c1]));\
    }\
    return ret;\
}\
template<class obj,int N> inline auto func(const iMatrix<obj,N> &z) -> iMatrix<obj,N>\
{\
    iMatrix<obj,N> ret;\
    for(int c1=0;c1<N;c1++){\
    for(int c2=0;c2<N;c2++){\
      ret._internal[c1][c2] = func( (z._internal[c1][c2]));\
    }}\
    return ret;\
}


#define BINARY_RSCALAR(func,scal)					\
template<class obj> inline iScalar<obj> func(const iScalar<obj> &z,scal y)	\
{\
    iScalar<obj> ret;\
    ret._internal = func(z._internal,y);	\
    return ret;\
}\
 template<class obj,int N> inline iVector<obj,N> func(const iVector<obj,N> &z,scal y) \
{\
    iVector<obj,N> ret;\
    for(int c1=0;c1<N;c1++){\
      ret._internal[c1] = func(z._internal[c1],y);	\
    }\
    return ret;\
}\
 template<class obj,int N> inline  iMatrix<obj,N> func(const iMatrix<obj,N> &z, scal y) \
{\
    iMatrix<obj,N> ret;\
    for(int c1=0;c1<N;c1++){\
    for(int c2=0;c2<N;c2++){\
      ret._internal[c1][c2] = func(z._internal[c1][c2],y);	\
    }}\
    return ret;\
}

UNARY_REAL(sqrt);
UNARY_REAL(rsqrt);
UNARY_REAL(sin);
UNARY_REAL(cos);

BINARY_RSCALAR(mod,Integer);
BINARY_RSCALAR(pow,RealD);


}
#endif
