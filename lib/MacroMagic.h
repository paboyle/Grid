#ifndef GRID_MACRO_MAGIC_H
#define GRID_MACRO_MAGIC_H

#define strong_inline __attribute__((always_inline)) inline

#ifndef MAX
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)>(y)?(y):(x))
#endif

#define GRID_MACRO_FIRST(a, ...) a
#define GRID_MACRO_SECOND(a, b, ...) b

#define GRID_MACRO_EMPTY()

#define GRID_MACRO_EVAL(...)     GRID_MACRO_EVAL1024(__VA_ARGS__)
#define GRID_MACRO_EVAL1024(...) GRID_MACRO_EVAL512(GRID_MACRO_EVAL512(__VA_ARGS__))
#define GRID_MACRO_EVAL512(...)  GRID_MACRO_EVAL256(GRID_MACRO_EVAL256(__VA_ARGS__))
#define GRID_MACRO_EVAL256(...)  GRID_MACRO_EVAL128(GRID_MACRO_EVAL128(__VA_ARGS__))
#define GRID_MACRO_EVAL128(...)  GRID_MACRO_EVAL64(GRID_MACRO_EVAL64(__VA_ARGS__))
#define GRID_MACRO_EVAL64(...)   GRID_MACRO_EVAL32(GRID_MACRO_EVAL32(__VA_ARGS__))
#define GRID_MACRO_EVAL32(...)   GRID_MACRO_EVAL16(GRID_MACRO_EVAL16(__VA_ARGS__))
#define GRID_MACRO_EVAL16(...)   GRID_MACRO_EVAL8(GRID_MACRO_EVAL8(__VA_ARGS__))
#define GRID_MACRO_EVAL8(...)    GRID_MACRO_EVAL4(GRID_MACRO_EVAL4(__VA_ARGS__))
#define GRID_MACRO_EVAL4(...)    GRID_MACRO_EVAL2(GRID_MACRO_EVAL2(__VA_ARGS__))
#define GRID_MACRO_EVAL2(...)    GRID_MACRO_EVAL1(GRID_MACRO_EVAL1(__VA_ARGS__))
#define GRID_MACRO_EVAL1(...) __VA_ARGS__

#define GRID_MACRO_DEFER1(m) m GRID_MACRO_EMPTY()
#define GRID_MACRO_DEFER2(m) m GRID_MACRO_EMPTY GRID_MACRO_EMPTY()()
#define GRID_MACRO_DEFER3(m) m GRID_MACRO_EMPTY GRID_MACRO_EMPTY GRID_MACRO_EMPTY()()()
#define GRID_MACRO_DEFER4(m) m GRID_MACRO_EMPTY GRID_MACRO_EMPTY GRID_MACRO_EMPTY GRID_MACRO_EMPTY()()()()

#define GRID_MACRO_IS_PROBE(...) GRID_MACRO_SECOND(__VA_ARGS__, 0)
#define GRID_MACRO_PROBE() ~, 1

#define GRID_MACRO_CAT(a,b) a ## b

#define GRID_MACRO_NOT(x) GRID_MACRO_IS_PROBE(GRID_MACRO_CAT(_GRID_MACRO_NOT_, x))
#define _GRID_MACRO_NOT_0 GRID_MACRO_PROBE()

#define GRID_MACRO_BOOL(x) GRID_MACRO_NOT(GRID_MACRO_NOT(x))

#define GRID_MACRO_IF_ELSE(condition) _GRID_MACRO_IF_ELSE(GRID_MACRO_BOOL(condition))
#define _GRID_MACRO_IF_ELSE(condition) GRID_MACRO_CAT(_GRID_MACRO_IF_, condition)

#define _GRID_MACRO_IF_1(...) __VA_ARGS__ _GRID_MACRO_IF_1_ELSE
#define _GRID_MACRO_IF_0(...)             _GRID_MACRO_IF_0_ELSE

#define _GRID_MACRO_IF_1_ELSE(...)
#define _GRID_MACRO_IF_0_ELSE(...) __VA_ARGS__

#define GRID_MACRO_HAS_ARGS(...) GRID_MACRO_BOOL(GRID_MACRO_FIRST(_GRID_MACRO_END_OF_ARGUMENTS_ __VA_ARGS__)())
#define _GRID_MACRO_END_OF_ARGUMENTS_() 0

#define GRID_MACRO_MAP(m, first, second, ...)   \
  m(first,second)                           \
  GRID_MACRO_IF_ELSE(GRID_MACRO_HAS_ARGS(__VA_ARGS__))(				       \
				 GRID_MACRO_DEFER4(_GRID_MACRO_MAP)()(m, __VA_ARGS__)   \
				     )(                                 \
				       /* Do nothing, just terminate */ \
									)
#define _GRID_MACRO_MAP() GRID_MACRO_MAP

#define GRID_MACRO_MEMBER(A,B)        A B;

#define GRID_MACRO_OS_WRITE_MEMBER(A,B) os<< #A <<" "#B <<" = "<< obj. B <<" ; " <<std::endl;

#define GRID_DECL_CLASS_MEMBERS(cname,...)		\
  GRID_MACRO_EVAL(GRID_MACRO_MAP(GRID_MACRO_MEMBER,__VA_ARGS__))		\
  friend std::ostream & operator << (std::ostream &os, const cname &obj ) {	\
    os<<"class "<<#cname<<" {"<<std::endl;\
    GRID_MACRO_EVAL(GRID_MACRO_MAP(GRID_MACRO_OS_WRITE_MEMBER,__VA_ARGS__))	\
      os<<"}";								\
    return os;\
  };  

#endif
