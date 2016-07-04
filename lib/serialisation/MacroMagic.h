#ifndef GRID_MACRO_MAGIC_H
#define GRID_MACRO_MAGIC_H

///////////////////////////////////////////
// Strong credit to :
//
// http://jhnet.co.uk/articles/cpp_magic
//
// 
// "The C Pre-Processor (CPP) is the somewhat basic macro system used by the C 
// programming language to implement features such as #include and #define 
// which allow very simple text-substitutions to be carried out at compile time.
// In this article we abuse the humble #define to implement if-statements and iteration.
//
// Before we begin, a disclaimer: these tricks, while perfectly valid C, should not be 
// considered good development practice and should almost certainly not be used for "real work".
// That said it can totally be used for fun home-automation projects...
//
// https://github.com/18sg/uSHET/blob/master/lib/cpp_magic.h
//
// The cpp_magic.h license (prior to my modifications):
/*
Copyright (c) 2014 Thomas Nixon, Jonathan Heathcote

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
///////////////////////////////////////////

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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// These are the Grid extensions to serialisation (beyond popping first AND second)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define GRID_MACRO_MEMBER(A,B)        A B;
#define GRID_MACRO_OS_WRITE_MEMBER(A,B) os<< #A <<" "#B <<" = "<< obj. B <<" ; " <<std::endl;
#define GRID_MACRO_READ_MEMBER(A,B) Grid::read(RD,#B,obj. B);
#define GRID_MACRO_WRITE_MEMBER(A,B) Grid::write(WR,#B,obj. B);

#define GRID_SERIALIZABLE_CLASS_MEMBERS(cname,...)		\
  \
  \
  GRID_MACRO_EVAL(GRID_MACRO_MAP(GRID_MACRO_MEMBER,__VA_ARGS__))		\
  \
  \
  template <typename T>\
  static inline void write(Writer<T> &WR,const std::string &s, const cname &obj){ \
    push(WR,s);\
    GRID_MACRO_EVAL(GRID_MACRO_MAP(GRID_MACRO_WRITE_MEMBER,__VA_ARGS__))	\
    pop(WR);\
  } \
  \
  \
  template <typename T>\
  static inline void read(Reader<T> &RD,const std::string &s, cname &obj){	\
    push(RD,s);\
    GRID_MACRO_EVAL(GRID_MACRO_MAP(GRID_MACRO_READ_MEMBER,__VA_ARGS__))	\
    pop(RD);\
  } \
  \
  \
  friend inline std::ostream & operator << (std::ostream &os, const cname &obj ) { \
    os<<"class "<<#cname<<" {"<<std::endl;\
    GRID_MACRO_EVAL(GRID_MACRO_MAP(GRID_MACRO_OS_WRITE_MEMBER,__VA_ARGS__))	\
      os<<"}";								\
    return os;\
  };  



#define GRID_ENUM_TYPE(obj) std::remove_reference<decltype(obj)>::type
#define GRID_MACRO_ENUMVAL(A,B) A = B,
#define GRID_MACRO_ENUMCASE(A,B) case GRID_ENUM_TYPE(obj)::A: Grid::write(WR,s,#A); break;
#define GRID_MACRO_ENUMTEST(A,B) else if (buf == #A) {obj = GRID_ENUM_TYPE(obj)::A;}
#define GRID_MACRO_ENUMCASEIO(A,B) case GRID_ENUM_TYPE(obj)::A: os << #A; break;

namespace Grid {
  template <typename U>
  class EnumIO {};
}

#define GRID_SERIALIZABLE_ENUM(name,undefname,...)\
  enum class name {\
      GRID_MACRO_EVAL(GRID_MACRO_MAP(GRID_MACRO_ENUMVAL,__VA_ARGS__))\
      undefname = -1\
  };\
  \
  template<>\
  class EnumIO<name> {\
    public:\
      template <typename T>\
      static inline void write(Writer<T> &WR,const std::string &s, const name &obj){ \
        switch (obj) {\
          GRID_MACRO_EVAL(GRID_MACRO_MAP(GRID_MACRO_ENUMCASE,__VA_ARGS__))\
          default: Grid::write(WR,s,#undefname); break;\
        }\
      }\
      \
      template <typename T>\
      static inline void read(Reader<T> &RD,const std::string &s, name &obj){ \
        std::string buf;\
        Grid::read(RD, s, buf);\
        if (buf == #undefname) {obj = name::undefname;}\
        GRID_MACRO_EVAL(GRID_MACRO_MAP(GRID_MACRO_ENUMTEST,__VA_ARGS__))\
        else {obj = name::undefname;}\
      }\
  };\
  \
  inline std::ostream & operator << (std::ostream &os, const name &obj ) { \
    switch (obj) {\
        GRID_MACRO_EVAL(GRID_MACRO_MAP(GRID_MACRO_ENUMCASEIO,__VA_ARGS__))\
        default: os << #undefname; break;\
    }\
    return os;\
  };

#endif
