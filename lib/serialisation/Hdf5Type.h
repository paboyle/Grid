#ifndef GRID_SERIALISATION_HDF5_TYPE_H
#define GRID_SERIALISATION_HDF5_TYPE_H

#include <H5Cpp.h>
#include <vector>

#ifndef H5_NO_NAMESPACE
#define H5NS H5
#endif

#define HDF5_NATIVE_TYPE(predType, cType)\
template <>\
struct Hdf5Type<cType>\
{\
static inline const H5NS::PredType & type(void)\
{\
  return H5NS::PredType::predType;\
}\
static constexpr bool       isNative = true;\
};

#define DEFINE_HDF5_NATIVE_TYPES \
HDF5_NATIVE_TYPE(NATIVE_B8,      bool);\
HDF5_NATIVE_TYPE(NATIVE_CHAR,    char);\
HDF5_NATIVE_TYPE(NATIVE_SCHAR,   signed char);\
HDF5_NATIVE_TYPE(NATIVE_UCHAR,   unsigned char);\
HDF5_NATIVE_TYPE(NATIVE_SHORT,   short);\
HDF5_NATIVE_TYPE(NATIVE_USHORT,  unsigned short);\
HDF5_NATIVE_TYPE(NATIVE_INT,     int);\
HDF5_NATIVE_TYPE(NATIVE_UINT,    unsigned int);\
HDF5_NATIVE_TYPE(NATIVE_LONG,    long);\
HDF5_NATIVE_TYPE(NATIVE_ULONG,   unsigned long);\
HDF5_NATIVE_TYPE(NATIVE_LLONG,   long long);\
HDF5_NATIVE_TYPE(NATIVE_ULLONG,  unsigned long long);\
HDF5_NATIVE_TYPE(NATIVE_FLOAT,   float);\
HDF5_NATIVE_TYPE(NATIVE_DOUBLE,  double);\
HDF5_NATIVE_TYPE(NATIVE_LDOUBLE, long double);

namespace Grid
{
  template <typename T> struct Hdf5Type
  {
    static constexpr bool isNative = false;
  };
  
  DEFINE_HDF5_NATIVE_TYPES;
}

#undef HDF5_NATIVE_TYPE

#endif /* GRID_SERIALISATION_HDF5_TYPE_H */
