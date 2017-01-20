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
static const H5NS::PredType *type;\
static constexpr bool       isNative = true;\
};

#define DEFINE_HDF5_NATIVE_TYPES \
HDF5_NATIVE_TYPE(STD_B8LE,   bool);\
HDF5_NATIVE_TYPE(STD_I8LE,   char);\
HDF5_NATIVE_TYPE(STD_U8LE,   unsigned char);\
HDF5_NATIVE_TYPE(STD_I16LE,  short);\
HDF5_NATIVE_TYPE(STD_U16LE,  unsigned short);\
HDF5_NATIVE_TYPE(STD_I32LE,  int);\
HDF5_NATIVE_TYPE(STD_U32LE,  unsigned int);\
HDF5_NATIVE_TYPE(STD_I64LE,  long);\
HDF5_NATIVE_TYPE(STD_U64LE,  unsigned long);\
HDF5_NATIVE_TYPE(STD_I64LE,  long long);\
HDF5_NATIVE_TYPE(STD_U64LE,  unsigned long long);\
HDF5_NATIVE_TYPE(IEEE_F32LE, float);\
HDF5_NATIVE_TYPE(IEEE_F64LE, double);\
HDF5_NATIVE_TYPE(IEEE_F64LE, long double);

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
