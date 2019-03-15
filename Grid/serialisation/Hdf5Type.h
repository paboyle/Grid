#ifndef GRID_SERIALISATION_HDF5_TYPE_H
#define GRID_SERIALISATION_HDF5_TYPE_H

#include <H5Cpp.h>
#include <complex>
#include <memory>

#ifndef H5_NO_NAMESPACE
#define H5NS H5
#endif

#define HDF5_NATIVE_TYPE(predType, cType)\
template <>\
class Hdf5Type<cType>\
{\
public:\
  static inline const H5NS::DataType & type(void)\
  {\
    return H5NS::PredType::predType;\
  }\
  static constexpr bool isNative = true;\
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
  template <typename T> class Hdf5Type
  {
  public:
    static constexpr bool isNative = false;
  };
  
  DEFINE_HDF5_NATIVE_TYPES;
  
  template <typename R>
  class Hdf5Type<std::complex<R>>
  {
  public:
    static inline const H5NS::DataType & type(void)
    {
      if (typePtr_ == nullptr)
      {
        typePtr_.reset(new H5NS::CompType(sizeof(std::complex<R>)));
        typePtr_->insertMember("re", 0,         Hdf5Type<R>::type());
        typePtr_->insertMember("im", sizeof(R), Hdf5Type<R>::type());
      }

      return *typePtr_;
    }
    static constexpr bool isNative = false;
  private:
    static std::unique_ptr<H5NS::CompType> typePtr_;
  };
  
  template <typename R>
  std::unique_ptr<H5NS::CompType> Hdf5Type<std::complex<R>>::typePtr_ = nullptr;
}

#undef HDF5_NATIVE_TYPE

#endif /* GRID_SERIALISATION_HDF5_TYPE_H */
