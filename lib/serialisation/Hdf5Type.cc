#include "Hdf5Type.h"

using namespace Grid;

#define HDF5_NATIVE_TYPE(predType, cType)\
const H5NS::PredType * Hdf5Type<cType>::type = &H5NS::PredType::predType;

DEFINE_HDF5_NATIVE_TYPES;
