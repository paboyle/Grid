    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/Test_serialisation.cc

    Copyright (C) 2015-2016

Author: Guido Cossu <guido.cossu@ed.ac.uk>
Author: Antonin Portelli <antonin.portelli@me.com>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
#include <Grid/Grid.h>

using namespace Grid;
using namespace Grid::QCD;

GRID_SERIALIZABLE_ENUM(myenum, undef, red, 1, blue, 2, green, 3);

class myclass: Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(myclass,
                          myenum, e,
                          std::vector<myenum>, ve,
                          std::string, name,
                          int, x,
                          double, y,
                          bool , b,
                          std::vector<double>, array,
                          std::vector<std::vector<double> >, twodimarray,
                          std::vector<std::vector<std::vector<Complex> > >, cmplx3darray,
                          SpinColourMatrix, scm
                          );
  myclass() {}
  myclass(int i)
  : array(4,5.1)
  , twodimarray(3,std::vector<double>(5, 1.23456))
  , cmplx3darray(3,std::vector<std::vector<Complex>>(5, std::vector<Complex>(7, Complex(1.2, 3.4))))
  , ve(2, myenum::blue)
  {
    e=myenum::red;
    x=i;
    y=2*i;
    b=true;
    name="bother said pooh";
    scm()(0, 1)(2, 1) = 2.356;
    scm()(3, 0)(1, 1) = 1.323;
    scm()(2, 1)(0, 1) = 5.3336;
    scm()(0, 2)(1, 1) = 6.336;
    scm()(2, 1)(2, 2) = 7.344;
    scm()(1, 1)(2, 0) = 8.3534;
  }
};

int16_t  i16 = 1;
uint16_t u16 = 2;
int32_t  i32 = 3;
uint32_t u32 = 4;
int64_t  i64 = 5;
uint64_t u64 = 6;
float    f   = M_PI;
double   d   = 2*M_PI;
bool     b   = false;

template <typename W, typename R, typename O>
void ioTest(const std::string &filename, const O &object, const std::string &name)
{
  // writer needs to be destroyed so that writing physically happens
  {
    W writer(filename);

    write(writer, "testobject", object);
  }

  R    reader(filename);
  O    buf;

#ifndef DEBUG
  read(reader, "testobject", buf);
  bool good = (object == buf);
  std::cout << name << " IO test: " << std::endl;
  if (!good) exit(EXIT_FAILURE);
#endif
}

#ifdef DEBUG
//typedef int TestScalar;
typedef std::complex<double> TestScalar;
typedef Eigen::Tensor<iMatrix<TestScalar,1>, 6> TestTensorSingle;
typedef Eigen::Tensor<TestScalar, 3, Eigen::StorageOptions::RowMajor> TestTensor;
typedef Eigen::TensorFixedSize<TestScalar, Eigen::Sizes<9,4,2>, Eigen::StorageOptions::RowMajor> TestTensorFixed;
typedef std::vector<TestTensorFixed> aTestTensorFixed;
typedef Eigen::TensorFixedSize<SpinColourVector, Eigen::Sizes<11,3,2>> LSCTensor;
typedef Eigen::TensorFixedSize<LorentzColourMatrix, Eigen::Sizes<5,7,2>> LCMTensor;
// From Test_serialisation.cc
class ETSerClass: Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(ETSerClass
                                  , SpinColourVector, scv
                                  , SpinColourMatrix, scm
                                  , TestTensor, Critter
                                  , TestTensorFixed, FixedCritter
                                  , aTestTensorFixed, aFixedCritter
                                  , LSCTensor, MyLSCTensor
                                  , LCMTensor, MyLCMTensor
                                  );
  ETSerClass() : Critter(7,3,2), aFixedCritter(3) {}
};

bool EigenIOTest(void) {
  constexpr TestScalar Inc{1,-1};
  TestTensorSingle ts(1,1,1,1,1,1);
  ts(0,0,0,0,0,0) = Inc * 3.1415927;
  ioTest<Hdf5Writer, Hdf5Reader, TestTensorSingle>("iotest_single.h5", ts, "Singlet");

  SpinColourVector scv, scv2;
  scv2 = scv;
  ioTest<Hdf5Writer, Hdf5Reader, SpinColourVector>("iotest_vector.h5", scv, "SpinColourVector");
  SpinColourMatrix scm;
  ioTest<Hdf5Writer, Hdf5Reader, SpinColourMatrix>("iotest_matrix.h5", scm, "SpinColourMatrix");
  
  TestTensor t(6,3,2);
  TestScalar Val{Inc};
  for( int i = 0 ; i < t.dimension(0) ; i++)
    for( int j = 0 ; j < t.dimension(1) ; j++)
      for( int k = 0 ; k < t.dimension(2) ; k++) {
        t(i,j,k) = Val;
        Val += Inc;
      }
  ioTest<Hdf5Writer, Hdf5Reader, TestTensor>("iotest_tensor.h5", t, "eigen_tensor_instance_name");

  // Now serialise a fixed size tensor
  using FixedTensor = Eigen::TensorFixedSize<TestScalar, Eigen::Sizes<8,4,3>>;
  FixedTensor tf;
  Val = Inc;
  for( int i = 0 ; i < tf.dimension(0) ; i++)
    for( int j = 0 ; j < tf.dimension(1) ; j++)
      for( int k = 0 ; k < tf.dimension(2) ; k++) {
        tf(i,j,k) = Val;
        Val += Inc;
      }
  ioTest<Hdf5Writer, Hdf5Reader, FixedTensor>("iotest_tensor_fixed.h5", tf, "eigen_tensor_fixed_name");
  
  ETSerClass o;
  ioTest<Hdf5Writer, Hdf5Reader, ETSerClass>("iotest_object.h5", o, "ETSerClass_object_instance_name");
  
  // Tensor of spin colour
  LSCTensor l;
  Val = 0;
  for( int i = 0 ; i < l.dimension(0) ; i++)
    for( int j = 0 ; j < l.dimension(1) ; j++)
      for( int k = 0 ; k < l.dimension(2) ; k++)
        for( int s = 0 ; s < Ns ; s++ )
          for( int c = 0 ; c < Nc ; c++ )
          {
            l(i,j,k)()(s)(c) = Val;
            Val += Inc;
          }
  ioTest<Hdf5Writer, Hdf5Reader, LSCTensor>("iotest_LSCTensor.h5", l, "LSCTensor_object_instance_name");
  std::cout << "t:";
  for_all( l, [](std::complex<double> &c, LSCTensor::Index index, const std::size_t * pDims ){
    std::cout << "  ";
    for( int i = 0 ; i < 5; i++ )
      std::cout << "[" << pDims[i] << "]";
    std::cout << " = " << c << std::endl;
  } );
  std::cout << std::endl;

  // Tensor of spin colour
  LCMTensor l2;
  Val = 0;
  for( int i = 0 ; i < l2.dimension(0) ; i++)
    for( int j = 0 ; j < l2.dimension(1) ; j++)
      for( int k = 0 ; k < l2.dimension(2) ; k++)
        for( int l = 0 ; l < Ns ; l++ )
          for( int c = 0 ; c < Nc ; c++ )
            for( int c2 = 0 ; c2 < Nc ; c2++ )
            {
              l2(i,j,k)(l)()(c,c2) = Val;
              Val += Inc;
            }
  ioTest<Hdf5Writer, Hdf5Reader, LCMTensor>("iotest_LCMTensor.h5", l2, "LCMTensor_object_instance_name");
  
  std::cout << "Wow!" << std::endl;
  
  return true;
}
#endif

template <typename T>
void tensorConvTestFn(GridSerialRNG &rng, const std::string label)
{
  T    t, ft;
  Real n;
  bool good;

  random(rng, t);
  auto tv = tensorToVec(t);
  vecToTensor(ft, tv);
  n    = norm2(t - ft);
  good = (n == 0);
  std::cout << label << " norm 2 diff: " << n << " -- " 
            << (good ? "success" : "failure") << std::endl;
}

#define tensorConvTest(rng, type) tensorConvTestFn<type>(rng, #type)

int main(int argc,char **argv)
{
  Grid_init(&argc,&argv);
  std::cout << std::boolalpha << "==== basic IO" << std::endl; // display true / false for boolean

#ifndef DEBUG
  GridSerialRNG    rng;

  rng.SeedFixedIntegers(std::vector<int>({42,10,81,9}));

  XmlWriter WR("bother.xml");

  // test basic type writing
  std::cout << "-- basic writing to 'bother.xml'..." << std::endl;
  push(WR,"BasicTypes");
  write(WR,std::string("i16"),i16);
  write(WR,"u16",u16);
  write(WR,"i32",i32);
  write(WR,"u32",u32);
  write(WR,"i64",i64);
  write(WR,"u64",u64);
  write(WR,"f",f);
  write(WR,"d",d);
  write(WR,"b",b);
  pop(WR);

  // test serializable class writing
  myclass              obj(1234); // non-trivial constructor
  std::vector<myclass> vec;
  std::pair<myenum, myenum> pair;

  std::cout << "-- serialisable class writing to 'bother.xml'..." << std::endl;
  write(WR,"obj",obj);
  WR.write("obj2", obj);
  vec.push_back(obj);
  vec.push_back(myclass(5678));
  vec.push_back(myclass(3838));
  pair = std::make_pair(myenum::red, myenum::blue);

  write(WR, "objvec", vec);
  std::cout << "-- serialisable class writing to std::cout:" << std::endl;
  std::cout << obj << std::endl;
  std::cout << "-- serialisable class comparison:" << std::endl;
  std::cout << "vec[0] == obj: " << (vec[0] == obj) << std::endl;
  std::cout << "vec[1] == obj: " << (vec[1] == obj) << std::endl;
  std::cout << "-- pair writing to std::cout:" << std::endl;
  std::cout << pair << std::endl;

  // read tests
  std::cout << "\n==== IO self-consistency tests" << std::endl;
  //// XML
  ioTest<XmlWriter, XmlReader>("iotest.xml", obj, "XML    (object)           ");
  ioTest<XmlWriter, XmlReader>("iotest.xml", vec, "XML    (vector of objects)");
  //// binary
  ioTest<BinaryWriter, BinaryReader>("iotest.bin", obj, "binary (object)           ");
  ioTest<BinaryWriter, BinaryReader>("iotest.bin", vec, "binary (vector of objects)");
  //// text
  ioTest<TextWriter, TextReader>("iotest.dat", obj, "text   (object)           ");
  ioTest<TextWriter, TextReader>("iotest.dat", vec, "text   (vector of objects)");
  //// text
  ioTest<JSONWriter, JSONReader>("iotest.json", obj,  "JSON   (object)           ");
  ioTest<JSONWriter, JSONReader>("iotest.json", vec,  "JSON   (vector of objects)");

  //// HDF5
#ifdef HAVE_HDF5
  ioTest<Hdf5Writer, Hdf5Reader>("iotest.h5", obj, "HDF5   (object)           ");
  ioTest<Hdf5Writer, Hdf5Reader>("iotest.h5", vec, "HDF5   (vector of objects)");
#endif

  std::cout << "\n==== vector flattening/reconstruction" << std::endl;
  typedef std::vector<std::vector<std::vector<double>>> vec3d;

  vec3d dv, buf;
  double d = 0.;

  dv.resize(4);
  for (auto &v1: dv)
  {
    v1.resize(3);
    for (auto &v2: v1)
    {
      v2.resize(5);
      for (auto &x: v2)
      {
        x = d++;
      }
    }
  }
  std::cout << "original 3D vector:" << std::endl;
  std::cout << dv << std::endl;

  Flatten<vec3d> flatdv(dv);

  std::cout << "\ndimensions:" << std::endl;
  std::cout << flatdv.getDim() << std::endl;
  std::cout << "\nflattened vector:" << std::endl;
  std::cout << flatdv.getFlatVector() << std::endl;

  Reconstruct<vec3d> rec(flatdv.getFlatVector(), flatdv.getDim());
  std::cout << "\nreconstructed vector:" << std::endl;
  std::cout << flatdv.getVector() << std::endl;
  std::cout << std::endl;

  std::cout << "==== Grid tensor to vector test" << std::endl;
  tensorConvTest(rng, SpinColourMatrix);
  tensorConvTest(rng, SpinColourVector);
  tensorConvTest(rng, ColourMatrix);
  tensorConvTest(rng, ColourVector);
  tensorConvTest(rng, SpinMatrix);
  tensorConvTest(rng, SpinVector);
#elif HAVE_HDF5
  if(! EigenIOTest() ) exit(EXIT_FAILURE);
#endif
  Grid_finalize();
}
