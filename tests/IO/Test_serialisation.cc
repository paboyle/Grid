    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/Test_serialisation.cc

    Copyright (C) 2015-2019

Author: Guido Cossu <guido.cossu@ed.ac.uk>
Author: Antonin Portelli <antonin.portelli@me.com>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Michael Marshall <michael.marshall@ed.ac.uk>

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
#include <Grid/util/EigenUtil.h>

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
void ioTest(const std::string &filename, const O &object, const std::string &name,
            const char * tag = "testobject", unsigned short Precision = 0 )
{
  std::cout << "IO test: " << name << " -> " << filename << " ...";
  // writer needs to be destroyed so that writing physically happens
  {
    W writer(filename);
    if( Precision )
      writer.setPrecision(Precision);
    write(writer, tag , object);
  }

  std::cout << " done. reading ...";
  R    reader(filename);
  std::unique_ptr<O> buf( new O ); // In case object too big for stack

  read(reader, tag, *buf);
  bool good = Serializable::CompareMember(object, *buf);
  if (!good) {
    std::cout << " failure!" << std::endl;
    if (EigenIO::is_tensor<O>::value)
      dump_tensor(*buf);
    exit(EXIT_FAILURE);
  }
  std::cout << " done." << std::endl;
}

typedef ComplexD TestScalar;
typedef Eigen::TensorFixedSize<unsigned short, Eigen::Sizes<5,4,3,2,1>> TensorRank5UShort;
typedef Eigen::TensorFixedSize<unsigned short, Eigen::Sizes<5,4,3,2,1>, Eigen::StorageOptions::RowMajor> TensorRank5UShortAlt;
typedef Eigen::Tensor<TestScalar, 3, Eigen::StorageOptions::RowMajor> TensorRank3;
typedef Eigen::TensorFixedSize<TestScalar, Eigen::Sizes<9,4,2>, Eigen::StorageOptions::RowMajor> Tensor_9_4_2;
typedef std::vector<Tensor_9_4_2> aTensor_9_4_2;
typedef Eigen::TensorFixedSize<SpinColourVector, Eigen::Sizes<6,5>> LSCTensor;

class PerambIOTestClass: Serializable {
  ComplexD Flag;
public:
  using PerambTensor = Eigen::Tensor<SpinColourVector, 6, Eigen::StorageOptions::RowMajor>;
  GRID_SERIALIZABLE_CLASS_MEMBERS(PerambIOTestClass
                                  , SpinColourVector, spinColourVector
                                  , SpinColourMatrix, spinColourMatrix
                                  , std::vector<std::string>, DistilParameterNames
                                  , std::vector<int>,         DistilParameterValues
                                  , PerambTensor,             Perambulator
                                  , PerambTensor,             Perambulator2
                                  , TensorRank5UShort,        tensorRank5UShort
                                  , TensorRank3,              tensorRank3
                                  , Tensor_9_4_2,             tensor_9_4_2
                                  , aTensor_9_4_2,            atensor_9_4_2
                                  , LSCTensor,                MyLSCTensor
                                  );
  PerambIOTestClass()
  : DistilParameterNames {"do", "androids", "dream", "of", "electric", "sheep?"}
  , DistilParameterValues{2,3,1,4,5,1}
  , Perambulator(2,3,1,4,5,1)
  , Perambulator2(7,1,6,1,5,1)
  , tensorRank3(7,3,2)
  , atensor_9_4_2(3)
  //, Flag(1,-3.1415927)
  , Flag(1,-1)
  {
    SequentialInit(Perambulator,  Flag);
    SequentialInit(Perambulator2, Flag);
    SequentialInit(tensorRank5UShort);
    SequentialInit(tensorRank3, Flag);
    SequentialInit(tensor_9_4_2, Flag);
    for( auto &t : atensor_9_4_2 ) SequentialInit(t, Flag);
    SequentialInit( MyLSCTensor, Flag );
  }
};

#define TEST_PARAMS( T ) #T, Flag, Precision, filename, pszExtension, TestNum

// Perform an I/O test for a single Eigen tensor (of any type)
template <typename WTR_, typename RDR_, typename T, typename... IndexTypes>
void EigenTensorTestSingle(const char * MyTypeName, typename EigenIO::Traits<T>::scalar_type Flag,
                           unsigned short Precision, std::string &filename, const char * pszExtension,
                           unsigned int &TestNum, IndexTypes... otherDims)
{
  using Traits = EigenIO::Traits<T>;
  using scalar_type = typename Traits::scalar_type;
  std::unique_ptr<T> pTensor{new T(otherDims...)};
  SequentialInit( * pTensor, Flag, Precision );
  filename = "iotest_" + std::to_string(++TestNum) + "_" + MyTypeName + pszExtension;
  ioTest<WTR_, RDR_, T>(filename, * pTensor, MyTypeName, MyTypeName);
}

// Perform a series of I/O tests for Eigen tensors, including a serialisable object
template <typename WTR_, typename RDR_>
void EigenTensorTest(const char * pszExtension, unsigned short Precision = 0)
{
  // Perform a series of tests on progressively more complex tensors
  unsigned int TestNum = 0;
  std::string filename;
  {
    int Flag = 7;
    using TensorSingle = Eigen::TensorFixedSize<Integer, Eigen::Sizes<1>>;
    EigenTensorTestSingle<WTR_, RDR_, TensorSingle>(TEST_PARAMS( TensorSingle ));
  }
  TestScalar Flag{1,-3.1415927};
  using TensorSimple = Eigen::Tensor<iMatrix<TestScalar,1>, 6>;
  using I = typename TensorSimple::Index;
  EigenTensorTestSingle<WTR_, RDR_, TensorSimple, I, I, I, I, I, I>( TEST_PARAMS( TensorSimple ), 1, 1, 1, 1, 1, 1 );
  EigenTensorTestSingle<WTR_, RDR_, TensorRank3, I, I, I>( TEST_PARAMS( TensorRank3 ), 6, 3, 2 );
  EigenTensorTestSingle<WTR_, RDR_, Tensor_9_4_2>(TEST_PARAMS( Tensor_9_4_2 ));
  EigenTensorTestSingle<WTR_, RDR_, LSCTensor>(TEST_PARAMS( LSCTensor ));
  // Now see whether we could write out a tensor in one memory order and read back in the other
  {
    unsigned short Flag = 1;
    EigenTensorTestSingle<WTR_, RDR_, TensorRank5UShort>(TEST_PARAMS( TensorRank5UShort ));
    std::cout << "    Testing alternate memory order read ... ";
    TensorRank5UShortAlt t2;
    RDR_ reader(filename);
    read(reader, "TensorRank5UShort", t2);
    bool good = true;
    using Index = typename TensorRank5UShortAlt::Index;
    // NB: I can't call
    for_all( t2, [&](unsigned short c, Index n,
                     const std::array<Index, TensorRank5UShortAlt::NumIndices> &TensorIndex,
                     const std::array<int, EigenIO::Traits<TensorRank5UShortAlt>::Rank> &GridTensorIndex ){
      good = good && ( c == n );
    } );
    if (!good) {
      std::cout << " failure!" << std::endl;
      dump_tensor(t2,"t2");
      exit(EXIT_FAILURE);
    }
    std::cout << " done." << std::endl;
  }
  // Now test a serialisable object containing a number of tensors
  {
    static const char MyTypeName[] = "PerambIOTestClass";
    std::unique_ptr<PerambIOTestClass> pObj{new PerambIOTestClass()};
    filename = "iotest_" + std::to_string(++TestNum) + "_" + MyTypeName + pszExtension;
    ioTest<WTR_, RDR_, PerambIOTestClass>(filename, * pObj, MyTypeName, MyTypeName);
  }
  // Stress test. Too large for the XML or text readers and writers!
#ifdef STRESS_TEST
  if( typeid( WTR_ ).name() == typeid( Hdf5Writer ).name() || typeid( WTR_ ).name() == typeid( BinaryWriter ).name() ) {
    using LCMTensor=Eigen::TensorFixedSize<iMatrix<iVector<iMatrix<iVector<LorentzColourMatrix,5>,2>,7>,3>,
        Eigen::Sizes<2,4,11,10,9>, Eigen::StorageOptions::RowMajor>;
    std::cout << "sizeof( LCMTensor ) = " << sizeof( LCMTensor ) / 1024 / 1024 << " MB" << std::endl;
    EigenTensorTestSingle<WTR_, RDR_, LCMTensor>(TEST_PARAMS( LCMTensor ));
  }
#endif
}

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

  std::cout << "-- serialisable class writing to 'bother.xml'..." << std::endl;
  write(WR,"obj",obj);
  WR.write("obj2", obj);
  vec.push_back(obj);
  vec.push_back(myclass(5678));
  vec.push_back(myclass(3838));

  write(WR, "objvec", vec);
  std::cout << "-- serialisable class writing to std::cout:" << std::endl;
  std::cout << obj << std::endl;
  std::cout << "-- serialisable class comparison:" << std::endl;
  std::cout << "vec[0] == obj: " << (vec[0] == obj) << std::endl;
  std::cout << "vec[1] == obj: " << (vec[1] == obj) << std::endl;
  std::cout << "-- pair writing to std::cout:" << std::endl;
  std::pair<myenum, myenum> pair = std::make_pair(myenum::red, myenum::blue);
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
  std::cout << "\n==== detailed Hdf5 tensor tests (Grid::EigenIO)" << std::endl;
  EigenTensorTest<Hdf5Writer, Hdf5Reader>(".h5");
#endif
  std::cout << "\n==== detailed binary tensor tests (Grid::EigenIO)" << std::endl;
  EigenTensorTest<BinaryWriter, BinaryReader>(".bin");
  std::cout << "\n==== detailed xml tensor tests (Grid::EigenIO)" << std::endl;
  EigenTensorTest<XmlWriter, XmlReader>(".xml", 6);
  std::cout << "\n==== detailed text tensor tests (Grid::EigenIO)" << std::endl;
  EigenTensorTest<TextWriter, TextReader>(".dat", 5);

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

  Grid_finalize();
}
