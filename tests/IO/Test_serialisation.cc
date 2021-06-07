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
#include <typeinfo>

using namespace Grid;

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
			  std::vector<std::vector<std::vector<std::complex<double>> > >, cmplx3darray,
                          std::vector<std::vector<std::vector<int> > >, ragged,
			  SpinColourMatrix, scm
                          );
  myclass() {}
  myclass(int i)
  : array(4,5.1)
  , twodimarray(3,std::vector<double>(5, 1.23456))
  , cmplx3darray(3,std::vector<std::vector<std::complex<double>>>(5, std::vector<std::complex<double>>(7, std::complex<double>(1.2, 3.4))))
  , ve(2, myenum::blue)
  , ragged( {{{i+1},{i+2,i+3}}, // ragged
            {{i+4,i+5,i+6,i+7},{i+8,i+9,i+10,i+11},{i+12,i+13,i+14,i+15}}, // block
            {{i+16,i+17},{i+18,i+19,i+20}}} ) //ragged
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
    exit(EXIT_FAILURE);
  }
  std::cout << " done." << std::endl;
}

// The only way I could get these iterators to work is to put the begin() and end() functions in the Eigen namespace
// So if Eigen ever defines these, we'll have a conflict and have to change this
namespace Eigen {
  template <typename ET>
  inline typename std::enable_if<EigenIO::is_tensor<ET>::value, typename EigenIO::Traits<ET>::scalar_type *>::type
  begin( ET & et ) { return reinterpret_cast<typename Grid::EigenIO::Traits<ET>::scalar_type *>(et.data()); }
  template <typename ET>
  inline typename std::enable_if<EigenIO::is_tensor<ET>::value, typename EigenIO::Traits<ET>::scalar_type *>::type
  end( ET & et ) { return begin(et) + et.size() * EigenIO::Traits<ET>::count; }
}

// Perform I/O tests on a range of tensor types
// Test coverage: scalars, complex and GridVectors in single, double and default precision
class TensorIO : public Serializable {
  using TestScalar = std::complex<double>;
  using SR3 = Eigen::Sizes<9,4,2>;
  using SR5 = Eigen::Sizes<5,4,3,2,1>;
  using ESO = Eigen::StorageOptions;
  using TensorRank3  = Eigen::Tensor<std::complex<float>, 3, ESO::RowMajor>;
  using TensorR5     = Eigen::TensorFixedSize<Real, SR5>;
  using TensorR5Alt  = Eigen::TensorFixedSize<Real, SR5, ESO::RowMajor>;
  using Tensor942    = Eigen::TensorFixedSize<TestScalar, SR3, ESO::RowMajor>;
  using aTensor942   = std::vector<Tensor942>;
  using Perambulator = Eigen::Tensor<SpinColourVector, 6, ESO::RowMajor>;
  using LSCTensor    = Eigen::TensorFixedSize<SpinColourMatrix, Eigen::Sizes<6,5>>;
  
  static const Real       FlagR;
  static const std::complex<double>   Flag;
  static const std::complex<float>    FlagF;
  static const TestScalar FlagTS;
  static const char * const pszFilePrefix;

  void Init(unsigned short Precision)
  {
    for( auto &s : Perambulator1 ) s = Flag;
    for( auto &s : Perambulator2 ) s = Flag;
    for( auto &s : tensorR5 )      s = FlagR;
    for( auto &s : tensorRank3 )   s = FlagF;
    for( auto &s : tensor_9_4_2 )  s = FlagTS;
    for( auto &t : atensor_9_4_2 )
      for( auto &s : t )           s = FlagTS;
    for( auto &s : MyLSCTensor )   s = Flag;
  }
  
  // Perform an I/O test for a single Eigen tensor (of any type)
  template <typename W, typename R, typename T, typename... IndexTypes>
  static void TestOne(const char * MyTypeName, unsigned short Precision, std::string &filename,
                      const char * pszExtension, unsigned int &TestNum,
                      typename EigenIO::Traits<T>::scalar_type Flag, IndexTypes... otherDims)
  {
    using Traits = EigenIO::Traits<T>;
    using scalar_type = typename Traits::scalar_type;
    std::unique_ptr<T> pTensor{new T(otherDims...)};
    for( auto &s : * pTensor ) s = Flag;
    filename = pszFilePrefix + std::to_string(++TestNum) + "_" + MyTypeName + pszExtension;
    ioTest<W, R, T>(filename, * pTensor, MyTypeName, MyTypeName);
  }
  
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(TensorIO
                                  , SpinColourVector,         spinColourVector
                                  , SpinColourMatrix,         spinColourMatrix
                                  , std::vector<std::string>, DistilParameterNames
                                  , std::vector<int>,         DistilParameterValues
                                  , Perambulator,             Perambulator1
                                  , Perambulator,             Perambulator2
                                  , TensorR5,                 tensorR5
                                  , TensorRank3,              tensorRank3
                                  , Tensor942,                tensor_9_4_2
                                  , aTensor942,               atensor_9_4_2
                                  , LSCTensor,                MyLSCTensor
                                  );
  TensorIO()
  : DistilParameterNames {"do", "androids", "dream", "of", "electric", "sheep?"}
  , DistilParameterValues{2,3,1,4,5,1}
  , Perambulator1(2,3,1,4,5,1)
  , Perambulator2(7,1,6,1,5,1)
  , tensorRank3(7,3,2)
  , atensor_9_4_2(3) {}
  
#define TEST_PARAMS( T ) #T, Precision, filename, pszExtension, TestNum
  
  // Perform a series of I/O tests for Eigen tensors, including a serialisable object
  template <typename WTR_, typename RDR_>
  static void Test(const char * pszExtension, unsigned short Precision = 0)
  {
    // Perform a series of tests on progressively more complex tensors
    unsigned int TestNum = 0;
    std::string filename;
    // Rank 1 tensor containing a single integer
    using TensorSingle = Eigen::TensorFixedSize<Integer, Eigen::Sizes<1>>;
    TestOne<WTR_, RDR_, TensorSingle>( TEST_PARAMS( TensorSingle ), 7 ); // lucky!
    // Rather convoluted way of defining four complex numbers
    using TensorSimple = Eigen::Tensor<iMatrix<TestScalar,2>, 6>;
    using I = typename TensorSimple::Index; // NB: Never specified, so same for all my test tensors
    // Try progressively more complicated tensors
    TestOne<WTR_, RDR_, TensorSimple, I,I,I,I,I,I>( TEST_PARAMS( TensorSimple ), FlagTS, 1,1,1,1,1,1 );
    TestOne<WTR_, RDR_, TensorRank3, I, I, I>( TEST_PARAMS( TensorRank3 ), FlagF, 6, 3, 2 );
    TestOne<WTR_, RDR_, Tensor942>(TEST_PARAMS( Tensor942 ), FlagTS);
    TestOne<WTR_, RDR_, LSCTensor>(TEST_PARAMS( LSCTensor ), Flag );
    TestOne<WTR_, RDR_, TensorR5>(TEST_PARAMS( TensorR5 ), FlagR);
    // Now test a serialisable object containing a number of tensors
    {
      static const char MyTypeName[] = "TensorIO";
      filename = pszFilePrefix + std::to_string(++TestNum) + "_" + MyTypeName + pszExtension;
      std::unique_ptr<TensorIO> pObj{new TensorIO()};
      pObj->Init(Precision);
      ioTest<WTR_, RDR_, TensorIO>(filename, * pObj, MyTypeName, MyTypeName, Precision);
    }
    // Stress test. Too large for the XML or text readers and writers!
#ifdef STRESS_TEST
    const std::type_info &tw = typeid( WTR_ );
    if( tw == typeid( Hdf5Writer ) || tw == typeid( BinaryWriter ) ) {
      using LCMTensor=Eigen::TensorFixedSize<iMatrix<iVector<iMatrix<iVector<LorentzColourMatrix,5>,2>,7>,3>,
      Eigen::Sizes<2,4,11,10,9>, Eigen::StorageOptions::RowMajor>;
      std::cout << "sizeof( LCMTensor ) = " << sizeof( LCMTensor ) / 1024 / 1024 << " MB" << std::endl;
      TestOne<WTR_, RDR_, LCMTensor>(TEST_PARAMS( LCMTensor ), Flag);
    }
#endif
  }
};

const Real                 TensorIO::FlagR {1};
const std::complex<double> TensorIO::Flag  {1,-1};
const std::complex<float>  TensorIO::FlagF {1,-1};
const TensorIO::TestScalar TensorIO::FlagTS{1,-1};
const char * const         TensorIO::pszFilePrefix = "tensor_";

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
  //  ioTest<JSONWriter, JSONReader>("iotest.json", obj,  "JSON   (object)           ");
  //  ioTest<JSONWriter, JSONReader>("iotest.json", vec,  "JSON   (vector of objects)");

  //// HDF5
#ifdef HAVE_HDF5
  ioTest<Hdf5Writer, Hdf5Reader>("iotest.h5", obj, "HDF5   (object)           ");
  ioTest<Hdf5Writer, Hdf5Reader>("iotest.h5", vec, "HDF5   (vector of objects)");
  std::cout << "\n==== detailed Hdf5 tensor tests (Grid::EigenIO)" << std::endl;
  TensorIO::Test<Hdf5Writer, Hdf5Reader>(".h5");
#endif
  std::cout << "\n==== detailed binary tensor tests (Grid::EigenIO)" << std::endl;
  TensorIO::Test<BinaryWriter, BinaryReader>(".bin");
  std::cout << "\n==== detailed xml tensor tests (Grid::EigenIO)" << std::endl;
  TensorIO::Test<XmlWriter, XmlReader>(".xml", 6);
  std::cout << "\n==== detailed text tensor tests (Grid::EigenIO)" << std::endl;
  TensorIO::Test<TextWriter, TextReader>(".dat", 5);

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

  {
    HMCparameters HMCparams;
    HMCparams.StartingType     =std::string("CheckpointStart");
    HMCparams.StartTrajectory  =7;
    HMCparams.Trajectories     =1000;
    HMCparams.NoMetropolisUntil=0;
    HMCparams.MD.name          =std::string("Force Gradient");
    HMCparams.MD.MDsteps       = 10;
    HMCparams.MD.trajL         = 1.0;

    XmlWriter HMCwr("HMCparameters.xml");
    write(HMCwr,"HMCparameters",HMCparams);
  }
  Grid_finalize();
}
