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
                          std::vector<std::vector<std::vector<Complex> > >, cmplx3darray
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
  bool good;
  
  read(reader, "testobject", buf);
  good = (object == buf);
  std::cout << name << " IO test: " << (good ? "success" : "failure");
  std::cout << std::endl;
  if (!good) exit(EXIT_FAILURE);
}

int main(int argc,char **argv)
{
  std::cout << "==== basic IO" << std::endl;
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
  vec.push_back(myclass(1234));
  vec.push_back(myclass(5678));
  vec.push_back(myclass(3838));
  pair = std::make_pair(myenum::red, myenum::blue);

  write(WR, "objvec", vec);
  std::cout << "-- serialisable class writing to std::cout:" << std::endl;
  std::cout << obj << std::endl;
  std::cout << "-- serialisable class comparison:" << std::endl;
  std::cout << "vec[0] == obj: " << ((vec[0] == obj) ? "true" : "false") << std::endl;
  std::cout << "vec[1] == obj: " << ((vec[1] == obj) ? "true" : "false") << std::endl;
  
  write(WR, "objpair", pair);
  std::cout << "-- pair writing to std::cout:" << std::endl;
  std::cout << pair << std::endl;
  
  // read tests
  std::cout << "\n==== IO self-consistency tests" << std::endl;
  //// XML
  ioTest<XmlWriter, XmlReader>("iotest.xml", obj, "XML    (object)           ");
  ioTest<XmlWriter, XmlReader>("iotest.xml", vec, "XML    (vector of objects)");
  ioTest<XmlWriter, XmlReader>("iotest.xml", pair, "XML    (pair of objects)");
  //// binary
  ioTest<BinaryWriter, BinaryReader>("iotest.bin", obj, "binary (object)           ");
  ioTest<BinaryWriter, BinaryReader>("iotest.bin", vec, "binary (vector of objects)");
  ioTest<BinaryWriter, BinaryReader>("iotest.bin", pair, "binary (pair of objects)");
  //// text
  ioTest<TextWriter, TextReader>("iotest.dat", obj, "text   (object)           ");
  ioTest<TextWriter, TextReader>("iotest.dat", vec, "text   (vector of objects)");
  ioTest<TextWriter, TextReader>("iotest.dat", pair, "text   (pair of objects)");
  //// HDF5
#undef HAVE_HDF5
#ifdef HAVE_HDF5
  ioTest<Hdf5Writer, Hdf5Reader>("iotest.h5", obj, "HDF5   (object)           ");
  ioTest<Hdf5Writer, Hdf5Reader>("iotest.h5", vec, "HDF5   (vector of objects)");
  ioTest<Hdf5Writer, Hdf5Reader>("iotest.h5", pair, "HDF5   (pair of objects)");
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


  std::cout << ".:::::: Testing JSON classes "<< std::endl;


  {
    JSONWriter JW("bother.json");
    
    // test basic type writing
    push(JW,"BasicTypes");
    write(JW,std::string("i16"),i16);
    write(JW,"u16",u16);
    write(JW,"i32",i32);
    write(JW,"u32",u32);
    write(JW,"i64",i64);
    write(JW,"u64",u64);
    write(JW,"f",f);
    write(JW,"d",d);
    write(JW,"b",b);
    pop(JW);
    
    // test serializable class writing
    myclass obj(1234); // non-trivial constructor
    std::cout << "-- serialisable class writing to 'bother.json'..." << std::endl;
    write(JW,"obj",obj);
    JW.write("obj2", obj);
    
    std::cout << obj << std::endl;
    
    std::vector<myclass> vec;
    vec.push_back(myclass(1234));
    vec.push_back(myclass(5678));
    vec.push_back(myclass(3838));
    write(JW, "objvec", vec);
    
  }

  {
    JSONReader RD("bother.json");
    myclass jcopy1;
    std::vector<myclass> jveccopy1;
    read(RD,"obj",jcopy1);
    read(RD,"objvec", jveccopy1);
    std::cout << "Loaded (JSON) -----------------" << std::endl;
    std::cout << jcopy1 << std::endl << jveccopy1 << std::endl;
  }
  
  { 
    ildgFormat format;
    format.version   =1.0;
    format.field     =std::string("su3gauge");
    format.precision =32;
    format.lx        =24;
    format.ly        =24;
    format.lz        =24;
    format.lt        =48;
    XmlWriter WR("ildg-format.xml","");
    XmlWriter WRs("","");
    write(WR,"ildgFormat",format);
    write(WRs,"ildgFormat",format);
    std::cout << " XmlString: " <<WRs.XmlString()<<std::endl;
  }
/* 
  // This is still work in progress
  {
    // Testing the next element function
    JSONReader RD("test.json");
    RD.push("grid");
    RD.push("Observable");
    std::string name;
    read(RD,"name", name);
  }
*/


}
