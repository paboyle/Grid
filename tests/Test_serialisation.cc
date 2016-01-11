    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_serialisation.cc

    Copyright (C) 2015

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
#include <Grid.h>

namespace Grid {
  
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
                            std::vector<std::vector<double>>, twodimarray,
                            );
    
    myclass() {}
    myclass(int i)
    : array(4,5.1), twodimarray(3,std::vector<double>(2,1.23456)), ve(2, myenum::blue)
    {
      e=myenum::red;
      x=i;
      y=2*i;
      b=true;
      name="bother said pooh";
    }
  };
  
}

using namespace Grid;

int16_t i16 = 1;
uint16_t u16 = 2;
int32_t i32 = 3;
uint32_t u32 = 4;
int64_t i64 = 5;
uint64_t u64 = 6;
float    f = M_PI;
double   d = 2*M_PI;
bool     b = false;

int main(int argc,char **argv)
{
  {
    XmlWriter WR("bother.xml");
    
    // test basic type writing
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
    myclass obj(1234); // non-trivial constructor
    write(WR,"obj",obj);
    WR.write("obj2", obj);
    std::cout << obj << std::endl;
    
    std::vector<myclass> vec;
    vec.push_back(myclass(1234));
    vec.push_back(myclass(5678));
    vec.push_back(myclass(3838));
    write(WR, "objvec", vec);
  };
  
  // read tests
  myclass copy1, copy2, copy3;
  std::vector<myclass> veccopy1, veccopy2, veccopy3;
  //// XML
  {
    XmlReader RD("bother.xml");
    read(RD,"obj",copy1);
    read(RD,"objvec", veccopy1);
    std::cout << "Loaded (XML) -----------------" << std::endl;
    std::cout << copy1 << std::endl << veccopy1 << std::endl;
  }
  //// binary
  {
    BinaryWriter BWR("bother.bin");
    write(BWR,"discard",copy1 );
    write(BWR,"discard",veccopy1 );
  }
  {
    BinaryReader BRD("bother.bin");
    read (BRD,"discard",copy2 );
    read (BRD,"discard",veccopy2 );
    std::cout << "Loaded (bin) -----------------" << std::endl;
    std::cout << copy2 << std::endl << veccopy2 << std::endl;
  }
  //// text
  {
    TextWriter TWR("bother.txt");
    write(TWR,"discard",copy1 );
    write(TWR,"discard",veccopy1 );
  }
  {
    TextReader TRD("bother.txt");
    read (TRD,"discard",copy3 );
    read (TRD,"discard",veccopy3 );
    std::cout << "Loaded (txt) -----------------" << std::endl;
    std::cout << copy3 << std::endl << veccopy3 << std::endl;
  }
}
