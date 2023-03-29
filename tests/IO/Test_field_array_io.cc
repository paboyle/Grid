    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/IO/Test_field_array_io.cc

    Copyright (C) 2015

Author: Christopher Kelly <ckelly@bnl.gov>
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

using namespace std;
using namespace Grid;

//This test demonstrates and checks a single-file write of an arbitrary array of fields

uint64_t writeHeader(const uint32_t size, const uint32_t checksum, const std::string &format, const std::string &file){
  std::ofstream fout(file,std::ios::out|std::ios::in);
  fout.seekp(0,std::ios::beg);
  fout << std::setw(10) << size << std::endl;
  fout << std::hex << std::setw(10) << checksum << std::endl;
  fout << format << std::endl;
  return fout.tellp();
}
 
uint64_t readHeader(uint32_t &size, uint32_t &checksum, std::string &format, const std::string &file){
  std::ifstream fin(file);
  std::string line;
  getline(fin,line);
  {
    std::stringstream ss; ss <<line ; ss >> size;
  }
  getline(fin,line);
  {
    std::stringstream ss; ss <<line ; ss >> std::hex >> checksum;
  }
  getline(fin,format);
  removeWhitespace(format);
      
  return fin.tellg();
}
 
template<typename FieldType>
void writeFieldArray(const std::string &file, const std::vector<FieldType> &data){
  typedef typename FieldType::vector_object vobj;
  typedef typename FieldType::scalar_object sobj;
  GridBase* grid = data[0].Grid(); //assume all fields have the same Grid
  BinarySimpleMunger<sobj, sobj> munge; //straight copy

  //We need a 2-pass header write, first to establish the size, the second pass writes the checksum
  std::string format = getFormatString<typename FieldType::vector_object>();

  uint64_t offset; //leave 64 bits for header
  if ( grid->IsBoss() ) { 
    NerscIO::truncate(file);
    offset = writeHeader(data.size(), 0, format, file);
  }
  grid->Broadcast(0,(void *)&offset,sizeof(offset)); //use as a barrier

  std::cout << "Data offset write " << offset << std::endl;
  std::cout << "Data size write " << data.size() << std::endl;
  uint64_t field_size = uint64_t(grid->gSites()) * sizeof(sobj);
  std::cout << "Field size = " << field_size << " B" << std::endl;

  uint32_t checksum = 0;
  for(int i=0;i<data.size();i++){
    std::cout << "Data field write " << i << " offset " << offset << std::endl;
    uint32_t nersc_csum,scidac_csuma,scidac_csumb;
    BinaryIO::writeLatticeObject<vobj,sobj>(const_cast<FieldType &>(data[i]),file,munge,offset,format,
					    nersc_csum,scidac_csuma,scidac_csumb);
    offset += field_size;
    checksum ^= nersc_csum + 0x9e3779b9 + (checksum<<6) + (checksum>>2);
  }
  std::cout << "Write checksum " << checksum << std::endl;

  if ( grid->IsBoss() ) { 
    writeHeader(data.size(), checksum, format, file);
  }
}


template<typename FieldType>
void readFieldArray(std::vector<FieldType> &data, const std::string &file){
  typedef typename FieldType::vector_object vobj;
  typedef typename FieldType::scalar_object sobj;
  assert(data.size() > 0);
  GridBase* grid = data[0].Grid(); //assume all fields have the same Grid
  BinarySimpleUnmunger<sobj, sobj> munge; //straight copy
  
  uint32_t hdr_checksum, hdr_size;
  std::string format;
  uint64_t offset = readHeader(hdr_size, hdr_checksum, format, file);
  
  std::cout << "Data offset read " << offset << std::endl;  
  std::cout << "Data size read " << hdr_size << std::endl;
  assert(data.size() == hdr_size);

  uint64_t field_size = uint64_t(grid->gSites()) * sizeof(sobj);

  uint32_t checksum = 0;

  for(int i=0;i<data.size();i++){
    std::cout << "Data field read " << i << " offset " << offset << std::endl;
    uint32_t nersc_csum,scidac_csuma,scidac_csumb;
    BinaryIO::readLatticeObject<vobj,sobj>(data[i],file,munge,offset,format,
					   nersc_csum,scidac_csuma,scidac_csumb);
    offset += field_size;
    checksum ^= nersc_csum + 0x9e3779b9 + (checksum<<6) + (checksum>>2);
  }

  std::cout << "Header checksum " << hdr_checksum << std::endl;    
  std::cout << "Read checksum " << checksum << std::endl;
    

  assert( hdr_checksum == checksum );
}




int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Coordinate latt   = GridDefaultLatt();
  Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();

  const int Ls=8;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(latt, simd_layout, mpi_layout);
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  typedef DomainWallFermionD::FermionField FermionField;

  int nfield = 20;
  std::vector<FermionField> data(nfield, FGrid);

  for(int i=0;i<data.size();i++)
    gaussian(RNG5, data[i]);
  
  std::string file = "test_field_array_io.0";
  writeFieldArray(file, data);

  std::vector<FermionField> data_r(nfield, FGrid);
  readFieldArray(data_r, file);
  
  for(int i=0;i<nfield;i++){
    FermionField diff = data_r[i] - data[i];
    RealD norm_diff = norm2(diff);
    std::cout << "Norm2 of difference between stored and loaded data index " << i << " : " << norm_diff << std::endl;
  }
  
  std::cout << "Done" << std::endl;

  Grid_finalize();
}
