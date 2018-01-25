    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_cshift.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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
 ;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();

  GridCartesian        Fine(latt_size,simd_layout,mpi_layout);

  GridParallelRNG      FineRNG(&Fine);  FineRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));

  LatticeComplex U(&Fine);
  LatticeComplex ShiftU(&Fine);

  LatticeComplex lex(&Fine);
  lex=zero;
  Integer stride =1;
  {
    LatticeComplex coor(&Fine);

    for(int d=0;d<4;d++){
      LatticeCoordinate(coor,d);
      lex = lex + coor*stride;
      stride=stride*latt_size[d];
    }
    U=lex;
  }

  std::stringstream ss;
  ss<<"error";
  for(int d=0;d<Fine._ndimension;d++){
    ss<<"."<<Fine._processor_coor[d];
  }
  ss<<"_wr_"<<Fine._processor;
  std::string fname(ss.str());
  std::ofstream ferr(fname);
  
  TComplex cm;
  
  for(int dir=0;dir<4;dir++){
    for(int shift=0;shift<latt_size[dir];shift++){
      if ( Fine.IsBoss() ) 
	std::cout<<GridLogMessage<<"Shifting by "<<shift<<" in direction"<<dir<<std::endl;

	ShiftU  = Cshift(U,dir,shift);    // Shift everything

	std::vector<int> coor(4);

	for(coor[3]=0;coor[3]<latt_size[3];coor[3]++){
	for(coor[2]=0;coor[2]<latt_size[2];coor[2]++){
	for(coor[1]=0;coor[1]<latt_size[1];coor[1]++){
	for(coor[0]=0;coor[0]<latt_size[0];coor[0]++){
	  
	  peekSite(cm,ShiftU,coor);

	  double nrm=norm2(U);

	  std::vector<int> scoor(coor);
	  scoor[dir] = (scoor[dir]+shift)%latt_size[dir];
	  
	  Integer slex = scoor[0]
	    + latt_size[0]*scoor[1]
	    + latt_size[0]*latt_size[1]*scoor[2]
	    + latt_size[0]*latt_size[1]*latt_size[2]*scoor[3];

	  Complex scm(slex);
	  
	  nrm = abs(scm-cm()()());
	  std::vector<int> peer(4);
	  Complex tmp  =cm;
	  Integer index=real(tmp);
	  Lexicographic::CoorFromIndex(peer,index,latt_size);

	  if (nrm > 0){


	    ferr<<"FAIL shift "<< shift<<" in dir "<< dir<<" ["<<coor[0]<<","<<coor[1]<<","<<coor[2]<<","<<coor[3]<<"] = "<< cm()()()<<" expect "<<scm<<"  "<<nrm<<std::endl;
	    ferr<<"Got    "<<index<<" " << peer[0]<<","<<peer[1]<<","<<peer[2]<<","<<peer[3]<<std::endl;
	    index=real(scm);
	    Lexicographic::CoorFromIndex(peer,index,latt_size);
	    ferr<<"Expect "<<index<<" " << peer[0]<<","<<peer[1]<<","<<peer[2]<<","<<peer[3]<<std::endl;
	  }
	}}}}
    }
  }

  Grid_finalize();
}
