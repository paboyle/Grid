    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_cshift_red_black.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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

#define POWER10

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Coordinate latt_size   = GridDefaultLatt();
  int Nd = latt_size.size();
  Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();

  Coordinate mask(Nd,1);
  mask[0]=0;

  GridCartesian         Fine  (latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian RBFine(&Fine,mask,1);

  GridParallelRNG      FineRNG(&Fine);  FineRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));

  LatticeComplex U(&Fine);
  LatticeComplex ShiftU(&Fine);
  LatticeComplex rbShiftU(&Fine);
  LatticeComplex err(&Fine);
  LatticeComplex Ue(&RBFine); 
  LatticeComplex Uo(&RBFine);
  LatticeComplex ShiftUe(&RBFine);
  LatticeComplex ShiftUo(&RBFine);
  LatticeComplex lex(&Fine);
  lex=Zero();
  Integer stride =1;
  {
    LatticeComplex coor(&Fine);

    for(int d=0;d<Nd;d++){
      //      Integer i=10000;
      Integer i=0;
      LatticeCoordinate(coor,d);
      lex = lex + coor*stride+i;
#ifndef POWER10
      stride=stride*latt_size[d];
#else 
      stride=stride*10;
#endif
    }
    U=lex;
  }

  pickCheckerboard(Even,Ue,U);
  pickCheckerboard(Odd,Uo,U);

  //  std::cout<<GridLogMessage << U<<std::endl;
  std::cout<<GridLogMessage<< U <<std::endl;
  std::cout<<GridLogMessage << "U " <<norm2(U)<<std::endl;
  std::cout<<GridLogMessage << "Ue " <<norm2(Ue)<<std::endl;
  std::cout<<GridLogMessage << "Uo " <<norm2(Uo)<<std::endl;

  TComplex cm;
  TComplex cmeo;
  for(int dir=0;dir<Nd;dir++){
    //    if ( dir!=1 ) continue;
    for(int shift=0;shift<latt_size[dir];shift++){

      std::cout<<GridLogMessage<<"Shifting by "<<shift<<" in direction"<<dir<< "  ";

	//	std::cout<<GridLogMessage<<"Even grid"<<std::endl;
	ShiftUe = Cshift(Ue,dir,shift);    // Shift everything cb by cb
	//	std::cout<<GridLogMessage << "\tShiftUe " <<norm2(ShiftUe)<<std::endl;

	//	std::cout<<GridLogMessage<<"Odd grid"<<std::endl;
	ShiftUo = Cshift(Uo,dir,shift);    
	//	std::cout<<GridLogMessage << "\tShiftUo " <<norm2(ShiftUo)<<std::endl;

	//	std::cout<<GridLogMessage<<"Recombined Even/Odd grids"<<std::endl;
	setCheckerboard(rbShiftU,ShiftUe);
	setCheckerboard(rbShiftU,ShiftUo);
	//std::cout<<GridLogMessage << "\trbShiftU " <<norm2(rbShiftU)<<std::endl;

	//	std::cout<<GridLogMessage<<"Full grid shift"<<std::endl;
	ShiftU  = Cshift(U,dir,shift);    // Shift everything
	//	std::cout<<GridLogMessage << "\tShiftU " <<norm2(rbShiftU)<<std::endl;

	err = ShiftU - rbShiftU;
	std::cout<< "\terror " <<norm2(err)<<std::endl;

	Coordinate coor(4);

	std::cout<<GridLogMessage << "  Checking the non-checkerboard shift "<<shift <<" dir "<<dir <<" ... ";
	for(coor[3]=0;coor[3]<latt_size[3];coor[3]++){
	for(coor[2]=0;coor[2]<latt_size[2];coor[2]++){
	for(coor[1]=0;coor[1]<latt_size[1];coor[1]++){
	for(coor[0]=0;coor[0]<latt_size[0];coor[0]++){
	  
	  peekSite(cm,ShiftU,coor);

	  /////////	  double nrm=norm2(U);

	  Coordinate scoor(coor);
	  scoor[dir] = (scoor[dir]+shift)%latt_size[dir];
	  
#ifndef POWER10
	  Coordinate powers=latt_size;
	  Integer slex = scoor[0]
	    + latt_size[0]*scoor[1]
	    + latt_size[0]*latt_size[1]*scoor[2]
	    + latt_size[0]*latt_size[1]*latt_size[2]*scoor[3];
#else
	  Coordinate powers({1,10,100,1000});
	  Integer slex = scoor[0]
	    + 10        *scoor[1]
	    + 100       *scoor[2]
	    + 1000      *scoor[3];
#endif
	  Complex scm(slex);
	  
	  double nrm = abs(scm-cm()()());
	  Coordinate peer(4);
	  Complex ctmp = cm;
	  Integer index=real(ctmp);
	  Lexicographic::CoorFromIndex(peer,index,powers);

	  if (nrm > 0){
	    std::cout<<"FAIL shift "<< shift<<" in dir "<< dir
		     <<" ["<<coor[0]<<","<<coor[1]<<","<<coor[2]<<","<<coor[3]<<"] = "
		     << cm()()()<<" expect "<<scm<<"  "<<nrm<<std::endl;
	    std::cout<<"Got    "<<index<<" " << peer[0]<<","<<peer[1]<<","<<peer[2]<<","<<peer[3]<<std::endl;
	    index=real(scm);
	    Lexicographic::CoorFromIndex(peer,index,latt_size);
	    std::cout<<"Expect "<<index<<" " << peer[0]<<","<<peer[1]<<","<<peer[2]<<","<<peer[3]<<std::endl;
	    exit(-1);
	  }
	}}}}
	std::cout << " OK !"<<std::endl;

	int exx=0;
	std::cout<<GridLogMessage << "  Checking the     checkerboard shift "<< shift << " dir " << dir <<" ... ";
	for(coor[3]=0;coor[3]<latt_size[3];coor[3]++){
	for(coor[2]=0;coor[2]<latt_size[2];coor[2]++){
	for(coor[1]=0;coor[1]<latt_size[1];coor[1]++){
	for(coor[0]=0;coor[0]<latt_size[0];coor[0]++){
	  
	  peekSite(cm,rbShiftU,coor);

	  Integer checkerboard = RBFine.CheckerBoard(coor);

	  if ( checkerboard == ShiftUo.Checkerboard() ) {
	    peekSite(cmeo,ShiftUo,coor);
	  } else { 
	    peekSite(cmeo,ShiftUe,coor);
	  }

	  Coordinate scoor(coor);
	  scoor[dir] = (scoor[dir]+shift)%latt_size[dir];
	  
#ifndef POWER10
	  Coordinate powers=latt_size;
	  Integer slex = scoor[0]
	    + latt_size[0]*scoor[1]
	    + latt_size[0]*latt_size[1]*scoor[2]
	    + latt_size[0]*latt_size[1]*latt_size[2]*scoor[3];
#else 
	  Coordinate powers({1,10,100,1000});
	  Integer slex = scoor[0]
	    + 10        *scoor[1]
	    + 100       *scoor[2]
	    + 1000      *scoor[3];
#endif
	  Complex scm(slex);

	  Coordinate peer(4);
	  Complex ctmp=cmeo;
	  Integer index=real(ctmp);
	  Lexicographic::CoorFromIndex(peer,index,powers);

	  double nrm = abs(cmeo()()()-scm);
	  if (nrm != 0) {

	    std::cout << " coor "<<" ["<<coor[0]<<","<<coor[1]<<","<<coor[2]<<","<<coor[3]<<"] \n ";
	    std::cout << "shift "<< shift <<" dir "<<dir<< " checker board "<< checkerboard << " ";
	    std::cout << "Uo cb = "   << ShiftUo.Checkerboard() << " Ue cb= "<<ShiftUe.Checkerboard()<<std::endl;

	    std::cout<<"EOFAIL shift "<< shift<<" in dir "<< dir
		     <<" ["<<coor[0]<<","<<coor[1]<<","<<coor[2]<<","<<coor[3]<<"] = "
		     <<" cm " << cm()()()<<" cmeo "
		     << cmeo()()()<<" expect "<<scm<<"  "<<nrm<<std::endl;
	    std::cout<<"Got    "<<index<<" " << peer[0]<<","<<peer[1]<<","<<peer[2]<<","<<peer[3]<<std::endl;
	    index=real(scm);
	    Lexicographic::CoorFromIndex(peer,index,powers);
	    std::cout<<"Expect "<<index<<" " << peer[0]<<","<<peer[1]<<","<<peer[2]<<","<<peer[3]<<std::endl;
	    exx=1;

	  }

	  ctmp=cm;
	  index=real(ctmp);
	  nrm = abs(scm-cm()()());

	  if (nrm > 0){
	    std::cout<<"FAIL shift "<< shift<<" in dir "<< dir
		     <<" ["<<coor[0]<<","<<coor[1]<<","<<coor[2]<<","<<coor[3]<<"] = "
		     << cm()()()<<" expect "<<scm<<"  "<<nrm<<std::endl;
	    std::cout<<"Got    "<<index<<" " << peer[0]<<","<<peer[1]<<","<<peer[2]<<","<<peer[3]<<std::endl;
	    index=real(scm);
	    Lexicographic::CoorFromIndex(peer,index,powers);
	    std::cout<<"Expect "<<index<<" " << peer[0]<<","<<peer[1]<<","<<peer[2]<<","<<peer[3]<<std::endl;
	    exx=1;
	  } else if (0) { 
	    std::cout<<GridLogMessage<<"PASS shift "<< shift<<" in dir "<< dir
		     <<" ["<<coor[0]<<","<<coor[1]<<","<<coor[2]<<","<<coor[3]<<"] = "
		     << cm()()()<<" expect "<<scm<<"  "<<nrm<<std::endl;
	  }
	}}}}
	if (exx) exit(-1);
	std::cout << " OK !"<<std::endl;
    }
  }

  Grid_finalize();
}
