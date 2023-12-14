    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_gfield_shift.cc

    Copyright (C) 2015

Author: Christopher Kelly <ckelly@bnl.gov>
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

//Test the shifting of the gauge field that respects the boundary conditions
#include <Grid/Grid.h>

using namespace Grid;
 ;

typedef ConjugateGimplR Gimpl; //can choose periodic / charge conjugate directions at wil
typedef Gimpl::GaugeField GaugeField;
typedef Gimpl::GaugeLinkField GaugeLinkField;
typedef Gimpl::SiteGaugeField SiteGaugeField;
typedef Gimpl::SiteGaugeLink SiteGaugeLink;

GaugeField CshiftGaugeField(const GaugeField &U, const int dir, const int shift){
  GridBase *Grid = U.Grid();

  GaugeField out(Grid);
  GaugeLinkField Umu(Grid);
  for(int mu=0;mu<Grid->Nd();mu++){
    Umu = PeekIndex<LorentzIndex>(U, mu);
    Umu = Gimpl::CshiftLink(Umu,dir,shift);
    PokeIndex<LorentzIndex>(out,Umu,mu);
  }
  return out;
}


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  auto latt_size   = GridDefaultLatt();
  auto simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();

  std::vector<int> conj_dirs = {1,1,0,0}; 
  Gimpl::setDirections(conj_dirs);

  GridCartesian        Fine(latt_size,simd_layout,mpi_layout);

  GridParallelRNG      FineRNG(&Fine);  FineRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));


  GaugeField U(&Fine);
  GaugeField ShiftU(&Fine);

  GaugeLinkField link_field(&Fine), link_field_2(&Fine);

  //Like Test_cshift we put the lex coordinate index on each site but make it imaginary
  //so we can tell when it was complex conjugated
  LatticeComplex lex(&Fine);
  lex=Zero();
  U = Zero();
  {
    LatticeComplex coor(&Fine);
    Integer stride =1;
    for(int d=0;d<4;d++){
      LatticeCoordinate(coor,d);
      lex = lex + coor*stride;
      stride=stride*latt_size[d];
    }
    PokeIndex<ColourIndex>(link_field, lex, 0,0); //place on 0,0 element of link

    for(int mu=0;mu<Nd;mu++){
      link_field_2 = link_field + mu*stride; //add in lex-mapping of mu
      link_field_2 = ComplexD(0,1) * link_field_2; //make imaginary
      PokeIndex<LorentzIndex>(U, link_field_2, mu);
    }
  }

  std::stringstream ss;
  ss<<"error";
  for(int d=0;d<Fine._ndimension;d++){
    ss<<"."<<Fine._processor_coor[d];
  }
  ss<<"_wr_"<<Fine._processor;
  std::string fname(ss.str());
  std::ofstream ferr(fname);
  
  Integer vol4d = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];

  bool fail = false;
  typename SiteGaugeField::scalar_object um;
  TComplex cm;
  
  for(int dir=0;dir<4;dir++){
    for(int shift=-latt_size[dir]+1;shift<latt_size[dir];shift++){
      if ( Fine.IsBoss() )
	std::cout<<GridLogMessage<<"Shifting by "<<shift<<" in direction "<<dir
		 << " dir is conj ? " << conj_dirs[dir] << std::endl;

      ShiftU = CshiftGaugeField(U,dir,shift);

      Coordinate coor(4);
      
      for(coor[3]=0;coor[3]<latt_size[3];coor[3]++){
      for(coor[2]=0;coor[2]<latt_size[2];coor[2]++){
      for(coor[1]=0;coor[1]<latt_size[1];coor[1]++){
      for(coor[0]=0;coor[0]<latt_size[0];coor[0]++){  
	peekSite(um,ShiftU,coor);

	Coordinate scoor(coor);
	scoor[dir] = (scoor[dir]+shift + latt_size[dir])%latt_size[dir];
	
	Integer slex = scoor[0]
	  + latt_size[0]*scoor[1]
	  + latt_size[0]*latt_size[1]*scoor[2]
	  + latt_size[0]*latt_size[1]*latt_size[2]*scoor[3];
	
	for(int mu = 0 ; mu < 4; mu++){
	  Integer slex_mu = slex + vol4d*mu;
	  Complex scm(0,slex_mu); //imaginary
	  if(
	     ( shift > 0 && coor[dir] >= latt_size[dir]-shift && conj_dirs[dir] ) 
	     ||
	     ( shift < 0 && coor[dir] <= -shift-1 && conj_dirs[dir] )
	     )
	    scm = conjugate(scm); //CC if pulled over boundary

	  cm = um(mu)()(0,0);

	  RealD nrm = abs(scm-cm()()());
	  //std::cout << cm << " " << scm << std::endl;

	  Coordinate peer(4);
	  Complex tmp  =cm;
	  Integer index=real(tmp);

	  Integer cm_mu = index / vol4d;
	  index = index % vol4d;
	  Lexicographic::CoorFromIndex(peer,index,latt_size);

	  if (nrm > 0){
	    ferr<<"FAIL mu " << mu << " shift "<< shift<<" in dir "<< dir<<" ["<<coor[0]<<","<<coor[1]<<","<<coor[2]<<","<<coor[3]<<"] = "<< cm()()()<<" expect "<<scm<<"  "<<nrm<<std::endl;
	    ferr<<"Got mu "<< cm_mu << " site " <<index<<" : " << peer[0]<<","<<peer[1]<<","<<peer[2]<<","<<peer[3]<<std::endl;

	    index=real(scm);
	    Integer scm_mu = index / vol4d;
	    index = index % vol4d;
	    Lexicographic::CoorFromIndex(peer,index,latt_size);
	    ferr<<"Expect mu " << scm_mu << " site " <<index<<": " << peer[0]<<","<<peer[1]<<","<<peer[2]<<","<<peer[3]<<std::endl;
	    fail = true;
	  }
	}

	}}}}
    }
  }

  if(fail) std::cout << "Test FAILED : see " << fname << " for more details" << std::endl;
  else std::cout << "Test Passed" << std::endl;

  Grid_finalize();
}
