    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_padded_cell.cc

    Copyright (C) 2015

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
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>

using namespace std;
using namespace Grid;

template<class vobj> void gpermute(vobj & inout,int perm){
  vobj tmp=inout;
  if (perm & 0x1 ) { permute(inout,tmp,0); tmp=inout;}
  if (perm & 0x2 ) { permute(inout,tmp,1); tmp=inout;}
  if (perm & 0x4 ) { permute(inout,tmp,2); tmp=inout;}
  if (perm & 0x8 ) { permute(inout,tmp,3); tmp=inout;}
}
  
int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Coordinate latt_size  = GridDefaultLatt();
  Coordinate simd_layout= GridDefaultSimd(Nd,vComplexD::Nsimd());
  Coordinate mpi_layout = GridDefaultMpi();
  std::cout << " mpi "<<mpi_layout<<std::endl;
  std::cout << " simd "<<simd_layout<<std::endl;
  std::cout << " latt "<<latt_size<<std::endl;
  GridCartesian GRID(latt_size,simd_layout,mpi_layout);

  GridParallelRNG   pRNG(&GRID);
  pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
  LatticeGaugeField Umu(&GRID);

  SU<Nc>::HotConfiguration(pRNG,Umu);

  Real plaq=WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu);
  LatticeComplex trplaq(&GRID);

  std::vector<LatticeColourMatrix> U(Nd, Umu.Grid());
  for (int mu = 0; mu < Nd; mu++) {
    U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
  }

  std::cout << GridLogMessage << " Average plaquette "<<plaq<<std::endl;

  LatticeComplex cplaq(&GRID); cplaq=Zero();

  /////////////////////////////////////////////////
  // Create a padded cell of extra padding depth=1
  /////////////////////////////////////////////////
  int depth = 1;
  PaddedCell Ghost(depth,&GRID);
  LatticeGaugeField Ughost = Ghost.Exchange(Umu);

  ///////////////////////////////////////////////////////////////////
  // Temporary debug Hack for single rank sim:
  // Check the contents of the cell are periodcally replicated
  // In future ONLY pad those dimensions that are not local to node
  ///////////////////////////////////////////////////////////////////
#if 0
  {
    double diff=0;
    double n=0;
  {
    autoView( Ug_v , Ughost, CpuRead);
    autoView( Ul_v , Umu   , CpuRead);
  for(int x=0;x<latt_size[0]+2;x++){
  for(int y=0;y<latt_size[1]+2;y++){
  for(int z=0;z<latt_size[2]+2;z++){
  for(int t=0;t<latt_size[3]+2;t++){
    int lx=(x-1+latt_size[0])%latt_size[0];
    int ly=(y-1+latt_size[1])%latt_size[1];
    int lz=(z-1+latt_size[2])%latt_size[2];
    int lt=(t-1+latt_size[3])%latt_size[3];
    Coordinate gcoor({x,y,z,t});
    Coordinate lcoor({lx,ly,lz,lt});
    LorentzColourMatrix g;
    LorentzColourMatrix l;
    peekLocalSite(g,Ug_v,gcoor);
    peekLocalSite(l,Ul_v,lcoor);
    g=g-l;
    assert(norm2(g)==0);
    diff = diff + norm2(g);
    n = n + norm2(l);
  }}}}
  }
  std::cout << "padded field check diff "<< diff <<" / "<< n<<std::endl;
  std::cout << norm2(Ughost)<< " " << norm2(Umu)<<std::endl;
  }
#endif

  ///// Array for the site plaquette
  GridBase *GhostGrid = Ughost.Grid();
  LatticeComplex gplaq(GhostGrid); 
  
  std::vector<Coordinate> shifts;
  for(int mu=0;mu<Nd;mu++){
    for(int nu=mu+1;nu<Nd;nu++){
  
      //    Umu(x) Unu(x+mu) Umu^dag(x+nu) Unu^dag(x)
      Coordinate shift_0(Nd,0);
      Coordinate shift_mu(Nd,0); shift_mu[mu]=1;
      Coordinate shift_nu(Nd,0); shift_nu[nu]=1;
      shifts.push_back(shift_0);
      shifts.push_back(shift_mu);
      shifts.push_back(shift_nu);
      shifts.push_back(shift_0);
    }
  }
  GeneralLocalStencil gStencil(GhostGrid,shifts);

  gplaq=Zero();
  {
    autoView( gp_v , gplaq, CpuWrite);
    autoView( t_v , trplaq, CpuRead);
    autoView( U_v , Ughost, CpuRead);
    for(int ss=0;ss<gp_v.size();ss++){
      int s=0;
      for(int mu=0;mu<Nd;mu++){
	for(int nu=mu+1;nu<Nd;nu++){

	  auto SE0 = gStencil.GetEntry(s+0,ss);
	  auto SE1 = gStencil.GetEntry(s+1,ss);
	  auto SE2 = gStencil.GetEntry(s+2,ss);
	  auto SE3 = gStencil.GetEntry(s+3,ss);
	
	  int o0 = SE0->_offset;
	  int o1 = SE1->_offset;
	  int o2 = SE2->_offset;
	  int o3 = SE3->_offset;
	  
	  auto U0 = U_v[o0](mu);
	  auto U1 = U_v[o1](nu);
	  auto U2 = adj(U_v[o2](mu));
	  auto U3 = adj(U_v[o3](nu));

	  gpermute(U0,SE0->_permute);
	  gpermute(U1,SE1->_permute);
	  gpermute(U2,SE2->_permute);
	  gpermute(U3,SE3->_permute);
	  
	  gp_v[ss]() =gp_v[ss]() + trace( U0*U1*U2*U3 );
	  s=s+4;
	}
      }
    }
  }
  cplaq = Ghost.Extract(gplaq);
  RealD vol = cplaq.Grid()->gSites();
  RealD faces = (Nd * (Nd-1))/2;
  auto p = TensorRemove(sum(cplaq));
  auto result = p.real()/vol/faces/Nc;

  std::cout << GridLogMessage << " Average plaquette via padded cell "<<result<<std::endl;
  std::cout << GridLogMessage << " Diff "<<result-plaq<<std::endl;
  
  assert(fabs(result-plaq)<1.0e-8);
  Grid_finalize();
}
