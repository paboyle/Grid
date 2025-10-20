    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/debug/Test_icosahedron.cc

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
#include <Grid/stencil/IcosahedralStencil.h>

using namespace std;
using namespace Grid;

const int MyNd=3;
template<typename vtype> using iIcosahedralLorentzComplex = iVector<iScalar<iScalar<vtype> >, MyNd+1 > ;

typedef iIcosahedralLorentzComplex<ComplexD >  IcosahedralLorentzComplexD;
typedef iIcosahedralLorentzComplex<vComplexD> vIcosahedralLorentzComplexD;
typedef Lattice<vIcosahedralLorentzComplexD> LatticeIcosahedralLorentzComplexD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int L=8;
  const int Npatch = IcosahedralPatches;

  // Put SIMD all in time direction
  Coordinate latt_size = GridDefaultLatt();
  Coordinate simd_layout({1,1,vComplexD::Nsimd(),1});
  Coordinate mpi_layout = GridDefaultMpi();

  std::cout << GridLogMessage << " mpi "<<mpi_layout<<std::endl;
  std::cout << GridLogMessage << " simd "<<simd_layout<<std::endl;
  std::cout << GridLogMessage << " latt "<<latt_size<<std::endl;
  GridCartesianCrossIcosahedron EdgeGrid(latt_size,simd_layout,mpi_layout,IcosahedralEdges);
  std::cout << GridLogMessage << " Created edge grid "<<std::endl;
  GridCartesianCrossIcosahedron VertexGrid(latt_size,simd_layout,mpi_layout,IcosahedralVertices);

  std::cout << GridLogMessage << " Created vertex grid "<<std::endl;
  LatticeIcosahedralLorentzComplexD Umu(&EdgeGrid);
  LatticeComplex Phi(&VertexGrid);
  std::cout << GridLogMessage << " Created two fields "<<std::endl;

  Phi = Zero();
  Umu = Zero();
  std::cout << GridLogMessage << " Zeroed two fields "<<std::endl;

  ComplexD one (1.0);
  Phi = one;
  Umu = one;
  
  std::cout << GridLogMessage << " V = "<<norm2(Phi)<<std::endl;
  std::cout << GridLogMessage << " Expect "<<latt_size[0]*latt_size[1]*latt_size[2]*10+2*latt_size[2]<<std::endl;
  
  std::cout << GridLogMessage << " E = "<<norm2(Umu)<<std::endl;
  std::cout << GridLogMessage << " Expect "<<latt_size[0]*latt_size[1]*latt_size[2]*10*4<<std::endl;

  //  std::cout << " Umu "<<Umu<<std::endl;
  //  std::cout << " Phi "<<Phi<<std::endl;
  LatticePole(Phi,South);
  //  std::cout << " Phi South Pole set\n"<<Phi<<std::endl;

  LatticePole(Phi,North);
  //  std::cout << " Phi North Pole set\n"<<Phi<<std::endl;

  for(int mu=0;mu<VertexGrid._ndimension;mu++){
    std::cout << " Calling lattice coordinate mu="<<mu<<std::endl;
    LatticeCoordinate(Phi,mu);
    //    std::cout << " Phi coor mu="<<mu<<"\n"<<Phi<<std::endl;
  }
  IcosahedralStencil stencil(&EdgeGrid);
  
  stencil.TestGeometry();

  std::cout << "Creating face stencil"<<std::endl;
  
  stencil.FaceStencil();
  Umu=one;
  LatticeComplex plaq1(&EdgeGrid);
  LatticeComplex plaq2(&EdgeGrid);

  {
    autoView(Umu_v,Umu,AcceleratorRead);
    autoView(plaq1_v,plaq1,AcceleratorWrite);
    autoView(plaq2_v,plaq2,AcceleratorWrite);
    autoView(stencil_v,stencil,AcceleratorRead);
    accelerator_for(ss,EdgeGrid.oSites(),vComplexD::Nsimd(),{

	const int x = IcosahedronPatchX;
	const int y = IcosahedronPatchY;
	const int d = IcosahedronPatchDiagonal;
	
	auto Lx = Umu_v(ss)(x);
	auto Ly = Umu_v(ss)(y);
	auto Ld = Umu_v(ss)(d);

	// for trace [ U_x(z) U_y(z+\hat x) adj(U_d(z)) ]
	{
	  auto SE1 = stencil_v.GetEntry(0,ss);
	  auto doAdj = SE1->_adjoint;
	  auto pol   = SE1->_polarisation;
	  auto s1    = SE1->_offset;
	  auto L1 = Umu_v(s1)(pol);
	  if(doAdj)
	    L1 = adj(L1);

	  coalescedWrite(plaq1_v[ss](),trace(Lx*L1*adj(Ld) ) );
	}
	
	// for trace [  U_y(z) adj(U_d(z))  U_x(z+\hat y) ]
	{
	  auto SE2 = stencil_v.GetEntry(1,ss);
	  auto doAdj = SE2->_adjoint;
	  auto pol   = SE2->_polarisation;
	  auto s2    = SE2->_offset;
	  auto L2 = Umu_v(s2)(pol);
	  if(doAdj)
	    L2 = adj(L2);

	  coalescedWrite(plaq2_v[ss](),trace(Ly*adj(Ld)*L2 ) );
	}
      });
  }
  std::cout << " plaq1 "<< norm2(plaq1)<<std::endl;
  std::cout << " plaq2 "<< norm2(plaq2)<<std::endl;
  Grid_finalize();
}
