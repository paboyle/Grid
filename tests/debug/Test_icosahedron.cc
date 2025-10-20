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

template<typename vtype> using iIcosahedralLorentzComplex      = iVector<iScalar<iScalar<vtype> >, MyNd+1 > ;
template<typename vtype> using iIcosahedralLorentzColourMatrix = iVector<iScalar<iMatrix<vtype,Nc> >, MyNd+1 > ;
template<typename vtype> using iIcosahedralColourMatrix        = iScalar<iScalar<iMatrix<vtype,Nc> > > ;

typedef iIcosahedralLorentzComplex<Complex >  IcosahedralLorentzComplex;
typedef iIcosahedralLorentzComplex<vComplex> vIcosahedralLorentzComplex;
typedef Lattice<vIcosahedralLorentzComplex> LatticeIcosahedralLorentzComplex;

typedef iIcosahedralLorentzColourMatrix<Complex >  IcosahedralLorentzColourMatrix;
typedef iIcosahedralLorentzColourMatrix<vComplex> vIcosahedralLorentzColourMatrix;
typedef Lattice<vIcosahedralLorentzColourMatrix>  LatticeIcosahedralLorentzColourMatrix;

typedef iIcosahedralColourMatrix<Complex >  IcosahedralColourMatrix;
typedef iIcosahedralColourMatrix<vComplex> vIcosahedralColourMatrix;
typedef Lattice<vIcosahedralColourMatrix>  LatticeIcosahedralColourMatrix;


class IcosahedralGimpl
{
public:
  typedef LatticeIcosahedralLorentzColourMatrix  GaugeField;
  typedef LatticeIcosahedralColourMatrix         GaugeLinkField;
  typedef LatticeComplex                          ComplexField;
};
  
template< class Gimpl>
class IcosahedralEdgeSupport {
public:
  // Move to inherit types macro as before
  typedef typename Gimpl::GaugeField GaugeField;
  typedef typename Gimpl::GaugeLinkField GaugeLinkField;
  typedef typename Gimpl::ComplexField ComplexField;
  //
  GridBase *VertexGrid;
  GridBase *EdgeGrid;
  IcosahedralStencil FaceStencil;
  
  IcosahedralEdgeSupport(GridBase *_VertexGrid,GridBase *_EdgeGrid)
    : FaceStencil (EdgeGrid), VertexGrid(_VertexGrid), EdgeGrid(_EdgeGrid)
  {
    FaceStencil.FaceStencil();
  }

  ////////////////////////////////////////////////////////////////////////////////////
  // Fixme: will need a version of "Gimpl" and a wrapper class following "WilsonLoops" style.
  // Gauge Link field GT is the gauge transform and lives on the VERTEX field
  ////////////////////////////////////////////////////////////////////////////////////
  void ForwardTriangles(GaugeField &Umu,LatticeComplex &plaq1,LatticeComplex &plaq2)
  {
    {
      autoView(Umu_v,Umu,AcceleratorRead);
      autoView(plaq1_v,plaq1,AcceleratorWrite);
      autoView(plaq2_v,plaq2,AcceleratorWrite);
      autoView(stencil_v,FaceStencil,AcceleratorRead);

      accelerator_for(ss,EdgeGrid->oSites(),vComplex::Nsimd(),{

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
  }

  void  GaugeTransform(GaugeLinkField &gt, GaugeField &Umu)
  {
    assert(gt.Grid()==VertexGrid);
    assert(Umu.Grid()==EdgeGrid);
    
    GridBase * vgrid = VertexGrid;
    GridBase * grid  = EdgeGrid;

    int osites  = grid->oSites();

    uint64_t cart_sites = grid->CartesianOsites();
    uint64_t Npole_sites = grid->NorthPoleOsites();
    uint64_t Spole_sites = grid->SouthPoleOsites();
    Coordinate pcoor = grid->ThisProcessorCoor();
    Coordinate pgrid = grid->ProcessorGrid();
    /*
     * resize the stencil entries array and set npoints
     */
    autoView(g_v,gt,CpuRead);
    autoView(Umu_v,Umu,CpuWrite);
    for(uint64_t site=0;site<cart_sites; site ++) {

      Coordinate Coor;
      Coordinate NbrCoor;

      int nd = grid->Nd();
      int L  = grid->LocalDimensions()[0];
	
      ////////////////////////////////////////////////
      // Outer index of neighbour Offset calculation
      ////////////////////////////////////////////////
      grid->oCoorFromOindex(Coor,site);
      NbrCoor = Coor;
      assert( grid->LocalDimensions()[1]==grid->LocalDimensions()[0]);
      assert( grid->_simd_layout[0]==1); // Cannot vectorise in these dims
      assert( grid->_simd_layout[1]==1);
      assert( grid->_processors[0]==1); // Cannot mpi distribute in these dims
      assert( grid->_processors[1]==1);

      int Patch = Coor[nd-1];
      int HemiPatch = Patch%HemiPatches;
      int north = Patch/HemiPatches;
      int south = 1-north;
      int isPoleY;
      int isPoleX;
      
      assert(Patch<IcosahedralPatches);
      assert((north==1)||(south==1));

      Coordinate XpCoor;
      Coordinate YpCoor;
      Coordinate DpCoor;

      FaceStencil.GetNbrForPlusDiagonal(grid,Coor,DpCoor);
      FaceStencil.GetNbrForPlusX(grid,Coor,XpCoor,isPoleX);
      FaceStencil.GetNbrForPlusY(grid,Coor,YpCoor,isPoleY);

      int XpHemiPatch  = XpCoor[nd-1]%HemiPatches;
      int XpHemisphere = XpCoor[nd-1]/HemiPatches;

      int DpPatch  = DpCoor[nd-1];
      int DpHemiPatch  = DpCoor[nd-1]%HemiPatches;
      int DpHemisphere = DpCoor[nd-1]/HemiPatches;

      // Work out the pole_osite
      Coordinate rdims;
      Coordinate ocoor;
      int64_t pole_osite;
      int Ndm1 = grid->Nd()-1;
      for(int d=2;d<Ndm1;d++){
	int dd=d-2;
	rdims.push_back(grid->_rdimensions[d]);
	ocoor.push_back(Coor[d]%grid->_rdimensions[d]);
      }
      Lexicographic::IndexFromCoor(ocoor,pole_osite,rdims);
      
      uint64_t xp_idx;
      uint64_t yp_idx;
      uint64_t dp_idx;
      if ( isPoleX ) {
	assert(vgrid->ownsSouthPole());
	xp_idx = pole_osite + grid->SouthPoleOsite();
      } else {
	xp_idx = grid->oIndex(XpCoor);
      }
      if ( isPoleY ) {
	assert(vgrid->ownsNorthPole());
	yp_idx = pole_osite + grid->NorthPoleOsite();
      } else {
	yp_idx = grid->oIndex(YpCoor);
      }
      dp_idx = grid->oIndex(DpCoor);

      auto g  = g_v(site)();
      auto gx = g_v(xp_idx)();
      auto gy = g_v(yp_idx)();
      auto gd = g_v(dp_idx)();

      auto lx = Umu_v(site)(IcosahedronPatchX);
      auto ly = Umu_v(site)(IcosahedronPatchY);
      auto ld = Umu_v(site)(IcosahedronPatchDiagonal);
      
      lx = g*lx*adj(gx);
      ly = g*ly*adj(gy);
      ld = g*ld*adj(gd);

      coalescedWrite(Umu_v[site](IcosahedronPatchX),lx);
      coalescedWrite(Umu_v[site](IcosahedronPatchY),ly);
      coalescedWrite(Umu_v[site](IcosahedronPatchDiagonal),ld);
    };
  }

};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int L=8;
  const int Npatch = IcosahedralPatches;

  // Put SIMD all in time direction
  Coordinate latt_size = GridDefaultLatt();
  Coordinate simd_layout({1,1,vComplex::Nsimd(),1});
  Coordinate mpi_layout = GridDefaultMpi();

  std::cout << GridLogMessage << " mpi "<<mpi_layout<<std::endl;
  std::cout << GridLogMessage << " simd "<<simd_layout<<std::endl;
  std::cout << GridLogMessage << " latt "<<latt_size<<std::endl;
  GridCartesianCrossIcosahedron EdgeGrid(latt_size,simd_layout,mpi_layout,IcosahedralEdges);
  std::cout << GridLogMessage << " Created edge grid "<<std::endl;
  GridCartesianCrossIcosahedron VertexGrid(latt_size,simd_layout,mpi_layout,IcosahedralVertices);

  std::cout << GridLogMessage << " Created vertex grid "<<std::endl;
  LatticeIcosahedralLorentzColourMatrix Umu(&EdgeGrid);
  LatticeComplex Phi(&VertexGrid);
  std::cout << GridLogMessage << " Created two fields "<<std::endl;

  Phi = Zero();
  Umu = Zero();
  std::cout << GridLogMessage << " Zeroed two fields "<<std::endl;

  Complex one (1.0);
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

  std::cout << "Creating face stencil"<<std::endl;
  IcosahedralEdgeSupport<IcosahedralGimpl> Support(&VertexGrid,&EdgeGrid);

  std::cout << " Calling Test Geometry "<<std::endl;
  Support.FaceStencil.TestGeometry();


  Umu=one;
  LatticeComplex plaq1(&EdgeGrid);
  LatticeComplex plaq2(&EdgeGrid);

  Support.ForwardTriangles(Umu,plaq1,plaq2);
  std::cout << " plaq1 "<< norm2(plaq1)<<std::endl;
  std::cout << " plaq2 "<< norm2(plaq2)<<std::endl;

  // Random gauge xform
  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG  vRNG(&EdgeGrid);   vRNG.SeedFixedIntegers(seeds);
  LatticeIcosahedralColourMatrix g(&VertexGrid);
  
  //  SU<Nc>::LieRandomize(vRNG,g);

  LatticeReal gr(&VertexGrid);
  LatticeComplex gc(&VertexGrid);
  gaussian(vRNG,gr);
  Complex ci(0.0,1.0);
  gc = toComplex(gr);
  g=one;
  g = g * exp(ci*gc);
  
  std::cout << "applying gauge transform"<<std::endl;
  Support.GaugeTransform(g,Umu);
  std::cout << "applied gauge transform "<<Umu<<std::endl;
  
  Support.ForwardTriangles(Umu,plaq1,plaq2);
  std::cout << " plaq1 "<< norm2(plaq1)<<std::endl;
  std::cout << " plaq2 "<< norm2(plaq2)<<std::endl;

  
  Grid_finalize();
}
