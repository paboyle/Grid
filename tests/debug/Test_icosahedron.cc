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
class IcosahedralSupport {
public:
  // Move to inherit types macro as before
  typedef typename Gimpl::GaugeField GaugeField;
  typedef typename Gimpl::GaugeLinkField GaugeLinkField;
  typedef typename Gimpl::ComplexField ComplexField;
  //
  GridBase *VertexGrid;
  GridBase *EdgeGrid;

  IcosahedralStencil FaceStencil;
  IcosahedralStencil NNee; // edge neighbours with edge domain
  IcosahedralStencil NNev; // vertex neighbours but in edge domain
  IcosahedralStencil NNvv; // vertex neighbours with vertex domain
  
  IcosahedralSupport(GridBase *_VertexGrid,GridBase *_EdgeGrid)
    : FaceStencil (_EdgeGrid), NNee(_EdgeGrid), NNev(_VertexGrid), NNvv(_VertexGrid), VertexGrid(_VertexGrid), EdgeGrid(_EdgeGrid)
  {
    FaceStencil.FaceStencil();
    NNee.NearestNeighbourStencil(false);// Edge nearest neighbour
    NNev.NearestNeighbourStencil(false);// Edge result, vertex neighbour
    NNvv.NearestNeighbourStencil(true); // vertex result and neighbour
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

	// This was wrong after GT
	// Could be EITHER the GT or the the plaq / stencil
	// for trace [  U_y(z) U_x(z+\hat y) adj(U_d(z))   ]
	{
	  auto SE2 = stencil_v.GetEntry(1,ss);
	  auto doAdj = SE2->_adjoint;
	  auto pol   = SE2->_polarisation;
	  auto s2    = SE2->_offset;
	  auto L2 = Umu_v(s2)(pol);
	  if(doAdj)
	    L2 = adj(L2);
	  coalescedWrite(plaq2_v[ss](),trace(Ly*L2*adj(Ld) ) );
	  //	  std::cout << "site "<< ss<<" plaq "<< plaq2_v[ss] << " doAdj "<< (int) doAdj<<" pol "<<(int) pol <<std::endl;
	}
      });
    }
  }

  // Staples for gauge force
  void IcosahedralStaples(GaugeField &Umu,
			  GaugeLinkField &stapleXY,
			  GaugeLinkField &stapleYX,
			  GaugeLinkField &stapleXD,
			  GaugeLinkField &stapleDX,
			  GaugeLinkField &stapleYD,
			  GaugeLinkField &stapleDY)
  {
      autoView(Umu_v,Umu,AcceleratorRead);
      autoView(stapleXY_v,stapleXY,AcceleratorWrite);
      autoView(stapleYX_v,stapleYX,AcceleratorWrite);
      autoView(stapleXD_v,stapleXD,AcceleratorWrite);
      autoView(stapleDX_v,stapleDX,AcceleratorWrite);
      autoView(stapleYD_v,stapleYD,AcceleratorWrite);
      autoView(stapleDY_v,stapleDY,AcceleratorWrite);
      autoView(stencil_v,NNee,AcceleratorRead);

      const int np = NNee._npoints;
      
      accelerator_for(ss,EdgeGrid->oSites(),vComplex::Nsimd(),{

	const int ent_Xp = 0;
	const int ent_Yp = 1;
	const int ent_Xm = 3;
	const int ent_Ym = 4;

	Integer lexXp = ss*np+ent_Xp;
	Integer lexYp = ss*np+ent_Yp;
	Integer lexXm = ss*np+ent_Xm;
	Integer lexYm = ss*np+ent_Ym;
	//	Integer lexDp = ss*np+2; // Not touched by staples.
	//	Integer lexDm = ss*np+5;
	  
	const int x = IcosahedronPatchX;
	const int y = IcosahedronPatchY;
	const int d = IcosahedronPatchDiagonal;

	// Three forward links from this site
	auto Lx = Umu_v(ss)(x);
	auto Ly = Umu_v(ss)(y);
	auto Ld = Umu_v(ss)(d);

	///////////////////////////////////////////////////////////////////
	// Terms for the staple orthog to PlusDiagonal
	///////////////////////////////////////////////////////////////////
	// adj( U_y(z+\hat x)) adj(U_x(z)) 
	{
	  auto SE1 = stencil_v.GetEntry(ent_Xp,ss);
	  auto doAdj = SE1->_adjoint;
	  auto pol   = SE1->_polarisation;
	  auto s1    = SE1->_offset;
	  auto Ly_at_xp = Umu_v(s1)(pol);
	  if(doAdj)
	    Ly_at_xp = adj(Ly_at_xp);
	  coalescedWrite(stapleXY_v[ss](),adj(Ly_at_xp)*adj(Lx)  );
	}
	// adj( U_y(z) )  adj(U_x(z+\hat y))
	{
	  auto SE2 = stencil_v.GetEntry(ent_Yp,ss);
	  auto doAdj = SE2->_adjoint;
	  auto pol   = SE2->_polarisation;
	  auto s2    = SE2->_offset;
	  auto Lx_at_yp = Umu_v(s2)(pol);
	  if(doAdj)
	    Lx_at_yp = adj(Lx_at_yp);

	  coalescedWrite(stapleYX_v[ss](),adj(Lx_at_yp)*adj(Ly) );
	}
	///////////////////////////////////////////////////////////////////
	// Terms for the staple covering Xp : Dp Yp(x++) and Yp(y--) Dp(y--)
	///////////////////////////////////////////////////////////////////
	// U_y(z+\hat x)*adj(U_d(z)) 
	{
	  auto SE1 = stencil_v.GetEntry(ent_Xp,ss);
	  auto doAdj = SE1->_adjoint;
	  auto pol   = SE1->_polarisation;
	  auto s1    = SE1->_offset;
	  auto Ly_at_xp = Umu_v(s1)(pol);
	  if(doAdj)
	    Ly_at_xp = adj(Ly_at_xp);
	  coalescedWrite(stapleDY_v[ss](), Ly_at_xp *adj(Ld));
	}
	//  adj(U_d(z-\hat y)) U_y(z-\hat y)
	{
	  auto SE2 = stencil_v.GetEntry(ent_Ym,ss);
	  auto doAdj = SE2->_adjoint;
	  auto pol   = SE2->_polarisation;
	  auto s2    = SE2->_offset;
	  int pol1 = IcosahedronPatchY;
	  int pol2 = IcosahedronPatchDiagonal;
	  if ( pol != IcosahedronPatchDiagonal ) {
	    pol1 = IcosahedronPatchDiagonal;
	    pol2 = IcosahedronPatchX;
	  }
	  auto Ly_at_ym = Umu_v(s2)(pol1);
	  auto Ld_at_ym = Umu_v(s2)(pol2);
	  coalescedWrite(stapleYD_v[ss](),adj(Ld_at_ym)*Ly_at_ym );
	}
	///////////////////////////////////////////////////////////////////
	// Terms for the staple covering Yp : Dp Xp(y++) and Xp(x--) Dp(x--)
	///////////////////////////////////////////////////////////////////
	//  U_x(z+\hat y)adj(U_d(z))
	{
	  auto SE1 = stencil_v.GetEntry(ent_Yp,ss);
	  auto doAdj = SE1->_adjoint;
	  auto pol   = SE1->_polarisation;
	  auto s1    = SE1->_offset;
	  auto Lx_at_yp = Umu_v(s1)(pol);
	  if(doAdj)
	    Lx_at_yp = adj(Lx_at_yp);
	  coalescedWrite(stapleDX_v[ss](),Lx_at_yp*adj(Ld) );
	}

	//  adj(U_d(z-\hat x))U_x(z-\hat x)
	{
	  auto SE2 = stencil_v.GetEntry(ent_Xm,ss);
	  auto doAdj = SE2->_adjoint;
	  auto pol   = SE2->_polarisation;
	  auto s2    = SE2->_offset;
	  int pol1 = IcosahedronPatchX;
	  int pol2 = IcosahedronPatchDiagonal;
	  if ( pol != IcosahedronPatchDiagonal ) {
	    pol1 = IcosahedronPatchDiagonal;
	    pol2 = IcosahedronPatchY;
	  }
	  auto Ly = Umu_v(ss)(IcosahedronPatchY);
	  auto Lx_at_xm = Umu_v(s2)(pol1);
	  auto Ld_at_xm = Umu_v(s2)(pol2);
	  coalescedWrite(stapleXD_v[ss](),adj(Ld_at_xm)*Lx_at_xm );
	}
      });
  }
  
  /*
   * Should be able to use the Vertex based stencil to do the GT, picking forward hops
   */
  template<class MatterField>
  void Laplacian(MatterField &in,MatterField &out)
  {
    autoView(out_v,out,AcceleratorWrite);
    autoView(in_v,in,AcceleratorRead);
    autoView(stencil_v,NNvv,AcceleratorRead);

    const int np = NNvv._npoints;

    std::cout << GridLogMessage<< "Free Laplacian via STENCIL "<<std::endl;

    accelerator_for(ss,VertexGrid->oSites(),vComplex::Nsimd(),{

      const int ent_Xp = 0;
      const int ent_Yp = 1;
      const int ent_Dp = 2;
      const int ent_Xm = 3;
      const int ent_Ym = 4;
      const int ent_Dm = 5;
      
      Integer lexXp = ss*np+ent_Xp;
      Integer lexYp = ss*np+ent_Yp;
      Integer lexDp = ss*np+ent_Dp;
      Integer lexXm = ss*np+ent_Xm;
      Integer lexYm = ss*np+ent_Ym;
      Integer lexDm = ss*np+ent_Dm;
	  
      const int x = IcosahedronPatchX;
      const int y = IcosahedronPatchY;
      const int d = IcosahedronPatchDiagonal;

      //      // Three forward links from this site
      //      auto Lx = Umu_v(ss)(x);
      //      auto Ly = Umu_v(ss)(y);
      //      auto Ld = Umu_v(ss)(d);

      auto SE = stencil_v.GetEntry(ent_Xp,ss);
      uint64_t xp_idx = SE->_offset;

      SE = stencil_v.GetEntry(ent_Yp,ss);
      uint64_t yp_idx = SE->_offset;

      SE = stencil_v.GetEntry(ent_Dp,ss);
      uint64_t dp_idx = SE->_offset;

      SE = stencil_v.GetEntry(ent_Xm,ss);
      uint64_t xm_idx = SE->_offset;

      SE = stencil_v.GetEntry(ent_Ym,ss);
      uint64_t ym_idx = SE->_offset;

      SE = stencil_v.GetEntry(ent_Dm,ss);
      uint64_t dm_idx = SE->_offset;
      int missingLink = SE->_missing_link;
      
      auto i    = in_v(ss)();
      auto inxp = in_v(xp_idx)();
      auto inyp = in_v(yp_idx)();
      auto indp = in_v(dp_idx)();
      auto inxm = in_v(xm_idx)();
      auto inym = in_v(ym_idx)();
      auto indm = in_v(dm_idx)();

      auto o    = out_v(ss)();

      if ( missingLink ) {
	o = (1.0/5.0)*(inxp+inyp+indp+inxm+inym)-o;
      } else {
	o = (1.0/6.0)*(inxp+inyp+indp+inxm+inym+indm)-o;
      }
      
      coalescedWrite(out_v[ss](),o);
      });
  }
    
  
  void  GaugeTransform(GaugeLinkField &gt, GaugeField &Umu)
  {
    autoView(Umu_v,Umu,AcceleratorWrite);
    autoView(g_v,gt,AcceleratorRead);
    autoView(stencil_v,NNev,AcceleratorRead);

    const int np = NNev._npoints;

    std::cout << GridLogMessage<< "GaugeTransform via STENCIL "<<std::endl;

    accelerator_for(ss,EdgeGrid->oSites(),vComplex::Nsimd(),{

      const int ent_Xp = 0;
      const int ent_Yp = 1;
      const int ent_Dp = 2;
      
      Integer lexXp = ss*np+ent_Xp;
      Integer lexYp = ss*np+ent_Yp;
      Integer lexDp = ss*np+ent_Dp;
	  
      // Three forward links from this site
      auto Lx = Umu_v(ss)(IcosahedronPatchX);
      auto Ly = Umu_v(ss)(IcosahedronPatchY);
      auto Ld = Umu_v(ss)(IcosahedronPatchDiagonal);

      auto SE = stencil_v.GetEntry(ent_Xp,ss);
      uint64_t xp_idx = SE->_offset;

      SE = stencil_v.GetEntry(ent_Yp,ss);
      uint64_t yp_idx = SE->_offset;

      SE = stencil_v.GetEntry(ent_Dp,ss);
      uint64_t dp_idx = SE->_offset;

      auto g  = g_v(ss)();
      auto gx = g_v(xp_idx)();
      auto gy = g_v(yp_idx)();
      auto gd = g_v(dp_idx)();

      auto lx = Umu_v(ss)(IcosahedronPatchX);
      auto ly = Umu_v(ss)(IcosahedronPatchY);
      auto ld = Umu_v(ss)(IcosahedronPatchDiagonal);

      lx = g*lx*adj(gx);
      ly = g*ly*adj(gy);
      ld = g*ld*adj(gd);
      
      coalescedWrite(Umu_v[ss](IcosahedronPatchX),lx);
      coalescedWrite(Umu_v[ss](IcosahedronPatchY),ly);
      coalescedWrite(Umu_v[ss](IcosahedronPatchDiagonal),ld);
      });
  }
  /*
   * This routine is slow and single threaded on CPU
   */
  void  GaugeTransformCPU(GaugeLinkField &gt, GaugeField &Umu)
  {
    assert(gt.Grid()==VertexGrid);
    assert(Umu.Grid()==EdgeGrid);
    assert(VertexGrid->isIcosahedralVertex());
    assert(EdgeGrid->isIcosahedralEdge());
    
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
	xp_idx = pole_osite + vgrid->SouthPoleOsite();
      } else {
	xp_idx = grid->oIndex(XpCoor);
      }
      if ( isPoleY ) {
	assert(vgrid->ownsNorthPole());
	yp_idx = pole_osite + vgrid->NorthPoleOsite();
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
      
      /*
	std::cout << "site "<<site<<std::endl;
	std::cout << " xp_idx "<<xp_idx<<std::endl;
	std::cout << " yp_idx "<<yp_idx<<std::endl;
	std::cout << " dp_idx "<<dp_idx<<std::endl;
	std::cout << " g  "<<g<<std::endl;
	std::cout << " gx "<<gx<<std::endl;
	std::cout << " gy "<<gy<<std::endl;
	std::cout << " gd "<<gd<<std::endl;
	std::cout << " lx "<<lx<<std::endl;
	std::cout << " ly "<<ly<<std::endl;
	std::cout << " ld "<<ld<<std::endl;
      */
      
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
  LatticeIcosahedralLorentzColourMatrix Umuck(&EdgeGrid);
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

  //  std::cout << " Umu "<<Umu<<std::endl; // debugged, so comment out to reduce verbose
  //  std::cout << " Phi "<<Phi<<std::endl; // debugged, so comment out to reduce verbose
  LatticePole(Phi,South);
  //  std::cout << " Phi South Pole set\n"<<Phi<<std::endl; // debugged, so comment out to reduce verbose

  LatticePole(Phi,North);
  //  std::cout << " Phi North Pole set\n"<<Phi<<std::endl; // debugged, so comment out to reduce verbose

  for(int mu=0;mu<VertexGrid._ndimension;mu++){
    std::cout << GridLogMessage << " Calling lattice coordinate mu="<<mu<<std::endl;
    LatticeCoordinate(Phi,mu);
    //    std::cout << GridLogMessage << " Phi coor mu="<<mu<<"\n"<<Phi<<std::endl; // debugged, so comment out to reduce verbose
  }

  std::cout << GridLogMessage << "Creating face stencil"<<std::endl;
  IcosahedralSupport<IcosahedralGimpl> Support(&VertexGrid,&EdgeGrid);

  std::cout << GridLogMessage << " Calling Test Geometry "<<std::endl;
  Support.FaceStencil.TestGeometry();


  Umu=one;
  LatticeComplex plaq1(&EdgeGrid);
  LatticeComplex plaq2(&EdgeGrid);

  LatticeComplex plaq_ref(&EdgeGrid);
  plaq_ref=1.0;

  Support.ForwardTriangles(Umu,plaq1,plaq2);
  std::cout << GridLogMessage << " plaq1 "<< norm2(plaq1)<<std::endl;
  std::cout << GridLogMessage << " plaq2 "<< norm2(plaq2)<<std::endl;

  std::cout << GridLogMessage << " plaq1 err "<< norm2(plaq1-plaq_ref)<<std::endl;
  std::cout << GridLogMessage << " plaq2 err "<< norm2(plaq2-plaq_ref)<<std::endl;

  // Random gauge xform
  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG  vRNG(&EdgeGrid);   vRNG.SeedFixedIntegers(seeds);
  
  //  SU<Nc>::LieRandomize(vRNG,g);
  LatticeIcosahedralColourMatrix g(&VertexGrid);
  LatticeReal gr(&VertexGrid);
  LatticeComplex gc(&VertexGrid);
  gr = 1.0;
  gaussian(vRNG,gr);
  Complex ci(0.0,1.0);
  gc = toComplex(gr);
  g=one;
  g = g * exp(ci*gc);
  
  std::cout << GridLogMessage << "****************************************"<<std::endl;
  std::cout << GridLogMessage << " Check plaquette is gauge invariant "<<std::endl;
  std::cout << GridLogMessage << "****************************************"<<std::endl;
  std::cout << GridLogMessage << " applying gauge transform"<<std::endl;
  Umuck = Umu;
  Support.GaugeTransform    (g,Umuck);
  Support.GaugeTransformCPU(g,Umu);
  Umuck = Umuck - Umu;
  std::cout << GridLogMessage <<"Diff between reference GT and stencil GT: "<<norm2(Umuck) <<std::endl;
  std::cout << GridLogMessage << " applied gauge transform "<<std::endl;
  //  std::cout << "Umu\n"<< Umu << std::endl;
  std::cout << GridLogMessage << " recalculating plaquette "<<std::endl;
  Support.ForwardTriangles(Umu,plaq1,plaq2);
  std::cout << GridLogMessage << " plaq1 "<< norm2(plaq1)<<std::endl;
  std::cout << GridLogMessage << " plaq2 "<< norm2(plaq2)<<std::endl;

  //  std::cout << " plaq1 "<< plaq1<<std::endl;
  //  std::cout << " plaq2 "<< plaq2<<std::endl;

  std::cout << GridLogMessage << " plaq1 err "<< norm2(plaq1-plaq_ref)<<std::endl;
  std::cout << GridLogMessage << " plaq2 err "<< norm2(plaq2-plaq_ref)<<std::endl;

  typedef IcosahedralGimpl::GaugeLinkField GaugeLinkField;
  typedef IcosahedralGimpl::GaugeField     GaugeField;

  GaugeLinkField stapleXY(&EdgeGrid);
  GaugeLinkField stapleYX(&EdgeGrid);
  GaugeLinkField stapleXD(&EdgeGrid);
  GaugeLinkField stapleDX(&EdgeGrid);
  GaugeLinkField stapleYD(&EdgeGrid);
  GaugeLinkField stapleDY(&EdgeGrid);

  GaugeLinkField linkX(&EdgeGrid);
  GaugeLinkField linkY(&EdgeGrid);
  GaugeLinkField linkD(&EdgeGrid);
  std::cout << GridLogMessage << "****************************************"<<std::endl;
  std::cout << GridLogMessage << " Check triangular staples match plaquette "<<std::endl;
  std::cout << GridLogMessage << "****************************************"<<std::endl;

  Support.IcosahedralStaples(Umu,stapleXY,stapleYX,
			     stapleXD,stapleDX,
			     stapleYD,stapleDY);

  linkX = peekLorentz(Umu,IcosahedronPatchX);
  linkY = peekLorentz(Umu,IcosahedronPatchY);
  linkD = peekLorentz(Umu,IcosahedronPatchDiagonal);

  // OK
  std::cout << GridLogMessage << " trace D*StapleXY "<<norm2(trace(linkD * stapleXY))<<std::endl;
  std::cout << GridLogMessage << " err " << norm2(trace(linkD * stapleXY)-plaq_ref)<<std::endl;

  // BAD
  std::cout << GridLogMessage << " trace D*StapleYX "<<norm2(trace(linkD * stapleYX))<<std::endl;
  std::cout << GridLogMessage << " err " << norm2(trace(linkD * stapleYX)-plaq_ref)<<std::endl;

  std::cout << GridLogMessage << " trace X*StapleYD "<<norm2(trace(linkX * stapleYD))<<std::endl;
  std::cout << GridLogMessage << " err " << norm2(trace(linkX * stapleYD)-plaq_ref)<<std::endl;
  std::cout << GridLogMessage << " trace X*StapleDY "<<norm2(trace(linkX * stapleDY))<<std::endl;
  std::cout << GridLogMessage << " err " << norm2(trace(linkX * stapleDY)-plaq_ref)<<std::endl;

  std::cout << GridLogMessage << " trace Y*StapleXD "<<norm2(trace(linkY * stapleXD))<<std::endl;
  std::cout << GridLogMessage << " err " << norm2(trace(linkY * stapleXD)-plaq_ref)<<std::endl;
  std::cout << GridLogMessage << " trace Y*StapleDX "<<norm2(trace(linkY * stapleDX))<<std::endl;
  std::cout << GridLogMessage << " err " << norm2(trace(linkY * stapleDX)-plaq_ref)<<std::endl;

  //  std::cout << " D " << linkD<<std::endl;
  //  std::cout << " X " << linkX<<std::endl;
  //  std::cout << " Y " << linkY<<std::endl;
  //  std::cout << " DXY\n " << closure(linkD * stapleYX) <<std::endl;
  //  std::cout << " YXD\n " << closure(linkY * stapleXD) <<std::endl;

  
  std::cout << GridLogMessage<< "Calling Laplacian" <<std::endl;
  LatticeComplex in(&VertexGrid);
  LatticeComplex out(&VertexGrid);
  gaussian(vRNG,in);
  Support.Laplacian(in,out);
  
  Grid_finalize();
}
