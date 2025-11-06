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


// FIXME -- the ROTATE code is CPU vector only
// Need coalescedReadPermute to support rotations

///////////////////////////////////////////////////////////
// 3 links per edge site in icos-plane; 3.10.L^2xLt
// 1 temporal per vertex perp to plane; (10.L^2+2)xLt
//
// This makes the HMC complicated -- the field type integrated is not simply U but a collection of fields
// and momenta with different support. Easiest to view as a "bunch of links", but breaks current code designs
// unless we have a way to split a bigger gauge field with "empty" or unused slots; probably the cleanest.
///////////////////////////////////////////////////////////
template<typename vtype> using iIcosahedralLorentzColourMatrix       = iVector<iScalar<iMatrix<vtype, Nc> >, NumIcosahedralPolarizations> ;
typedef iIcosahedralLorentzColourMatrix<Complex> IcosahedralLorentzColourMatrix;
typedef iIcosahedralLorentzColourMatrix<vComplex >  vIcosahedralLorentzColourMatrix;
typedef Lattice<vIcosahedralLorentzColourMatrix>   LatticeIcosahedralLorentzColourMatrix;

class IcosahedralGimpl
{
public:
  typedef LatticeIcosahedralLorentzColourMatrix  GaugeField;
  typedef LatticeDoubledGaugeField               DoubledGaugeField;
  typedef LatticeColourMatrix                    GaugeLinkField;
  typedef LatticeComplex              ComplexField;
};
  
template< class Gimpl>
class IcosahedralSupport {
public:
  // Move to inherit types macro as before
  typedef typename Gimpl::GaugeField GaugeField;
  typedef typename Gimpl::GaugeLinkField GaugeLinkField;
  typedef typename Gimpl::ComplexField ComplexField;
  typedef typename Gimpl::DoubledGaugeField DoubledGaugeField;
  //
  GridBase *VertexGrid;
  GridBase *EdgeGrid;

  IcosahedralStencil FaceStencil;
  IcosahedralStencil NNee; // edge neighbours with edge domain
  IcosahedralStencil NNev; // vertex neighbours but in edge domain
  IcosahedralStencil NNvv; // vertex neighbours with vertex domain
  IcosahedralStencil NNve; // edge neighbours with vertex domain

  IcosahedralStencil TimePlaqStencil;
  IcosahedralStencil TimeStapleStencil;
  
  IcosahedralSupport(GridBase *_VertexGrid,GridBase *_EdgeGrid)
    : FaceStencil (_EdgeGrid,_VertexGrid),
      NNee(_EdgeGrid,_VertexGrid),
      NNev(_EdgeGrid,_VertexGrid),
      NNve(_EdgeGrid,_VertexGrid),
      NNvv(_EdgeGrid,_VertexGrid),
      TimePlaqStencil(_EdgeGrid,_VertexGrid),
      TimeStapleStencil(_EdgeGrid,_VertexGrid),
      VertexGrid(_VertexGrid), EdgeGrid(_EdgeGrid)
  {
    FaceStencil.FaceStencil();

    // true/false is "vertexInput, VertexOutput"
    NNee.NearestNeighbourStencil(false,false);// edge input + output      ; used by face stencil
    NNev.NearestNeighbourStencil(true,false); // vertex input, edge ouput ; used by gauge transform
    NNvv.NearestNeighbourStencil(true,true);  // vertex input + output    ; used by Laplacian
    NNve.NearestNeighbourStencil(false,true); // edge input, vertex output; used by double store
    TimePlaqStencil.TemporalPlaquetteStencil();
    TimeStapleStencil.TemporalStapleStencil();
  }

  ////////////////////////////////////////////////////////////////////////////////////
  // Gauge Link field GT is the gauge transform and lives on the VERTEX field
  ////////////////////////////////////////////////////////////////////////////////////
  void ForwardTriangles(GaugeField &Umu,ComplexField &plaq1,ComplexField &plaq2)
  {
    {
      autoView(Umu_v,Umu,AcceleratorRead);
      autoView(plaq1_v,plaq1,AcceleratorWrite);
      autoView(plaq2_v,plaq2,AcceleratorWrite);
      autoView(stencil_v,FaceStencil,AcceleratorRead);

      accelerator_for(ss,EdgeGrid->oSites(),vComplex::Nsimd(),{
	
	auto Lx = Umu_v(ss)(IcosahedronPatchX);
	auto Ly = Umu_v(ss)(IcosahedronPatchY);
	auto Ld = Umu_v(ss)(IcosahedronPatchDiagonal);

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
      
      const int ent_Xp = 0;
      const int ent_Yp = 1;
      const int ent_Xm = 3;
      const int ent_Ym = 4;

      accelerator_for(ss,EdgeGrid->oSites(),vComplex::Nsimd(),{

	  // Three forward links from this site
	auto Lx = Umu_v(ss)(IcosahedronPatchX);
	auto Ly = Umu_v(ss)(IcosahedronPatchY);
	auto Ld = Umu_v(ss)(IcosahedronPatchDiagonal);

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
  
  template<class MatterField>
  void Laplacian(MatterField &in,MatterField &out)
  {
    autoView(out_v,out,AcceleratorWrite);
    autoView(in_v,in,AcceleratorRead);
    autoView(stencil_v,NNvv,AcceleratorRead);

    const int np = NNvv._npoints;

    const int ent_Xp = 0;
    const int ent_Yp = 1;
    const int ent_Dp = 2;
    const int ent_Xm = 3;
    const int ent_Ym = 4;
    const int ent_Dm = 5;
    
    accelerator_for(ss,VertexGrid->oSites(),vComplex::Nsimd(),{
	  
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
  /*
   * Temporarily compute all without summing p and m directions, and XYD for the temporal links
   * Do this so we can compute and check all link x staple loops
  void TemporalStaplesOld(GaugeLinkField &Ut,
			  GaugeField &Umu,
			  GaugeLinkField &stapleXTp,
			  GaugeLinkField &stapleYTp,
			  GaugeLinkField &stapleDTp,
			  GaugeLinkField &stapleXTm,
			  GaugeLinkField &stapleYTm,
			  GaugeLinkField &stapleDTm,
			  GaugeLinkField &stapleTXp, // these 5 are on vertex graph
			  GaugeLinkField &stapleTYp, // and contain the north south pole temporal link staples
			  GaugeLinkField &stapleTDp,
			  GaugeLinkField &stapleTXm,
			  GaugeLinkField &stapleTYm,
			  GaugeLinkField &stapleTDm)
  {
  
      autoView(Umu_v,Umu,AcceleratorRead);
      autoView(Ut_v,Ut,AcceleratorRead);

      autoView(stapleXTp_v,stapleXTp,AcceleratorWrite);
      autoView(stapleYTp_v,stapleYTp,AcceleratorWrite);
      autoView(stapleDTp_v,stapleDTp,AcceleratorWrite);

      autoView(stapleXTm_v,stapleXTm,AcceleratorWrite);
      autoView(stapleYTm_v,stapleYTm,AcceleratorWrite);
      autoView(stapleDTm_v,stapleDTm,AcceleratorWrite);

      autoView(stapleTXp_v,stapleTXp,AcceleratorWrite);
      autoView(stapleTYp_v,stapleTYp,AcceleratorWrite);
      autoView(stapleTDp_v,stapleTDp,AcceleratorWrite);

      autoView(stapleTXm_v,stapleTXm,AcceleratorWrite);
      autoView(stapleTYm_v,stapleTYm,AcceleratorWrite);
      autoView(stapleTDm_v,stapleTDm,AcceleratorWrite);

      autoView(stencil_v,TimeStapleStencil,AcceleratorRead);

      const int np = TimeStapleStencil._npoints;
      
      Integer ent_Tp   = 0;
      Integer ent_Xp   = 1;
      Integer ent_Yp   = 2;
      Integer ent_Dp   = 3;
      Integer ent_Tm   = 4;
      Integer ent_TmXp = 5;
      Integer ent_TmYp = 6;
      Integer ent_TmDp = 7;
      Integer ent_Xm   = 8;
      Integer ent_Ym   = 9;
      Integer ent_Dm   = 10;
      Integer ent_XmTp = 11;
      Integer ent_YmTp = 12;
      Integer ent_DmTp = 13; 

      accelerator_for(ss,EdgeGrid->oSites(),vComplex::Nsimd(),{

	// Four forward links from this site
	auto Lx = Umu_v(ss)(IcosahedronPatchX);
	auto Ly = Umu_v(ss)(IcosahedronPatchY);
	auto Ld = Umu_v(ss)(IcosahedronPatchDiagonal);
	auto Lt = Ut_v(ss)();
	IcosahedralStencilEntry *SE;
	int64_t s;
	int p;
	///////////////////////////////////////////////////////////////////
	// Load all links required and multiple round loops
	///////////////////////////////////////////////////////////////////
	SE = stencil_v.GetEntry(ent_Tp,ss);
	s  = SE->_offset;
	p  = SE->_permute;
	auto Lx_at_tp = Umu_v(s)(IcosahedronPatchX);
	auto Ly_at_tp = Umu_v(s)(IcosahedronPatchY);
	auto Ld_at_tp = Umu_v(s)(IcosahedronPatchDiagonal);
	if ( p ) {
	  #ifdef GRID_SIMT
	  #error GPU not supported yet
	  #endif
	  rotate(Lx_at_tp,Lx_at_tp,1);
	  rotate(Ly_at_tp,Ly_at_tp,1);
	  rotate(Ld_at_tp,Ld_at_tp,1);
	}
	SE = stencil_v.GetEntry(ent_Xp,ss);
	s  = SE->_offset;
	auto Lt_at_xp = Ut_v(s)();

	SE = stencil_v.GetEntry(ent_Yp,ss);
	s  = SE->_offset;
	auto Lt_at_yp = Ut_v(s)();

	SE = stencil_v.GetEntry(ent_Dp,ss);
	s  = SE->_offset;
	auto Lt_at_dp = Ut_v(s)();

	Coordinate Coor;
	EdgeGrid->oCoorFromOindex(Coor,ss);
	
	// Lines start at origin and progress to top of link
	// close as gauge invariant if multiplied from left by link
	coalescedWrite(stapleXTp_v[ss](),Lt_at_xp*adj(Lx_at_tp)*adj(Lt));
	coalescedWrite(stapleYTp_v[ss](),Lt_at_yp*adj(Ly_at_tp)*adj(Lt));
	coalescedWrite(stapleDTp_v[ss](),Lt_at_dp*adj(Ld_at_tp)*adj(Lt));

	SE = stencil_v.GetEntry(ent_Tm,ss);
	s = SE->_offset;
	p = SE->_permute;
	auto Lt_at_tm = Ut_v(s)();
	auto Lx_at_tm = Umu_v(s)(IcosahedronPatchX);
	auto Ly_at_tm = Umu_v(s)(IcosahedronPatchY);
	auto Ld_at_tm = Umu_v(s)(IcosahedronPatchDiagonal);
	const int Nsimdm1 = vComplex::Nsimd()-1;

	if ( p ) {
	  rotate(Lt_at_tm,Lt_at_tm,Nsimdm1);
	  rotate(Lx_at_tm,Lx_at_tm,Nsimdm1);
	  rotate(Ly_at_tm,Ly_at_tm,Nsimdm1);
	  rotate(Ld_at_tm,Ld_at_tm,Nsimdm1);
	}

	SE = stencil_v.GetEntry(ent_TmXp,ss);
	s  = SE->_offset;
	p  = SE->_permute;
	auto Lt_at_tmxp = Ut_v(s)();
	
	SE = stencil_v.GetEntry(ent_TmYp,ss);
	s  = SE->_offset;
	auto Lt_at_tmyp = Ut_v(s)();
	
	SE = stencil_v.GetEntry(ent_TmDp,ss);
	s  = SE->_offset;
	auto Lt_at_tmdp = Ut_v(s)();
	
	if ( p ) {
	  rotate(Lt_at_tmxp,Lt_at_tmxp,Nsimdm1);
	  rotate(Lt_at_tmyp,Lt_at_tmyp,Nsimdm1);
	  rotate(Lt_at_tmdp,Lt_at_tmdp,Nsimdm1);
	}

	coalescedWrite(stapleXTm_v[ss](),adj(Lt_at_tmxp)*adj(Lx_at_tm)*Lt_at_tm);
	coalescedWrite(stapleYTm_v[ss](),adj(Lt_at_tmyp)*adj(Ly_at_tm)*Lt_at_tm);
	coalescedWrite(stapleDTm_v[ss](),adj(Lt_at_tmdp)*adj(Ld_at_tm)*Lt_at_tm);
	
	coalescedWrite(stapleTXp_v[ss](),Lx_at_tp*adj(Lt_at_xp)*adj(Lx));
	coalescedWrite(stapleTYp_v[ss](),Ly_at_tp*adj(Lt_at_yp)*adj(Ly));
	coalescedWrite(stapleTDp_v[ss](),Ld_at_tp*adj(Lt_at_dp)*adj(Ld));

	int pol;
	SE = stencil_v.GetEntry(ent_Xm,ss);
	s  = SE->_offset;
	pol= SE->_polarisation;
	auto Lt_at_xm = Ut_v(s)();
	auto Lx_at_xm = Umu_v(s)(pol);

	SE = stencil_v.GetEntry(ent_XmTp,ss);
	s  = SE->_offset;
	p  = SE->_permute;
	pol= SE->_polarisation;
	auto Lx_at_xmtp = Umu_v(s)(pol);
	
	SE = stencil_v.GetEntry(ent_Ym,ss);
	s  = SE->_offset;
	pol= SE->_polarisation;
	auto Lt_at_ym = Ut_v(s)();
	auto Ly_at_ym = Umu_v(s)(pol);

	SE = stencil_v.GetEntry(ent_YmTp,ss);
	s  = SE->_offset;
	p  = SE->_permute;
	pol= SE->_polarisation;
	auto Ly_at_ymtp = Umu_v(s)(pol);

	if ( p ) {
	  rotate(Ly_at_ymtp,Ly_at_ymtp,1);
	  rotate(Lx_at_xmtp,Lx_at_xmtp,1);
	}

	coalescedWrite(stapleTXm_v[ss](),adj(Lx_at_xmtp)*adj(Lt_at_xm)*Lx_at_xm);
	coalescedWrite(stapleTYm_v[ss](),adj(Ly_at_ymtp)*adj(Lt_at_ym)*Ly_at_ym);
	
	SE = stencil_v.GetEntry(ent_Dm,ss);
	int missingLink = SE->_missing_link;
	s  = SE->_offset;
	pol= SE->_polarisation;
	if ( missingLink ) {
	  auto zero = Lx;
	  zero = Zero();
	  coalescedWrite(stapleTDm_v[ss](),zero);
	} else {
	  auto Lt_at_dm = Ut_v(s)();
	  auto Ld_at_dm = Umu_v(s)(pol);
	  SE = stencil_v.GetEntry(ent_DmTp,ss);
	  s  = SE->_offset;
	  auto Ld_at_dmtp = Umu_v(s)(pol);
	  if ( p ) {
	    rotate(Ld_at_dmtp,Ld_at_dmtp,1);
	  }
	  coalescedWrite(stapleTDm_v[ss](),adj(Ld_at_dmtp)*adj(Lt_at_dm)*Ld_at_dm);
	}
	
      });

      const int ent_h0   = 0;
      const int ent_h0tp = 1;
      const int ent_h1   = 2;
      const int ent_h1tp = 3;
      const int ent_h2   = 4;
      const int ent_h2tp = 5;
      const int ent_h3   = 6;
      const int ent_h3tp = 7;
      const int ent_h4   = 8;
      const int ent_h4tp = 9;
      
      uint64_t cart_sites = EdgeGrid->oSites();
      uint64_t pole_sites = VertexGrid->NorthPoleOsites()
	                  + VertexGrid->SouthPoleOsites();
      
#define LINK(ENT_0,ENT_1) \
	{SE = stencil_v.GetEntry(ENT_0,ss);\
	  s   = SE->_offset;		   \
	pol = SE->_polarisation;	   \
					   \
	auto Lt    = Ut_v(ss)();	   \
	auto Li    = Umu_v(s)(pol);	   \
	auto Lth   = Ut_v(s)();		   \
					   \
	SE=stencil_v.GetEntry(ENT_1,ss);   \
	s   = SE->_offset;		   \
	p   = SE->_permute;		   \
	auto Litp  = Umu_v(s)(pol);	   \
	if ( p ) {			   \
	  rotate(Litp,Litp,1);		   \
	}				   \
					   \
	auto stap = adj(Litp)*adj(Lth)*Li;		\
	std::cout << ss<<" "<<ENT_0<<" plaq "<<Lt*stap<<std::endl;}

      accelerator_for(sss,pole_sites,vComplex::Nsimd(),{

	IcosahedralStencilEntry *SE;
	int64_t s;
	int p;
	int pol;

	uint64_t ss = sss+cart_sites;

	LINK(ent_h0,ent_h0tp);
	LINK(ent_h1,ent_h1tp);
	LINK(ent_h2,ent_h2tp);
	LINK(ent_h3,ent_h3tp);
	LINK(ent_h4,ent_h4tp);
	
	});
#undef LINK
  }
   */


  void TemporalStaples(GaugeLinkField &Ut,
		       GaugeField &Umu,
		       GaugeLinkField &stapleX,
		       GaugeLinkField &stapleY,
		       GaugeLinkField &stapleD,
		       GaugeLinkField &stapleT)
  {
      autoView(Umu_v,Umu,AcceleratorRead);
      autoView(Ut_v,Ut,AcceleratorRead);

      autoView(stapleX_v,stapleX,AcceleratorWrite);
      autoView(stapleY_v,stapleY,AcceleratorWrite);
      autoView(stapleD_v,stapleD,AcceleratorWrite);
      autoView(stapleT_v,stapleT,AcceleratorWrite);

      autoView(stencil_v,TimeStapleStencil,AcceleratorRead);

      const int np = TimeStapleStencil._npoints;
      
      Integer ent_Tp   = 0;
      Integer ent_Xp   = 1;
      Integer ent_Yp   = 2;
      Integer ent_Dp   = 3;
      Integer ent_Tm   = 4;
      Integer ent_TmXp = 5;
      Integer ent_TmYp = 6;
      Integer ent_TmDp = 7;
      Integer ent_Xm   = 8;
      Integer ent_Ym   = 9;
      Integer ent_Dm   = 10;
      Integer ent_XmTp = 11;
      Integer ent_YmTp = 12;
      Integer ent_DmTp = 13; 

      accelerator_for(ss,EdgeGrid->oSites(),vComplex::Nsimd(),{
	
	// Four forward links from this site
	auto Lx = Umu_v(ss)(IcosahedronPatchX);
	auto Ly = Umu_v(ss)(IcosahedronPatchY);
	auto Ld = Umu_v(ss)(IcosahedronPatchDiagonal);
	auto Lt = Ut_v(ss)();
	IcosahedralStencilEntry *SE;
	int64_t s;
	int p;
	///////////////////////////////////////////////////////////////////
	// Load all links required and multiple round loops
	///////////////////////////////////////////////////////////////////
	SE = stencil_v.GetEntry(ent_Tp,ss);
	s  = SE->_offset;
	p  = SE->_permute;
	auto Lx_at_tp = Umu_v(s)(IcosahedronPatchX);
	auto Ly_at_tp = Umu_v(s)(IcosahedronPatchY);
	auto Ld_at_tp = Umu_v(s)(IcosahedronPatchDiagonal);
	if ( p ) {
	  rotate(Lx_at_tp,Lx_at_tp,1);
	  rotate(Ly_at_tp,Ly_at_tp,1);
	  rotate(Ld_at_tp,Ld_at_tp,1);
	}
	SE = stencil_v.GetEntry(ent_Xp,ss);
	s  = SE->_offset;
	auto Lt_at_xp = Ut_v(s)();

	SE = stencil_v.GetEntry(ent_Yp,ss);
	s  = SE->_offset;
	auto Lt_at_yp = Ut_v(s)();

	SE = stencil_v.GetEntry(ent_Dp,ss);
	s  = SE->_offset;
	auto Lt_at_dp = Ut_v(s)();

	Coordinate Coor;
	EdgeGrid->oCoorFromOindex(Coor,ss);
	

	SE = stencil_v.GetEntry(ent_Tm,ss);
	s = SE->_offset;
	p = SE->_permute;
	auto Lt_at_tm = Ut_v(s)();
	auto Lx_at_tm = Umu_v(s)(IcosahedronPatchX);
	auto Ly_at_tm = Umu_v(s)(IcosahedronPatchY);
	auto Ld_at_tm = Umu_v(s)(IcosahedronPatchDiagonal);
	const int Nsimdm1 = vComplex::Nsimd()-1;

	if ( p ) {
	  rotate(Lt_at_tm,Lt_at_tm,Nsimdm1);
	  rotate(Lx_at_tm,Lx_at_tm,Nsimdm1);
	  rotate(Ly_at_tm,Ly_at_tm,Nsimdm1);
	  rotate(Ld_at_tm,Ld_at_tm,Nsimdm1);
	}

	SE = stencil_v.GetEntry(ent_TmXp,ss);
	s  = SE->_offset;
	p  = SE->_permute;
	auto Lt_at_tmxp = Ut_v(s)();
	
	SE = stencil_v.GetEntry(ent_TmYp,ss);
	s  = SE->_offset;
	auto Lt_at_tmyp = Ut_v(s)();
	
	SE = stencil_v.GetEntry(ent_TmDp,ss);
	s  = SE->_offset;
	auto Lt_at_tmdp = Ut_v(s)();
	
	if ( p ) {
	  rotate(Lt_at_tmxp,Lt_at_tmxp,Nsimdm1);
	  rotate(Lt_at_tmyp,Lt_at_tmyp,Nsimdm1);
	  rotate(Lt_at_tmdp,Lt_at_tmdp,Nsimdm1);
	}

	// Lines start at origin and progress to top of link
	// close as gauge invariant if multiplied from left by link
	
	coalescedWrite(stapleX_v[ss](),Lt_at_xp*adj(Lx_at_tp)*adj(Lt) + adj(Lt_at_tmxp)*adj(Lx_at_tm)*Lt_at_tm);
	coalescedWrite(stapleY_v[ss](),Lt_at_yp*adj(Ly_at_tp)*adj(Lt) + adj(Lt_at_tmyp)*adj(Ly_at_tm)*Lt_at_tm);
	coalescedWrite(stapleD_v[ss](),Lt_at_dp*adj(Ld_at_tp)*adj(Lt) + adj(Lt_at_tmdp)*adj(Ld_at_tm)*Lt_at_tm);

	int pol;
	SE = stencil_v.GetEntry(ent_Xm,ss);
	s  = SE->_offset;
	pol= SE->_polarisation;
	auto Lt_at_xm = Ut_v(s)();
	auto Lx_at_xm = Umu_v(s)(pol);

	SE = stencil_v.GetEntry(ent_XmTp,ss);
	s  = SE->_offset;
	p  = SE->_permute;
	pol= SE->_polarisation;
	auto Lx_at_xmtp = Umu_v(s)(pol);
	
	SE = stencil_v.GetEntry(ent_Ym,ss);
	s  = SE->_offset;
	pol= SE->_polarisation;
	auto Lt_at_ym = Ut_v(s)();
	auto Ly_at_ym = Umu_v(s)(pol);

	SE = stencil_v.GetEntry(ent_YmTp,ss);
	s  = SE->_offset;
	p  = SE->_permute;
	pol= SE->_polarisation;
	auto Ly_at_ymtp = Umu_v(s)(pol);

	if ( p ) {
	  rotate(Ly_at_ymtp,Ly_at_ymtp,1);
	  rotate(Lx_at_xmtp,Lx_at_xmtp,1);
	}

	auto stapT = Lx;
	stapT= Lx_at_tp*adj(Lt_at_xp)*adj(Lx);
	stapT= stapT + Ly_at_tp*adj(Lt_at_yp)*adj(Ly);
	stapT= stapT + Ld_at_tp*adj(Lt_at_dp)*adj(Ld);
	stapT= stapT + adj(Lx_at_xmtp)*adj(Lt_at_xm)*Lx_at_xm;
	stapT= stapT + adj(Ly_at_ymtp)*adj(Lt_at_ym)*Ly_at_ym;

	SE = stencil_v.GetEntry(ent_Dm,ss);
	int missingLink = SE->_missing_link;
	s  = SE->_offset;
	pol= SE->_polarisation;
	if ( !missingLink ) {

	  auto Lt_at_dm = Ut_v(s)();
	  auto Ld_at_dm = Umu_v(s)(pol);
	  SE = stencil_v.GetEntry(ent_DmTp,ss);
	  s  = SE->_offset;
	  auto Ld_at_dmtp = Umu_v(s)(pol);
	  if ( p ) {
	    rotate(Ld_at_dmtp,Ld_at_dmtp,1);
	  }
	  stapT = stapT + adj(Ld_at_dmtp)*adj(Lt_at_dm)*Ld_at_dm;
	}

        coalescedWrite(stapleT_v[ss](),stapT);
	
	});

      const int ent_h0   = 0;
      const int ent_h0tp = 1;
      const int ent_h1   = 2;
      const int ent_h1tp = 3;
      const int ent_h2   = 4;
      const int ent_h2tp = 5;
      const int ent_h3   = 6;
      const int ent_h3tp = 7;
      const int ent_h4   = 8;
      const int ent_h4tp = 9;
      
      uint64_t cart_sites = EdgeGrid->oSites();
      uint64_t pole_sites = VertexGrid->NorthPoleOsites()
	                  + VertexGrid->SouthPoleOsites();
      

#define LINK(ENT_0,ENT_1) \
	SE = stencil_v.GetEntry(ENT_0,ss); \
	  s   = SE->_offset;		   \
	pol = SE->_polarisation;	   \
					   \
	Li    = Umu_v(s)(pol);		   \
	Lth   = Ut_v(s)();		   \
					   \
	SE=stencil_v.GetEntry(ENT_1,ss);   \
	s   = SE->_offset;		   \
	p   = SE->_permute;		   \
	Litp  = Umu_v(s)(pol);		   \
	if ( p ) {			   \
	  rotate(Litp,Litp,1);		   \
	}
      
#define LINK0(ENT_0,ENT_1) LINK(ENT_0,ENT_1); auto stapT = adj(Litp)*adj(Lth)*Li;
#define LINKn(ENT_0,ENT_1) LINK(ENT_0,ENT_1); stapT = stapT + adj(Litp)*adj(Lth)*Li;

      accelerator_for(sss,pole_sites,vComplex::Nsimd(),{

	IcosahedralStencilEntry *SE;
	int64_t s;
	int p;
	int pol;

	uint64_t ss = sss+cart_sites;

	auto Lt  = Ut_v(ss)();
	auto Lth = Lt;
	auto Li  = Lt;
	auto Litp= Lt;

	LINK0(ent_h0,ent_h0tp);
	LINKn(ent_h1,ent_h1tp);
	LINKn(ent_h2,ent_h2tp);
	LINKn(ent_h3,ent_h3tp);
	LINKn(ent_h4,ent_h4tp);
	
	coalescedWrite(stapleT_v[ss](),stapT);
	
	});
#undef LINK
#undef LINK0
#undef LINKn
  }

  
  template<class Field>
  void CopyNonPoles(Field &in,Field &out)
  {
    autoView(out_v,out,AcceleratorWrite);
    autoView(in_v,in,AcceleratorRead);
    
    accelerator_for(ss,EdgeGrid->oSites(),vComplex::Nsimd(),{
	  
    auto i    = in_v(ss)();
     coalescedWrite(out_v[ss](),i);
    });
  }

  template<class MatterField>
  void CovariantLaplacian(MatterField &in,MatterField &out,DoubledGaugeField &Uds)
  {
    autoView(out_v,out,AcceleratorWrite);
    autoView(in_v,in,AcceleratorRead);
    autoView(U_v,Uds,AcceleratorRead);
    autoView(stencil_v,NNvv,AcceleratorRead);

    const int np = NNvv._npoints;

    const int ent_Xp = 0;
    const int ent_Yp = 1;
    const int ent_Dp = 2;
    const int ent_Xm = 3;
    const int ent_Ym = 4;
    const int ent_Dm = 5;
    const int ent_Tp = 6;
    const int ent_Tm = 7;

    
    accelerator_for(ss,VertexGrid->oSites(),vComplex::Nsimd(),{
	  
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

      inxp = U_v(ss)(0)*inxp;
      inyp = U_v(ss)(1)*inyp;
      indp = U_v(ss)(2)*indp;
      inxm = U_v(ss)(3)*inxm;
      inym = U_v(ss)(4)*inym;
      
      auto o    = i;

      if ( missingLink ) {
	o = (1.0/5.0)*(inxp+inyp+indp+inxm+inym)-i;
      } else {
	indm = U_v(ss)(5)*indm;
	o = (1.0/6.0)*(inxp+inyp+indp+inxm+inym+indm)-i;
      }
      
      coalescedWrite(out_v[ss](),o);
      });
  }
  
  void DoubleStore(GaugeField &U,DoubledGaugeField &Uds)
  {
    assert(U.Grid()==EdgeGrid);
    assert(Uds.Grid()==VertexGrid);
    autoView(Uds_v,Uds,AcceleratorWrite);
    autoView(U_v,U,AcceleratorRead);
    autoView(stencil_v,NNvv,AcceleratorRead);

    // Vertex result
    // Edge valued input
    // Might need an extra case (?)
    const int np = NNvv._npoints;

    const int ent_Xm = 3;
    const int ent_Ym = 4;
    const int ent_Dm = 5;
    
    const int ent_Tp = 6;
    const int ent_Tm = 7;
    
    accelerator_for(ss,VertexGrid->CartesianOsites(),vComplex::Nsimd(),{
	
      // Three local links
      {
	auto Lx = U_v(ss)(IcosahedronPatchX);
	auto Ly = U_v(ss)(IcosahedronPatchY);
	auto Ld = U_v(ss)(IcosahedronPatchDiagonal);
	coalescedWrite(Uds_v[ss](IcosahedronPatchX),Lx);
	coalescedWrite(Uds_v[ss](IcosahedronPatchY),Ly);
	coalescedWrite(Uds_v[ss](IcosahedronPatchDiagonal),Ld);
      }	
      // Three backwards links
      {
	auto SE = stencil_v.GetEntry(ent_Xm,ss);
	auto pol = SE->_polarisation;
	auto s   = SE->_offset;
	int pol1 = IcosahedronPatchX;
	//
	// Xm link is given diagonal unless HemiPatch changes in Northern hemisphere
	// But here for double storing we need the reverse link so
	// return either the Xplus or Diag plus direction
	//
	if ( pol != IcosahedronPatchDiagonal ) {
	  pol1 = IcosahedronPatchDiagonal;
	}
	auto Lx_at_xm = U_v(s)(pol1);
	Lx_at_xm = adj(Lx_at_xm);
	coalescedWrite(Uds_v[ss](3),Lx_at_xm );
      }
      {
	auto SE   = stencil_v.GetEntry(ent_Ym,ss);
	auto pol  = SE->_polarisation;
	auto s    = SE->_offset;
	//
	// Ym link is given diagonal unless HemiPatch changes in the Southern hemisphere
	// But here for double storing we need the reverse link so
	// return either the Yplus or Diag plus direction
	// 
	int pol1 = IcosahedronPatchY;
	if ( pol != IcosahedronPatchDiagonal ) {
	  pol1 = IcosahedronPatchDiagonal;
	}
	auto Ly_at_ym = U_v(s)(pol1);
	Ly_at_ym = adj(Ly_at_ym);
	coalescedWrite(Uds_v[ss](4),Ly_at_ym );
      }
      int missingLink;
      auto Link =  adj(U_v(ss)(0));
      Link=Zero();
      {      
	auto SE = stencil_v.GetEntry(ent_Dm,ss);
	auto s  = SE->_offset;
	// Dm link given the correct value for this use case
	auto pol= SE->_polarisation;
	missingLink = SE->_missing_link;
	if ( ! missingLink ) {
	  auto Ld_at_dm  =  U_v(s)(pol);
	  Ld_at_dm = adj(Ld_at_dm);
	  coalescedWrite(Uds_v[ss](5),Ld_at_dm );
	} else {
	  coalescedWrite(Uds_v[ss](5),Link);
	}
      }
      coalescedWrite(Uds_v[ss](6),Link);
      coalescedWrite(Uds_v[ss](7),Link);
    });
    auto pole_sites = VertexGrid->oSites() - VertexGrid->CartesianOsites();
    auto pole_offset= VertexGrid->CartesianOsites();
    
    accelerator_for(ss,pole_sites,vComplex::Nsimd(),{
	for(int p=0;p<HemiPatches;p++){ // Neighbours of each pole
	auto SE = stencil_v.GetEntry(p,pole_offset+ss);
	auto s  = SE->_offset;
	auto pol= SE->_polarisation;
	auto Link =  adj(U_v(s)(pol));
	coalescedWrite(Uds_v[pole_offset+ss](p),Link);
	}
	auto Link =  adj(U_v(0)(0));
	Link=Zero();
	coalescedWrite(Uds_v[pole_offset+ss](5),Link);
	coalescedWrite(Uds_v[pole_offset+ss](6),Link);
	coalescedWrite(Uds_v[pole_offset+ss](7),Link);
    });
  }
  void  GaugeTransform(GaugeLinkField &gt, GaugeField &Umu, GaugeLinkField &Ut)
  {
    autoView(Umu_v,Umu,AcceleratorWrite);
    autoView(Ut_v,Ut,AcceleratorWrite);
    autoView(g_v,gt,AcceleratorRead);
    autoView(stencil_v ,NNev,AcceleratorRead);
    autoView(stencil_vv,NNvv,AcceleratorRead);
    
    const int ent_Xp = 0;
    const int ent_Yp = 1;
    const int ent_Dp = 2;

    accelerator_for(ss,EdgeGrid->oSites(),vComplex::Nsimd(),{

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

    const int ent_Tp = 6;
    accelerator_for(ss,VertexGrid->oSites(),vComplex::Nsimd(),{
      auto lt = Ut_v(ss)();
      auto g  = g_v(ss)();
      auto SE = stencil_vv.GetEntry(ent_Tp,ss);
      uint64_t tp_idx = SE->_offset;
      auto permute    = SE->_permute;
      auto gtp = g_v(tp_idx)();
      //      std::cout << ss<<" tp "<<tp_idx<<" permute "<<(int)SE->_permute<<" : g"<<g<<" l "<<lt<<" gtp "<<gtp<<std::endl;
      if ( permute ) rotate(gtp,gtp,1);
      lt = g*lt*adj(gtp);
      coalescedWrite(Ut_v[ss](),lt);
    });    
  }
  void  TemporalPlaquette(GaugeLinkField &Ut, GaugeField &Uj,ComplexField &plaq1,ComplexField &plaq2,ComplexField &plaq3)
  {
    autoView(Uj_v,Uj,AcceleratorRead);
    autoView(Ut_v,Ut,AcceleratorRead);
    autoView(plaq1_v,plaq1,AcceleratorWrite);
    autoView(plaq2_v,plaq2,AcceleratorWrite);
    autoView(plaq3_v,plaq3,AcceleratorWrite);
    autoView(stencil_v,TimePlaqStencil,AcceleratorRead);

    const int np = TimePlaqStencil._npoints;
    
    const int ent_Tp = 0;
    const int ent_Xp = 1;
    const int ent_Yp = 2;
    const int ent_Dp = 3;

    accelerator_for(ss,EdgeGrid->oSites(),vComplex::Nsimd(),{

      // Three forward links from this site
      auto Lx = Uj_v(ss)(IcosahedronPatchX);
      auto Ly = Uj_v(ss)(IcosahedronPatchY);
      auto Ld = Uj_v(ss)(IcosahedronPatchDiagonal);
      auto Lt = Ut_v(ss)();

      // Three forward links from this site at Tp
      auto SE = stencil_v.GetEntry(ent_Tp,ss);
      uint64_t tp_idx = SE->_offset;
      uint64_t permute= SE->_permute;

      auto Ltp_x = Uj_v(tp_idx)(IcosahedronPatchX);
      auto Ltp_y = Uj_v(tp_idx)(IcosahedronPatchY);
      auto Ltp_d = Uj_v(tp_idx)(IcosahedronPatchDiagonal);

      if ( permute ) {
	rotate(Ltp_x,Ltp_x,1); // does in=out work?
	rotate(Ltp_y,Ltp_y,1);
	rotate(Ltp_d,Ltp_d,1);
      }
      
      SE = stencil_v.GetEntry(ent_Xp,ss);
      uint64_t xp_idx = SE->_offset;

      SE = stencil_v.GetEntry(ent_Yp,ss);
      uint64_t yp_idx = SE->_offset;

      SE = stencil_v.GetEntry(ent_Dp,ss);
      uint64_t dp_idx = SE->_offset;

      auto Lxp_t = Ut_v(xp_idx)();
      auto Lyp_t = Ut_v(yp_idx)();
      auto Ldp_t = Ut_v(dp_idx)();
      
      coalescedWrite(plaq1_v[ss](),trace(Lx*Lxp_t*adj(Ltp_x)*adj(Lt)));
      coalescedWrite(plaq2_v[ss](),trace(Ly*Lyp_t*adj(Ltp_y)*adj(Lt)));
      coalescedWrite(plaq3_v[ss](),trace(Ld*Ldp_t*adj(Ltp_d)*adj(Lt)));

      /*
      std::cout << ss <<"----------------"<<std::endl;
      std::cout << ss <<" 1 "<<plaq1_v[ss]<<std::endl;
      std::cout << ss <<" 2 "<<plaq2_v[ss]<<std::endl;
      std::cout << ss <<" 3 "<<plaq3_v[ss]<<std::endl;
      std::cout << ss <<" xp_idx "<<xp_idx<<std::endl;
      std::cout << ss <<" yp_idx "<<yp_idx<<std::endl;
      std::cout << ss <<" dp_idx "<<dp_idx<<std::endl;
      std::cout << ss <<" tp_idx "<<dp_idx<<std::endl;
      std::cout << ss <<" Lx "<<Lx<<std::endl;
      std::cout << ss <<" Ly "<<Ly<<std::endl;
      std::cout << ss <<" Ld "<<Ld<<std::endl;
      std::cout << ss <<" Lt    "<<Lt<<std::endl;
      std::cout << ss <<" Ldp_t "<<Ldp_t<<std::endl;
      std::cout << ss <<" Lxp_t "<<Lxp_t<<std::endl;
      std::cout << ss <<" Lyp_t "<<Lyp_t<<std::endl;
      std::cout << ss <<" Ltp_x "<<Ltp_x<<std::endl;
      std::cout << ss <<" Ltp_y "<<Ltp_y<<std::endl;
      std::cout << ss <<" Ltp_d "<<Ltp_d<<std::endl;
      */
      });
  }
};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  //  const int L=8;
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
  typedef IcosahedralGimpl::ComplexField ComplexField;
  typedef IcosahedralGimpl::GaugeLinkField GaugeLinkField;
  typedef IcosahedralGimpl::GaugeField        GaugeField;
  typedef IcosahedralGimpl::DoubledGaugeField DoubledGaugeField;
  GaugeLinkField Ut(&VertexGrid);
  GaugeLinkField Utnp(&EdgeGrid);
  GaugeField Umu(&EdgeGrid);
  GaugeField Umuck(&EdgeGrid);
  ComplexField Phi(&VertexGrid);
  std::cout << GridLogMessage << " Created two fields "<<std::endl;

  Phi = Zero();
  Umu = Zero();
  std::cout << GridLogMessage << " Zeroed two fields "<<std::endl;

  Complex one (1.0);
  Phi = one;
  Umu = one;
  Ut  = one;
  
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
  ComplexField plaq1(&EdgeGrid);
  ComplexField plaq2(&EdgeGrid);
  ComplexField plaq3(&EdgeGrid);

  ComplexField plaq_ref(&EdgeGrid);
  plaq_ref=1.0;

  Support.ForwardTriangles(Umu,plaq1,plaq2);
  std::cout << GridLogMessage << " plaq1 "<< norm2(plaq1)<<std::endl;
  std::cout << GridLogMessage << " plaq2 "<< norm2(plaq2)<<std::endl;

  std::cout << GridLogMessage << " plaq1 err "<< norm2(plaq1-plaq_ref)<<std::endl;
  std::cout << GridLogMessage << " plaq2 err "<< norm2(plaq2-plaq_ref)<<std::endl;

  // Random gauge xform
  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG  vRNG(&EdgeGrid);   vRNG.SeedFixedIntegers(seeds);
  
  GaugeLinkField g(&VertexGrid);
  LatticeReal gr(&VertexGrid);
  ComplexField gc(&VertexGrid);
  //  gr = Zero();
  gaussian(vRNG,gr);
  Complex ci(0.0,1.0);
  gc = toComplex(gr);
  g=one;
  g = g * exp(ci*gc);
 
  std::cout << GridLogMessage << "****************************************"<<std::endl;
  std::cout << GridLogMessage << " Check plaquette is gauge invariant "<<std::endl;
  std::cout << GridLogMessage << "****************************************"<<std::endl;
  std::cout << GridLogMessage << " applying gauge transform"<<std::endl;
  Support.GaugeTransform    (g,Umu,Ut);
  std::cout << GridLogMessage << " applied gauge transform "<<std::endl;
  std::cout << GridLogMessage << " recalculating plaquette "<<std::endl;
  Support.ForwardTriangles(Umu,plaq1,plaq2);
  std::cout << GridLogMessage << " lower plaq1 "<< norm2(plaq1)<<std::endl;
  std::cout << GridLogMessage << " upper plaq2 "<< norm2(plaq2)<<std::endl;

  std::cout << GridLogMessage << " lower plaq1 err "<< norm2(plaq1-plaq_ref)<<std::endl;
  std::cout << GridLogMessage << " upper plaq2 err "<< norm2(plaq2-plaq_ref)<<std::endl;


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

  std::cout << GridLogMessage << " trace D*StapleXY "<<norm2(trace(linkD * stapleXY))<<std::endl;
  std::cout << GridLogMessage << " err " << norm2(trace(linkD * stapleXY)-plaq_ref)<<std::endl;

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

  std::cout << GridLogMessage<< "Calling Laplacian" <<std::endl;
  LatticeColourVector in(&VertexGrid);
  LatticeColourVector out(&VertexGrid);
  LatticeColourVector gout(&VertexGrid);
  LatticeColourVector gin(&VertexGrid);
  gaussian(vRNG,in);
  Support.Laplacian(in,out);

  std::cout << GridLogMessage<< "Calling double storing gauge field" <<std::endl;
  DoubledGaugeField Uds(&VertexGrid);
  Support.DoubleStore(Umu,Uds);
  
  Support.CovariantLaplacian(in,out,Uds);

  auto ip = innerProduct(out, in);
  std::cout << GridLogMessage<< "Applied covariant laplacian !" <<std::endl;
  /*
   * CovariantLaplacian testing -- check the laplacian is gauge invariant
   *
   * D[U_gt](gF) = g D[U] g^dag g F = g D[U] F
   */
  gout = g*out;
  gin  = g*in;

  Support.GaugeTransform(g,Umu,Ut);
  Support.DoubleStore(Umu,Uds);
  Support.CovariantLaplacian(gin,out,Uds);
  std::cout << GridLogMessage<< "Applied gauge transformed covariant laplacian to transformed vector !" <<std::endl;
  auto ipgt = innerProduct(out, gin);

  std::cout <<GridLogMessage<< "Testing D[U_gt](gF) = g D[U] F : defect is "<<norm2(out-gout)<<std::endl;

  ip = ip - ipgt;
  std::cout <<GridLogMessage<< "Testing F D[U](F) = (gF) D[U_gt] gF : defect is "<<ip<<std::endl;

  /*
   * Temporal plaquettes
   */
  std::cout << GridLogMessage<< " Computing temporal plaquettes"<<std::endl;
  Support.TemporalPlaquette(Ut, Umu,plaq1,plaq2,plaq3);

  
  std::cout << GridLogMessage<<" Computed temporal plaquettes"<<std::endl;
  std::cout << GridLogMessage<<" XT plaq1 "<<norm2(plaq1)<<std::endl;
  std::cout << GridLogMessage<<" YT plaq2 "<<norm2(plaq2)<<std::endl;
  std::cout << GridLogMessage<<" DT plaq3 "<<norm2(plaq3)<<std::endl;

  std::cout << GridLogMessage << " XT plaq1 err "<< norm2(plaq1-plaq_ref)<<std::endl;
  std::cout << GridLogMessage << " YT plaq2 err "<< norm2(plaq2-plaq_ref)<<std::endl;
  std::cout << GridLogMessage << " DT plaq3 err "<< norm2(plaq3-plaq_ref)<<std::endl;

  GaugeLinkField stapleX(&EdgeGrid);
  GaugeLinkField stapleY(&EdgeGrid);
  GaugeLinkField stapleD(&EdgeGrid);
  GaugeLinkField stapleT(&VertexGrid); // Note live on different Grids

  std::cout << GridLogMessage<<" Computing temporal staples"<<std::endl;

  Support.TemporalStaples(Ut,Umu,stapleX,stapleY,stapleD,stapleT);
  
  std::cout << GridLogMessage<<" Computed temporal staples"<<std::endl;
  linkX = peekLorentz(Umu,IcosahedronPatchX);
  linkY = peekLorentz(Umu,IcosahedronPatchY);
  linkD = peekLorentz(Umu,IcosahedronPatchDiagonal);

  std::cout << GridLogMessage <<"------------------------------------------------"<<std::endl;

  std::cout << GridLogMessage << " trace X*StapleX "<<norm2(trace(linkX * stapleX))<<std::endl;
  std::cout << GridLogMessage << " err " << norm2(trace(linkX * stapleX)-2.0*plaq_ref)<<std::endl;

  std::cout << GridLogMessage << " trace Y*StapleY "<<norm2(trace(linkY * stapleY))<<std::endl;
  std::cout << GridLogMessage << " err " << norm2(trace(linkY * stapleY)-2.0*plaq_ref)<<std::endl;

  std::cout << GridLogMessage << " trace D*StapleD "<<norm2(trace(linkD * stapleD))<<std::endl;
  std::cout << GridLogMessage << " err " << norm2(trace(linkD * stapleD)-2.0*plaq_ref)<<std::endl;

  std::cout << GridLogMessage <<"------------------------------------------------"<<std::endl;

  int L = latt_size[0];
  int T = latt_size[2];

  double T_staple_norm = 10*(L*L-1)*T*36.0 + 12*T*25.0;

  double tr = norm2(trace(Ut * stapleT));
  std::cout << GridLogMessage << " trace T*StapleT "<<tr<<std::endl;
  std::cout << GridLogMessage << " got " << tr<<" expect "<<T_staple_norm<<" diff "<<tr-T_staple_norm <<std::endl;

  
  /*
  GaugeLinkField stapleXTp(&EdgeGrid);
  GaugeLinkField stapleYTp(&EdgeGrid);
  GaugeLinkField stapleDTp(&EdgeGrid);

  GaugeLinkField stapleXTm(&EdgeGrid);
  GaugeLinkField stapleYTm(&EdgeGrid);
  GaugeLinkField stapleDTm(&EdgeGrid);

  GaugeLinkField stapleTXp(&EdgeGrid);
  GaugeLinkField stapleTYp(&EdgeGrid);
  GaugeLinkField stapleTDp(&EdgeGrid);
  GaugeLinkField stapleTXm(&EdgeGrid);
  GaugeLinkField stapleTYm(&EdgeGrid);
  GaugeLinkField stapleTDm(&EdgeGrid);

  std::cout << GridLogMessage<<" Computing temporal staples"<<std::endl;

  Support.TemporalStaplesOld(Ut,Umu,
			  stapleXTp,stapleYTp,stapleDTp,
			  stapleXTm,stapleYTm,stapleDTm,
			  stapleTXp,stapleTYp,stapleTDp,
			  stapleTXm,stapleTYm,stapleTDm);

  std::cout << GridLogMessage<<" Computed temporal staples"<<std::endl;
  linkX = peekLorentz(Umu,IcosahedronPatchX);
  linkY = peekLorentz(Umu,IcosahedronPatchY);
  linkD = peekLorentz(Umu,IcosahedronPatchDiagonal);
  
  Support.CopyNonPoles(Ut,Utnp);
  
  std::cout << GridLogMessage <<"------------------------------------------------"<<std::endl;
  
  std::cout << GridLogMessage << " trace X*StapleXTp "<<norm2(trace(linkX * stapleXTp))<<std::endl;
  std::cout << GridLogMessage << " err " << norm2(trace(linkX * stapleXTp)-plaq_ref)<<std::endl;

  std::cout << GridLogMessage << " trace Y*StapleYTp "<<norm2(trace(linkY * stapleYTp))<<std::endl;
  std::cout << GridLogMessage << " err " << norm2(trace(linkY * stapleYTp)-plaq_ref)<<std::endl;

  std::cout << GridLogMessage << " trace D*StapleDTp "<<norm2(trace(linkD * stapleDTp))<<std::endl;
  std::cout << GridLogMessage << " norm  StapleDTp "<<norm2(stapleDTp)<<std::endl;
  std::cout << GridLogMessage << " err " << norm2(trace(linkD * stapleDTp)-plaq_ref)<<std::endl;

  std::cout << GridLogMessage <<"------------------------------------------------"<<std::endl;

  std::cout << GridLogMessage << " trace T*StapleTXp "<<norm2(trace(Utnp * stapleTXp))<<std::endl;
  std::cout << GridLogMessage << " err " << norm2(trace(Utnp * stapleTXp)-plaq_ref)<<std::endl;

  std::cout << GridLogMessage << " trace T*StapleTYp "<<norm2(trace(Utnp * stapleTYp))<<std::endl;
  std::cout << GridLogMessage << " err " << norm2(trace(Utnp * stapleTYp)-plaq_ref)<<std::endl;

  std::cout << GridLogMessage << " trace T*StapleTDp "<<norm2(trace(Utnp * stapleTDp))<<std::endl;
  std::cout << GridLogMessage << " err " << norm2(trace(Utnp * stapleTDp)-plaq_ref)<<std::endl;

  std::cout << GridLogMessage <<"------------------------------------------------"<<std::endl;
  
  std::cout << GridLogMessage << " trace X*StapleXTm "<<norm2(trace(linkX * stapleXTm))<<std::endl;
  std::cout << GridLogMessage << " err " << norm2(trace(linkX * stapleXTm)-plaq_ref)<<std::endl;

  std::cout << GridLogMessage << " trace Y*StapleYTm "<<norm2(trace(linkY * stapleYTm))<<std::endl;
  std::cout << GridLogMessage << " err " << norm2(trace(linkY * stapleYTm)-plaq_ref)<<std::endl;

  std::cout << GridLogMessage << " trace D*StapleDTm "<<norm2(trace(linkD * stapleDTm))<<std::endl;
  std::cout << GridLogMessage << " err " << norm2(trace(linkD * stapleDTm)-plaq_ref)<<std::endl;

  std::cout << GridLogMessage <<"------------------------------------------------"<<std::endl;

  std::cout << GridLogMessage << " trace T*StapleTXm "<<norm2(trace(Utnp * stapleTXm))<<std::endl;
  std::cout << GridLogMessage << " err " << norm2(trace(Utnp * stapleTXm)-plaq_ref)<<std::endl;

  std::cout << GridLogMessage << " trace T*StapleTYm "<<norm2(trace(Utnp * stapleTYm))<<std::endl;
  std::cout << GridLogMessage << " err " << norm2(trace(Utnp * stapleTYm)-plaq_ref)<<std::endl;

  std::cout << GridLogMessage << " trace T*StapleTDm "<<norm2(trace(Utnp * stapleTDm))<<std::endl;
  std::cout << GridLogMessage << " diff expected to be number of exceptional points " << norm2(trace(Utnp * stapleTDm)-plaq_ref)<<std::endl;

  std::cout << GridLogMessage <<"------------------------------------------------"<<std::endl;
  //  std::cout << closure(Utnp*stapleTXm) << std::endl;
*/
  
  Grid_finalize();
}

