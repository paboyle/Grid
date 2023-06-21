    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_padded_cell_staple.cc

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
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>

using namespace std;
using namespace Grid;

template <class Gimpl> class WilsonLoopsTest : public Gimpl {
public:
  INHERIT_GIMPL_TYPES(Gimpl);

  typedef typename Gimpl::GaugeLinkField GaugeMat;
  typedef typename Gimpl::GaugeField GaugeLorentz;


  //Original implementation
  static void StapleOrig(GaugeMat &staple, const GaugeLorentz &Umu, int mu,
			 int nu) {

    GridBase *grid = Umu.Grid();

    std::vector<GaugeMat> U(Nd, grid);
    for (int d = 0; d < Nd; d++) {
      U[d] = PeekIndex<LorentzIndex>(Umu, d);
    }
    staple = Zero();

    if (nu != mu) {

      // mu
      // ^
      // |__>  nu

      //    __
      //      |
      //    __|
      //

      //Forward: Out(x) = Link(x)*field(x+mu)
      //Backward: Out(x) = Link^dag(x-mu)*field(x-mu)
      //ShiftStaple: Link(x) = Link(x+mu)

      //tmp1 = U^dag_nu(x-nu)
      //tmp2 = U^dag_mu(x-mu) tmp1(x-mu) = U^dag_mu(x-mu) U^dag_nu(x-nu-mu)
      //tmp3 = U_nu(x) tmp2(x+nu) = U_nu(x)U^dag_mu(x-mu+nu) U^dag_nu(x-mu)
      //tmp4 = tmp(x+mu) = U_nu(x+mu)U^dag_mu(x+nu) U^dag_nu(x)

      staple += Gimpl::ShiftStaple(
				   Gimpl::CovShiftForward(
							  U[nu], nu,
							  Gimpl::CovShiftBackward(
										  U[mu], mu, Gimpl::CovShiftIdentityBackward(U[nu], nu))),
				   mu);

      //  __
      // |
      // |__
      //
      //
    
      //tmp1 = U_mu^dag(x-mu) U_nu(x-mu)
      //tmp2 = U_nu^dag(x-nu) tmp1(x-nu) = U_nu^dag(x-nu) U_mu^dag(x-mu-nu) U_nu(x-mu-nu)
      //tmp3 = tmp2(x+mu) = U_nu^dag(x-nu+mu) U_mu^dag(x-nu) U_nu(x-nu)
      staple += Gimpl::ShiftStaple(
				   Gimpl::CovShiftBackward(U[nu], nu,
							   Gimpl::CovShiftBackward(U[mu], mu, U[nu])),
				   mu);
    }
  }

  static void StaplePadded(GaugeMat &staple, const GaugeLorentz &U, int mu,
			   int nu) {
    if(nu==mu){
      staple = Zero();
      return;
    }
    double peek = 0, construct = 0, exchange = 0, coord = 0, stencil =0, kernel = 0, extract = 0, total = 0;
    
    double tstart = usecond();
    double t=tstart;
    
    PaddedCell Ghost(1, (GridCartesian*)U.Grid());

    construct += usecond() - t;
      
    t=usecond();      
    GaugeMat U_mu = PeekIndex<LorentzIndex>(U, mu);
    GaugeMat U_nu = PeekIndex<LorentzIndex>(U, nu);
    peek += usecond() - t;

    t=usecond();
    CshiftImplGauge<Gimpl> cshift_impl;
    GaugeMat Ug_mu = Ghost.Exchange(U_mu, cshift_impl);
    GaugeMat Ug_nu = Ghost.Exchange(U_nu, cshift_impl);
    exchange += usecond() - t;
    
    GridBase *ggrid = Ug_mu.Grid();

    GaugeMat gStaple(ggrid);

    t=usecond();
    Coordinate shift_0(Nd,0);
    Coordinate shift_mu(Nd,0); shift_mu[mu]=1;
    Coordinate shift_nu(Nd,0); shift_nu[nu]=1;
    Coordinate shift_mnu(Nd,0); shift_mnu[nu]=-1;
    Coordinate shift_mnu_pmu(Nd,0); shift_mnu_pmu[nu]=-1; shift_mnu_pmu[mu]=1;

    std::vector<Coordinate> shifts;

    //U_nu(x+mu)U^dag_mu(x+nu) U^dag_nu(x)
    shifts.push_back(shift_0);
    shifts.push_back(shift_nu);
    shifts.push_back(shift_mu);

    //U_nu^dag(x-nu+mu) U_mu^dag(x-nu) U_nu(x-nu)
    shifts.push_back(shift_mnu);
    shifts.push_back(shift_mnu);
    shifts.push_back(shift_mnu_pmu);
    coord += usecond()-t;

    t=usecond();
    GeneralLocalStencil gStencil(ggrid,shifts);
    stencil += usecond() -t;

    t=usecond();
    {
      autoView( gStaple_v , gStaple, AcceleratorWrite);
      auto gStencil_v = gStencil.View();
      autoView( Ug_mu_v , Ug_mu, AcceleratorRead);
      autoView( Ug_nu_v , Ug_nu, AcceleratorRead);
  
      accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
	  GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
	  auto Udag_nu_x = adj(coalescedReadGeneralPermute(Ug_nu_v[e->_offset], e->_permute, Nd));
	  e = gStencil_v.GetEntry(1,ss);
	  auto Udag_mu_xpnu = adj(coalescedReadGeneralPermute(Ug_mu_v[e->_offset], e->_permute, Nd));
	  e = gStencil_v.GetEntry(2,ss);
	  auto U_nu_xpmu = coalescedReadGeneralPermute(Ug_nu_v[e->_offset], e->_permute, Nd);
      
	  auto stencil_ss = U_nu_xpmu * Udag_mu_xpnu * Udag_nu_x;

	  e = gStencil_v.GetEntry(3,ss);
	  auto U_nu_xmnu = coalescedReadGeneralPermute(Ug_nu_v[e->_offset], e->_permute, Nd);
	  e = gStencil_v.GetEntry(4,ss);
	  auto Udag_mu_xmnu = adj(coalescedReadGeneralPermute(Ug_mu_v[e->_offset], e->_permute, Nd));
	  e = gStencil_v.GetEntry(5,ss);
	  auto Udag_nu_xmnu_pmu = adj(coalescedReadGeneralPermute(Ug_nu_v[e->_offset], e->_permute, Nd));

	  stencil_ss = stencil_ss + Udag_nu_xmnu_pmu * Udag_mu_xmnu * U_nu_xmnu;
      
	  coalescedWrite(gStaple_v[ss],stencil_ss);
	}
	);
    } //ensure views are all closed!
    kernel += usecond() - t;

    t=usecond();
    staple = Ghost.Extract(gStaple);
    extract += usecond()-t;
    
    total += usecond() - tstart;
    std::cout << GridLogMessage << "StaplePadded timings peek:" << peek << " construct:" << construct << " exchange:" << exchange << " coord:" << coord << " stencil:" << stencil << " kernel:" << kernel << " extract:" << extract << " total:" << total << std::endl;
  }

  static void RectStapleOrig(GaugeMat &Stap, const GaugeLorentz &Umu,
			     int mu) {
    GridBase *grid = Umu.Grid();

    std::vector<GaugeMat> U(Nd, grid);
    for (int d = 0; d < Nd; d++) {
      U[d] = PeekIndex<LorentzIndex>(Umu, d);
    }

    Stap = Zero();

    for (int nu = 0; nu < Nd; nu++) {
      if (nu != mu) {
        //           __ ___
        //          |    __ |
        //
	//tmp1 = U_nu^dag(x-nu)
	//tmp2 = U_mu^dag(x-mu)tmp1(x-mu) = U_mu^dag(x-mu) U_nu^dag(x-nu-mu)
	//tmp3 = U_mu^dag(x-mu)tmp2(x-mu) = U_mu^dag(x-mu) U_mu^dag(x-2mu) U_nu^dag(x-nu-2mu)
	//tmp4 = U_nu(x)tmp3(x+nu) = U_nu(x)U_mu^dag(x-mu+nu) U_mu^dag(x-2mu+nu) U_nu^dag(x-2mu)
	//tmp5 = U_mu(x)tmp4(x+mu) = U_mu(x)U_nu(x+mu)U_mu^dag(x+nu) U_mu^dag(x-mu+nu) U_nu^dag(x-mu)
	//tmp6 = tmp5(x+mu) = U_mu(x+mu)U_nu(x+2mu)U_mu^dag(x+nu+mu) U_mu^dag(x+nu) U_nu^dag(x)
	
        Stap += Gimpl::ShiftStaple(
				   Gimpl::CovShiftForward(
							  U[mu], mu,
							  Gimpl::CovShiftForward(
										 U[nu], nu,
										 Gimpl::CovShiftBackward(
													 U[mu], mu,
													 Gimpl::CovShiftBackward(
																 U[mu], mu,
																 Gimpl::CovShiftIdentityBackward(U[nu], nu))))),
				   mu);

        //              __
        //          |__ __ |

	//tmp1 = U^dag_mu(x-mu)U_nu(x-mu)
	//tmp2 = U^dag_mu(x-mu)tmp1(x-mu) = U^dag_mu(x-mu)U^dag_mu(x-2mu)U_nu(x-2mu)
	//tmp3 = U^dag_nu(x-nu)tmp2(x-nu) = U^dag_nu(x-nu)U^dag_mu(x-mu-nu)U^dag_mu(x-2mu-nu)U_nu(x-2mu-nu)
	//tmp4 = U_mu(x)tmp3(x+mu) = U_mu(x)U^dag_nu(x-nu+mu)U^dag_mu(x-nu)U^dag_mu(x-mu-nu)U_nu(x-mu-nu)
	//tmp5 = tmp4(x+mu) = U_mu(x+mu)U^dag_nu(x-nu+2mu)U^dag_mu(x-nu+mu)U^dag_mu(x-nu)U_nu(x-nu)
	
        Stap += Gimpl::ShiftStaple(
				   Gimpl::CovShiftForward(
							  U[mu], mu,
							  Gimpl::CovShiftBackward(
										  U[nu], nu,
										  Gimpl::CovShiftBackward(
													  U[mu], mu, Gimpl::CovShiftBackward(U[mu], mu, U[nu])))),
				   mu);

        //           __
        //          |__ __ |
	//Forward: Out(x) = Link(x)*field(x+mu)
	//Backward: Out(x) = Link^dag(x-mu)*field(x-mu)
	//ShiftStaple: Link(x) = Link(x+mu)

	//tmp1 = U_nu(x)U_mu(x+nu)
	//tmp2 = U^dag_mu(x-mu)tmp1(x-mu) = U^dag_mu(x-mu)U_nu(x-mu)U_mu(x+nu-mu)
	//tmp3 = U^dag_mu(x-mu)tmp2(x-mu) = U^dag_mu(x-mu)U^dag_mu(x-2mu)U_nu(x-2mu)U_mu(x+nu-2mu)
	//tmp4 = U^dag_nu(x-nu)tmp3(x-nu) = U^dag_nu(x-nu)U^dag_mu(x-mu-nu)U^dag_mu(x-2mu-nu)U_nu(x-2mu-nu)U_mu(x-2mu)
	//tmp5 = tmp4(x+mu) = U^dag_nu(x-nu+mu)U^dag_mu(x-nu)U^dag_mu(x-mu-nu)U_nu(x-mu-nu)U_mu(x-mu)
        Stap += Gimpl::ShiftStaple(
				   Gimpl::CovShiftBackward(
							   U[nu], nu,
							   Gimpl::CovShiftBackward(
										   U[mu], mu,
										   Gimpl::CovShiftBackward(
													   U[mu], mu, Gimpl::CovShiftForward(U[nu], nu, U[mu])))),
				   mu);

        //           __ ___
        //          |__    |
	//tmp1 = U_nu^dag(x-nu)U_mu(x-nu)
	//tmp2 = U_mu^dag(x-mu)tmp1(x-mu) = U_mu^dag(x-mu)U_nu^dag(x-mu-nu)U_mu(x-mu-nu)
	//tmp3 = U_mu^dag(x-mu)tmp2(x-mu) = U_mu^dag(x-mu)U_mu^dag(x-2mu)U_nu^dag(x-2mu-nu)U_mu(x-2mu-nu)
	//tmp4 = U_nu(x)tmp3(x+nu) = U_nu(x)U_mu^dag(x-mu+nu)U_mu^dag(x-2mu+nu)U_nu^dag(x-2mu)U_mu(x-2mu)
	//tmp5 = tmp4(x+mu) = U_nu(x+mu)U_mu^dag(x+nu)U_mu^dag(x-mu+nu)U_nu^dag(x-mu)U_mu(x-mu)
        Stap += Gimpl::ShiftStaple(
				   Gimpl::CovShiftForward(
							  U[nu], nu,
							  Gimpl::CovShiftBackward(
										  U[mu], mu,
										  Gimpl::CovShiftBackward(
													  U[mu], mu, Gimpl::CovShiftBackward(U[nu], nu, U[mu])))),
				   mu);

        //       --
        //      |  |
        //
        //      |  |
	//tmp1 = U_nu^dag(x-nu)
	//tmp2 = U_nu^dag(x-nu)tmp1(x-nu) = U_nu^dag(x-nu)U_nu^dag(x-2nu)
	//tmp3 = U_mu^dag(x-mu)tmp2(x-mu) = U_mu^dag(x-mu)U_nu^dag(x-mu-nu)U_nu^dag(x-mu-2nu)
	//tmp4 = U_nu(x)tmp3(x+nu) = U_nu(x)U_mu^dag(x-mu+nu)U_nu^dag(x-mu)U_nu^dag(x-mu-nu)
	//tmp5 = U_nu(x)tmp4(x+nu) = U_nu(x)U_nu(x+nu)U_mu^dag(x-mu+2nu)U_nu^dag(x-mu+nu)U_nu^dag(x-mu)
	//tmp6 = tmp5(x+mu) = U_nu(x+mu)U_nu(x+mu+nu)U_mu^dag(x+2nu)U_nu^dag(x+nu)U_nu^dag(x)
        Stap += Gimpl::ShiftStaple(
				   Gimpl::CovShiftForward(
							  U[nu], nu,
							  Gimpl::CovShiftForward(
										 U[nu], nu,
										 Gimpl::CovShiftBackward(
													 U[mu], mu,
													 Gimpl::CovShiftBackward(
																 U[nu], nu,
																 Gimpl::CovShiftIdentityBackward(U[nu], nu))))),
				   mu);

        //      |  |
        //
        //      |  |
        //       --
	//tmp1 = U_nu(x)U_nu(x+nu)
	//tmp2 = U_mu^dag(x-mu)tmp1(x-mu) = U_mu^dag(x-mu)U_nu(x-mu)U_nu(x-mu+nu)
	//tmp3 = U_nu^dag(x-nu)tmp2(x-nu) = U_nu^dag(x-nu)U_mu^dag(x-mu-nu)U_nu(x-mu-nu)U_nu(x-mu)
	//tmp4 = U_nu^dag(x-nu)tmp3(x-nu) = U_nu^dag(x-nu)U_nu^dag(x-2nu)U_mu^dag(x-mu-2nu)U_nu(x-mu-2nu)U_nu(x-mu-nu)
	//tmp5 = tmp4(x+mu) = U_nu^dag(x+mu-nu)U_nu^dag(x+mu-2nu)U_mu^dag(x-2nu)U_nu(x-2nu)U_nu(x-nu)
        Stap += Gimpl::ShiftStaple(
				   Gimpl::CovShiftBackward(
							   U[nu], nu,
							   Gimpl::CovShiftBackward(
										   U[nu], nu,
										   Gimpl::CovShiftBackward(
													   U[mu], mu, Gimpl::CovShiftForward(U[nu], nu, U[nu])))),
				   mu);
      }
    }
  }


  static void RectStaplePadded(GaugeMat &Stap, const GaugeLorentz &U,
			       int mu) {
    PaddedCell Ghost(2,(GridCartesian*)U.Grid());
    GridBase *ggrid = Ghost.grids.back();
    
    CshiftImplGauge<Gimpl> cshift_impl;
    std::vector<GaugeMat> Ug_dirs(Nd,ggrid);
    for(int i=0;i<Nd;i++) Ug_dirs[i] = Ghost.Exchange(PeekIndex<LorentzIndex>(U, i), cshift_impl);

    GaugeMat gStaple(ggrid);

    std::vector<Coordinate> shifts;
    for (int nu = 0; nu < Nd; nu++) {
      if (nu != mu) {
	auto genShift = [&](int mushift,int nushift){
	  Coordinate out(Nd,0); out[mu]=mushift; out[nu]=nushift; return out;
	};

	//tmp6 = tmp5(x+mu) = U_mu(x+mu)U_nu(x+2mu)U_mu^dag(x+nu+mu) U_mu^dag(x+nu) U_nu^dag(x)
	shifts.push_back(genShift(0,0));
	shifts.push_back(genShift(0,+1));
	shifts.push_back(genShift(+1,+1));
	shifts.push_back(genShift(+2,0));
	shifts.push_back(genShift(+1,0));

	//tmp5 = tmp4(x+mu) = U_mu(x+mu)U^dag_nu(x-nu+2mu)U^dag_mu(x-nu+mu)U^dag_mu(x-nu)U_nu(x-nu)
	shifts.push_back(genShift(0,-1));
	shifts.push_back(genShift(0,-1));
	shifts.push_back(genShift(+1,-1));
	shifts.push_back(genShift(+2,-1));
	shifts.push_back(genShift(+1,0));

	//tmp5 = tmp4(x+mu) = U^dag_nu(x-nu+mu)U^dag_mu(x-nu)U^dag_mu(x-mu-nu)U_nu(x-mu-nu)U_mu(x-mu)
	shifts.push_back(genShift(-1,0));
	shifts.push_back(genShift(-1,-1));
	shifts.push_back(genShift(-1,-1));
	shifts.push_back(genShift(0,-1));
	shifts.push_back(genShift(+1,-1));

	//tmp5 = tmp4(x+mu) = U_nu(x+mu)U_mu^dag(x+nu)U_mu^dag(x-mu+nu)U_nu^dag(x-mu)U_mu(x-mu)
	shifts.push_back(genShift(-1,0));
	shifts.push_back(genShift(-1,0));
	shifts.push_back(genShift(-1,+1));
	shifts.push_back(genShift(0,+1));
	shifts.push_back(genShift(+1,0));

	//tmp6 = tmp5(x+mu) = U_nu(x+mu)U_nu(x+mu+nu)U_mu^dag(x+2nu)U_nu^dag(x+nu)U_nu^dag(x)
	shifts.push_back(genShift(0,0));
	shifts.push_back(genShift(0,+1));
	shifts.push_back(genShift(0,+2));
	shifts.push_back(genShift(+1,+1));
	shifts.push_back(genShift(+1,0));

	//tmp5 = tmp4(x+mu) = U_nu^dag(x+mu-nu)U_nu^dag(x+mu-2nu)U_mu^dag(x-2nu)U_nu(x-2nu)U_nu(x-nu)
	shifts.push_back(genShift(0,-1));
	shifts.push_back(genShift(0,-2));
	shifts.push_back(genShift(0,-2));
	shifts.push_back(genShift(+1,-2));
	shifts.push_back(genShift(+1,-1));
      }
    }
    size_t nshift = shifts.size();

    GeneralLocalStencil gStencil(ggrid,shifts);
    {
      autoView( gStaple_v , gStaple, AcceleratorWrite);
      auto gStencil_v = gStencil.View();

      typedef LatticeView<typename GaugeMat::vector_object> GaugeViewType;
      size_t vsize = Nd*sizeof(GaugeViewType);
      GaugeViewType* Ug_dirs_v_host = (GaugeViewType*)malloc(vsize);
      for(int i=0;i<Nd;i++) Ug_dirs_v_host[i] = Ug_dirs[i].View(AcceleratorRead);
      GaugeViewType* Ug_dirs_v = (GaugeViewType*)acceleratorAllocDevice(vsize);
      acceleratorCopyToDevice(Ug_dirs_v_host,Ug_dirs_v,vsize);

      accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
	  decltype(coalescedRead(Ug_dirs_v[0][0])) stencil_ss;
	  stencil_ss = Zero();
	  int s=0;
	  for(int nu=0;nu<Nd;nu++){
	    if(nu != mu){
	      //tmp6 = tmp5(x+mu) = U_mu(x+mu)U_nu(x+2mu)U_mu^dag(x+nu+mu) U_mu^dag(x+nu) U_nu^dag(x)
	      GeneralStencilEntry const* e = gStencil_v.GetEntry(s++,ss);
	      auto U0 = adj(coalescedReadGeneralPermute(Ug_dirs_v[nu][e->_offset], e->_permute, Nd));
	      e = gStencil_v.GetEntry(s++,ss);
	      auto U1 = adj(coalescedReadGeneralPermute(Ug_dirs_v[mu][e->_offset], e->_permute, Nd));
	      e = gStencil_v.GetEntry(s++,ss);
	      auto U2 = adj(coalescedReadGeneralPermute(Ug_dirs_v[mu][e->_offset], e->_permute, Nd));
	      e = gStencil_v.GetEntry(s++,ss);
	      auto U3 = coalescedReadGeneralPermute(Ug_dirs_v[nu][e->_offset], e->_permute, Nd);
	      e = gStencil_v.GetEntry(s++,ss);
	      auto U4 = coalescedReadGeneralPermute(Ug_dirs_v[mu][e->_offset], e->_permute, Nd);
	    
	      stencil_ss = stencil_ss + U4*U3*U2*U1*U0;

	      //tmp5 = tmp4(x+mu) = U_mu(x+mu)U^dag_nu(x-nu+2mu)U^dag_mu(x-nu+mu)U^dag_mu(x-nu)U_nu(x-nu)
	      e = gStencil_v.GetEntry(s++,ss);
	      U0 = coalescedReadGeneralPermute(Ug_dirs_v[nu][e->_offset], e->_permute, Nd);
	      e = gStencil_v.GetEntry(s++,ss);
	      U1 = adj(coalescedReadGeneralPermute(Ug_dirs_v[mu][e->_offset], e->_permute, Nd));
	      e = gStencil_v.GetEntry(s++,ss);
	      U2 = adj(coalescedReadGeneralPermute(Ug_dirs_v[mu][e->_offset], e->_permute, Nd));
	      e = gStencil_v.GetEntry(s++,ss);
	      U3 = adj(coalescedReadGeneralPermute(Ug_dirs_v[nu][e->_offset], e->_permute, Nd));
	      e = gStencil_v.GetEntry(s++,ss);
	      U4 = coalescedReadGeneralPermute(Ug_dirs_v[mu][e->_offset], e->_permute, Nd);

	      stencil_ss = stencil_ss + U4*U3*U2*U1*U0;

	      //tmp5 = tmp4(x+mu) = U^dag_nu(x-nu+mu)U^dag_mu(x-nu)U^dag_mu(x-mu-nu)U_nu(x-mu-nu)U_mu(x-mu)
	      e = gStencil_v.GetEntry(s++,ss);
	      U0 = coalescedReadGeneralPermute(Ug_dirs_v[mu][e->_offset], e->_permute, Nd);
	      e = gStencil_v.GetEntry(s++,ss);
	      U1 = coalescedReadGeneralPermute(Ug_dirs_v[nu][e->_offset], e->_permute, Nd);
	      e = gStencil_v.GetEntry(s++,ss);
	      U2 = adj(coalescedReadGeneralPermute(Ug_dirs_v[mu][e->_offset], e->_permute, Nd));
	      e = gStencil_v.GetEntry(s++,ss);
	      U3 = adj(coalescedReadGeneralPermute(Ug_dirs_v[mu][e->_offset], e->_permute, Nd));
	      e = gStencil_v.GetEntry(s++,ss);
	      U4 = adj(coalescedReadGeneralPermute(Ug_dirs_v[nu][e->_offset], e->_permute, Nd));

	      stencil_ss = stencil_ss + U4*U3*U2*U1*U0;

	      //tmp5 = tmp4(x+mu) = U_nu(x+mu)U_mu^dag(x+nu)U_mu^dag(x-mu+nu)U_nu^dag(x-mu)U_mu(x-mu)
	      e = gStencil_v.GetEntry(s++,ss);
	      U0 = coalescedReadGeneralPermute(Ug_dirs_v[mu][e->_offset], e->_permute, Nd);
	      e = gStencil_v.GetEntry(s++,ss);
	      U1 = adj(coalescedReadGeneralPermute(Ug_dirs_v[nu][e->_offset], e->_permute, Nd));
	      e = gStencil_v.GetEntry(s++,ss);
	      U2 = adj(coalescedReadGeneralPermute(Ug_dirs_v[mu][e->_offset], e->_permute, Nd));
	      e = gStencil_v.GetEntry(s++,ss);
	      U3 = adj(coalescedReadGeneralPermute(Ug_dirs_v[mu][e->_offset], e->_permute, Nd));
	      e = gStencil_v.GetEntry(s++,ss);
	      U4 = coalescedReadGeneralPermute(Ug_dirs_v[nu][e->_offset], e->_permute, Nd);

	      stencil_ss = stencil_ss + U4*U3*U2*U1*U0;

	      //tmp6 = tmp5(x+mu) = U_nu(x+mu)U_nu(x+mu+nu)U_mu^dag(x+2nu)U_nu^dag(x+nu)U_nu^dag(x)
	      e = gStencil_v.GetEntry(s++,ss);
	      U0 = adj(coalescedReadGeneralPermute(Ug_dirs_v[nu][e->_offset], e->_permute, Nd));
	      e = gStencil_v.GetEntry(s++,ss);
	      U1 = adj(coalescedReadGeneralPermute(Ug_dirs_v[nu][e->_offset], e->_permute, Nd));
	      e = gStencil_v.GetEntry(s++,ss);
	      U2 = adj(coalescedReadGeneralPermute(Ug_dirs_v[mu][e->_offset], e->_permute, Nd));
	      e = gStencil_v.GetEntry(s++,ss);
	      U3 = coalescedReadGeneralPermute(Ug_dirs_v[nu][e->_offset], e->_permute, Nd);
	      e = gStencil_v.GetEntry(s++,ss);
	      U4 = coalescedReadGeneralPermute(Ug_dirs_v[nu][e->_offset], e->_permute, Nd);

	      stencil_ss = stencil_ss + U4*U3*U2*U1*U0;   

	      //tmp5 = tmp4(x+mu) = U_nu^dag(x+mu-nu)U_nu^dag(x+mu-2nu)U_mu^dag(x-2nu)U_nu(x-2nu)U_nu(x-nu)
	      e = gStencil_v.GetEntry(s++,ss);
	      U0 = coalescedReadGeneralPermute(Ug_dirs_v[nu][e->_offset], e->_permute, Nd);
	      e = gStencil_v.GetEntry(s++,ss);
	      U1 = coalescedReadGeneralPermute(Ug_dirs_v[nu][e->_offset], e->_permute, Nd);
	      e = gStencil_v.GetEntry(s++,ss);
	      U2 = adj(coalescedReadGeneralPermute(Ug_dirs_v[mu][e->_offset], e->_permute, Nd));
	      e = gStencil_v.GetEntry(s++,ss);
	      U3 = adj(coalescedReadGeneralPermute(Ug_dirs_v[nu][e->_offset], e->_permute, Nd));
	      e = gStencil_v.GetEntry(s++,ss);
	      U4 = adj(coalescedReadGeneralPermute(Ug_dirs_v[nu][e->_offset], e->_permute, Nd));

	      stencil_ss = stencil_ss + U4*U3*U2*U1*U0;   

	    }
	  }
	  assert(s==nshift);
	  coalescedWrite(gStaple_v[ss],stencil_ss);
	}
	);
  
      for(int i=0;i<Nd;i++) Ug_dirs_v_host[i].ViewClose();
      free(Ug_dirs_v_host);
      acceleratorFreeDevice(Ug_dirs_v);
    }   
    Stap = Ghost.Extract(gStaple);    
  }



};  
  
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
  LatticeGaugeField U(&GRID);

  SU<Nc>::HotConfiguration(pRNG,U);

  //typedef PeriodicGimplD Gimpl;
  typedef ConjugateGimplD Gimpl;
  std::vector<int> conj_dirs(Nd,0); conj_dirs[0]=1; conj_dirs[3]=1;
  Gimpl::setDirections(conj_dirs);

  typedef typename WilsonLoopsTest<Gimpl>::GaugeMat GaugeMat;
  typedef typename WilsonLoopsTest<Gimpl>::GaugeLorentz GaugeLorentz;

  std::cout << GridLogMessage << "Checking Staple" << std::endl;
  int count = 0;
  double torig=0, tpadded=0;
  
  for(int mu=0;mu<Nd;mu++){
    for(int nu=0;nu<Nd;nu++){
      if(mu != nu){
	GaugeMat staple_orig(&GRID), staple_padded(&GRID);
	double t0 = usecond();
	WilsonLoopsTest<Gimpl>::StapleOrig(staple_orig,U,mu,nu);
	double t1 = usecond();
	WilsonLoopsTest<Gimpl>::StaplePadded(staple_padded,U,mu,nu);
	double t2 = usecond();
	torig += t1-t0;  tpadded += t2-t1;
	++count;
	
	GaugeMat diff = staple_orig - staple_padded;
	double n = norm2(diff);
	std::cout << GridLogMessage << mu << " " << nu << " " << n << std::endl;
	assert(n<1e-10);
      }
    }
  }
  std::cout << GridLogMessage << "Staple timings orig: " << torig/1000/count << "ms,  padded: " << tpadded/1000/count << "ms" << std::endl;
  count=0; torig=tpadded=0;
    
  std::cout << GridLogMessage << "Checking RectStaple" << std::endl;
  for(int mu=0;mu<Nd;mu++){
    GaugeMat staple_orig(&GRID), staple_padded(&GRID);
    double t0 = usecond();
    WilsonLoopsTest<Gimpl>::RectStapleOrig(staple_orig,U,mu);
    double t1 = usecond();
    WilsonLoopsTest<Gimpl>::RectStaplePadded(staple_padded,U,mu);
    double t2 = usecond();
    torig += t1-t0;  tpadded += t2-t1;
    ++count;
    
    GaugeMat diff = staple_orig - staple_padded;
    double n = norm2(diff);
    std::cout << GridLogMessage << mu << " " << n << std::endl;
    assert(n<1e-10);
  }
  std::cout << GridLogMessage << "RectStaple timings orig: " << torig/1000/count << "ms,  padded: " << tpadded/1000/count << "ms" << std::endl;
  
  Grid_finalize();
}
