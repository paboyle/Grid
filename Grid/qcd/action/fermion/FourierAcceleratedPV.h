
    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/FourierAcceleratedPV.h

    Copyright (C) 2015

Author: Christoph Lehner (lifted with permission by Peter Boyle, brought back to Grid)
Author: Peter Boyle <pabobyle@ph.ed.ac.uk>

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
#pragma once
namespace Grid {
namespace QCD {

  template<typename M>
    void get_real_const_bc(M& m, RealD& _b, RealD& _c) {
    ComplexD b,c;
    b=m.bs[0];
    c=m.cs[0];
    std::cout << GridLogMessage << "b=" << b << ", c=" << c << std::endl;
    for (size_t i=1;i<m.bs.size();i++) {
      assert(m.bs[i] == b);
      assert(m.cs[i] == c);
    }
    assert(b.imag() == 0.0);
    assert(c.imag() == 0.0);
    _b = b.real();
    _c = c.real();
  }


template<typename Vi, typename M, typename G>
class FourierAcceleratedPV {
 public:

  ConjugateGradient<Vi> &cg;
  M& dwfPV;
  G& Umu;
  GridCartesian* grid5D;
  GridRedBlackCartesian* gridRB5D;
  int group_in_s;

  FourierAcceleratedPV(M& _dwfPV, G& _Umu, ConjugateGradient<Vi> &_cg, int _group_in_s = 2) 
   : dwfPV(_dwfPV), Umu(_Umu), cg(_cg), group_in_s(_group_in_s) 
  {
    assert( dwfPV.FermionGrid()->_fdimensions[0] % (2*group_in_s) == 0);
    grid5D = QCD::SpaceTimeGrid::makeFiveDimGrid(2*group_in_s, (GridCartesian*)Umu._grid);
    gridRB5D = QCD::SpaceTimeGrid::makeFiveDimRedBlackGrid(2*group_in_s, (GridCartesian*)Umu._grid);
  }

  void rotatePV(const Vi& _src, Vi& dst, bool forward) const {

    GridStopWatch gsw1, gsw2;

    typedef typename Vi::scalar_type Coeff_t;
    int Ls = dst._grid->_fdimensions[0];

    Vi _tmp(dst._grid);
    double phase = M_PI / (double)Ls;
    Coeff_t bzero(0.0,0.0);

    FFT theFFT((GridCartesian*)dst._grid);

    if (!forward) {
      gsw1.Start();
      for (int s=0;s<Ls;s++) {
	Coeff_t a(::cos(phase*s),-::sin(phase*s));
	axpby_ssp(_tmp,a,_src,bzero,_src,s,s);
      }
      gsw1.Stop();

      gsw2.Start();
      theFFT.FFT_dim(dst,_tmp,0,FFT::forward);
      gsw2.Stop();

    } else {

      gsw2.Start();
      theFFT.FFT_dim(_tmp,_src,0,FFT::backward);
      gsw2.Stop();

      gsw1.Start();
      for (int s=0;s<Ls;s++) {
	Coeff_t a(::cos(phase*s),::sin(phase*s));
	axpby_ssp(dst,a,_tmp,bzero,_tmp,s,s);
      }
      gsw1.Stop();
    }

    std::cout << GridLogMessage << "Timing rotatePV: " << gsw1.Elapsed() << ", " << gsw2.Elapsed() << std::endl;

  }

  void pvInv(const Vi& _src, Vi& _dst) const {

    std::cout << GridLogMessage << "Fourier-Accelerated Outer Pauli Villars"<<std::endl;

    typedef typename Vi::scalar_type Coeff_t;
    int Ls = _dst._grid->_fdimensions[0];

    GridStopWatch gswT;
    gswT.Start();

    RealD b,c;
    get_real_const_bc(dwfPV,b,c);
    RealD M5 = dwfPV.M5;
    
    // U(true) Rightinv TMinv U(false) = Minv

    Vi _src_diag(_dst._grid);
    Vi _src_diag_slice(dwfPV.GaugeGrid());
    Vi _dst_diag_slice(dwfPV.GaugeGrid());
    Vi _src_diag_slices(grid5D);
    Vi _dst_diag_slices(grid5D);
    Vi _dst_diag(_dst._grid);

    rotatePV(_src,_src_diag,false);

    // now do TM solves
    Gamma G5(Gamma::Algebra::Gamma5);

    GridStopWatch gswA, gswB;

    gswA.Start();

    typedef typename M::Impl_t Impl;
    //WilsonTMFermion<Impl> tm(x.Umu,*x.UGridF,*x.UrbGridF,0.0,0.0,solver_outer.parent.par.wparams_f);
    std::vector<RealD> vmass(grid5D->_fdimensions[0],0.0);
    std::vector<RealD> vmu(grid5D->_fdimensions[0],0.0);

    WilsonTMFermion5D<Impl> tm(Umu,*grid5D,*gridRB5D,
			   *(GridCartesian*)dwfPV.GaugeGrid(),
			   *(GridRedBlackCartesian*)dwfPV.GaugeRedBlackGrid(),
			   vmass,vmu);
    
    //SchurRedBlackDiagTwoSolve<Vi> sol(cg);
    SchurRedBlackDiagMooeeSolve<Vi> sol(cg); // same performance as DiagTwo
    gswA.Stop();

    gswB.Start();

    for (int sgroup=0;sgroup<Ls/2/group_in_s;sgroup++) {

      for (int sidx=0;sidx<group_in_s;sidx++) {

	int s = sgroup*group_in_s + sidx;
	int sprime = Ls-s-1;

	RealD phase = M_PI / (RealD)Ls * (2.0 * s + 1.0);
	RealD cosp = ::cos(phase);
	RealD sinp = ::sin(phase);
	RealD denom = b*b + c*c + 2.0*b*c*cosp;
	RealD mass = -(b*b*M5 + c*(1.0 - cosp + c*M5) + b*(-1.0 + cosp + 2.0*c*cosp*M5))/denom;
	RealD mu = (b+c)*sinp/denom;

	vmass[2*sidx + 0] = mass;
	vmass[2*sidx + 1] = mass;
	vmu[2*sidx + 0] = mu;
	vmu[2*sidx + 1] = -mu;

      }

      tm.update(vmass,vmu);

      for (int sidx=0;sidx<group_in_s;sidx++) {

	int s = sgroup*group_in_s + sidx;
	int sprime = Ls-s-1;

	ExtractSlice(_src_diag_slice,_src_diag,s,0);
	InsertSlice(_src_diag_slice,_src_diag_slices,2*sidx + 0,0);

	ExtractSlice(_src_diag_slice,_src_diag,sprime,0);
	InsertSlice(_src_diag_slice,_src_diag_slices,2*sidx + 1,0);

      }

      GridStopWatch gsw;
      gsw.Start();
      _dst_diag_slices = zero; // zero guess
      sol(tm,_src_diag_slices,_dst_diag_slices);
      gsw.Stop();
      std::cout << GridLogMessage << "Solve[sgroup=" << sgroup << "] completed in " << gsw.Elapsed() << ", " << gswA.Elapsed() << std::endl;

      for (int sidx=0;sidx<group_in_s;sidx++) {

	int s = sgroup*group_in_s + sidx;
	int sprime = Ls-s-1;

	RealD phase = M_PI / (RealD)Ls * (2.0 * s + 1.0);
	RealD cosp = ::cos(phase);
	RealD sinp = ::sin(phase);

	// now rotate with inverse of
	Coeff_t pA = b + c*cosp;
	Coeff_t pB = - Coeff_t(0.0,1.0)*c*sinp;
	Coeff_t pABden = pA*pA - pB*pB;
	// (pA + pB * G5) * (pA - pB*G5) = (pA^2 - pB^2)
      
	ExtractSlice(_dst_diag_slice,_dst_diag_slices,2*sidx + 0,0);
	_dst_diag_slice = (pA/pABden) * _dst_diag_slice - (pB/pABden) * (G5 * _dst_diag_slice);
	InsertSlice(_dst_diag_slice,_dst_diag,s,0);
	
	ExtractSlice(_dst_diag_slice,_dst_diag_slices,2*sidx + 1,0);
	_dst_diag_slice = (pA/pABden) * _dst_diag_slice + (pB/pABden) * (G5 * _dst_diag_slice);
	InsertSlice(_dst_diag_slice,_dst_diag,sprime,0);
      }
    }
    gswB.Stop();

    rotatePV(_dst_diag,_dst,true);

    gswT.Stop();
    std::cout << GridLogMessage << "PV completed in " << gswT.Elapsed() << " (Setup: " << gswA.Elapsed() << ", s-loop: " << gswB.Elapsed() << ")" << std::endl;
  }

};
}}
