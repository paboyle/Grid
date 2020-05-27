/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/action/fermion/WilsonCloverFermion.cc

    Copyright (C) 2017

    Author: paboyle <paboyle@ph.ed.ac.uk>
    Author: Guido Cossu <guido.cossu@ed.ac.uk>

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
#include <Grid/qcd/spin/Dirac.h>
#include <Grid/qcd/action/fermion/WilsonCloverFermion.h>

NAMESPACE_BEGIN(Grid);

// *NOT* EO
template <class Impl>
void WilsonCloverFermion<Impl>::M(const FermionField &in, FermionField &out)
{
  FermionField temp(out.Grid());

  // Wilson term
  out.Checkerboard() = in.Checkerboard();
  this->Dhop(in, out, DaggerNo);

  // Clover term
  Mooee(in, temp);

  out += temp;
}

template <class Impl>
void WilsonCloverFermion<Impl>::Mdag(const FermionField &in, FermionField &out)
{
  FermionField temp(out.Grid());

  // Wilson term
  out.Checkerboard() = in.Checkerboard();
  this->Dhop(in, out, DaggerYes);

  // Clover term
  MooeeDag(in, temp);

  out += temp;
}

template <class Impl>
void WilsonCloverFermion<Impl>::ImportGauge(const GaugeField &_Umu)
{
  WilsonFermion<Impl>::ImportGauge(_Umu);
  GridBase *grid = _Umu.Grid();
  typename Impl::GaugeLinkField Bx(grid), By(grid), Bz(grid), Ex(grid), Ey(grid), Ez(grid);

  // Compute the field strength terms mu>nu
  WilsonLoops<Impl>::FieldStrength(Bx, _Umu, Zdir, Ydir);
  WilsonLoops<Impl>::FieldStrength(By, _Umu, Zdir, Xdir);
  WilsonLoops<Impl>::FieldStrength(Bz, _Umu, Ydir, Xdir);
  WilsonLoops<Impl>::FieldStrength(Ex, _Umu, Tdir, Xdir);
  WilsonLoops<Impl>::FieldStrength(Ey, _Umu, Tdir, Ydir);
  WilsonLoops<Impl>::FieldStrength(Ez, _Umu, Tdir, Zdir);

  // Compute the Clover Operator acting on Colour and Spin
  // multiply here by the clover coefficients for the anisotropy
  CloverTerm  = fillCloverYZ(Bx) * csw_r;
  CloverTerm += fillCloverXZ(By) * csw_r;
  CloverTerm += fillCloverXY(Bz) * csw_r;
  CloverTerm += fillCloverXT(Ex) * csw_t;
  CloverTerm += fillCloverYT(Ey) * csw_t;
  CloverTerm += fillCloverZT(Ez) * csw_t;
  CloverTerm += diag_mass;

  int lvol = _Umu.Grid()->lSites();
  int DimRep = Impl::Dimension;

  Eigen::MatrixXcd EigenCloverOp = Eigen::MatrixXcd::Zero(Ns * DimRep, Ns * DimRep);
  Eigen::MatrixXcd EigenInvCloverOp = Eigen::MatrixXcd::Zero(Ns * DimRep, Ns * DimRep);

  Coordinate lcoor;
  typename SiteCloverType::scalar_object Qx = Zero(), Qxinv = Zero();

  for (int site = 0; site < lvol; site++)
  {
    grid->LocalIndexToLocalCoor(site, lcoor);
    EigenCloverOp = Eigen::MatrixXcd::Zero(Ns * DimRep, Ns * DimRep);
    peekLocalSite(Qx, CloverTerm, lcoor);
    Qxinv = Zero();
    //if (csw!=0){
    for (int j = 0; j < Ns; j++)
      for (int k = 0; k < Ns; k++)
        for (int a = 0; a < DimRep; a++)
          for (int b = 0; b < DimRep; b++){
	    auto zz =  Qx()(j, k)(a, b);
            EigenCloverOp(a + j * DimRep, b + k * DimRep) = std::complex<double>(zz);
	  }
    //   if (site==0) std::cout << "site =" << site << "\n" << EigenCloverOp << std::endl;

    EigenInvCloverOp = EigenCloverOp.inverse();
    //std::cout << EigenInvCloverOp << std::endl;
    for (int j = 0; j < Ns; j++)
      for (int k = 0; k < Ns; k++)
        for (int a = 0; a < DimRep; a++)
          for (int b = 0; b < DimRep; b++)
            Qxinv()(j, k)(a, b) = EigenInvCloverOp(a + j * DimRep, b + k * DimRep);
    //    if (site==0) std::cout << "site =" << site << "\n" << EigenInvCloverOp << std::endl;
    //  }
    pokeLocalSite(Qxinv, CloverTermInv, lcoor);
  }

  // Separate the even and odd parts
  pickCheckerboard(Even, CloverTermEven, CloverTerm);
  pickCheckerboard(Odd, CloverTermOdd, CloverTerm);

  pickCheckerboard(Even, CloverTermInvEven, CloverTermInv);
  pickCheckerboard(Odd, CloverTermInvOdd, CloverTermInv);
}

template <class Impl>
void WilsonCloverFermion<Impl>::Mooee(const FermionField &in, FermionField &out)
{
  this->MooeeInternal(in, out, DaggerNo, InverseNo);
}

template <class Impl>
void WilsonCloverFermion<Impl>::MooeeDag(const FermionField &in, FermionField &out)
{
  this->MooeeInternal(in, out, DaggerYes, InverseNo);
}

template <class Impl>
void WilsonCloverFermion<Impl>::MooeeInv(const FermionField &in, FermionField &out)
{
  this->MooeeInternal(in, out, DaggerNo, InverseYes);
}

template <class Impl>
void WilsonCloverFermion<Impl>::MooeeInvDag(const FermionField &in, FermionField &out)
{
  this->MooeeInternal(in, out, DaggerYes, InverseYes);
}

template <class Impl>
void WilsonCloverFermion<Impl>::MooeeInternal(const FermionField &in, FermionField &out, int dag, int inv)
{
  out.Checkerboard() = in.Checkerboard();
  assert(in.Checkerboard() == Odd || in.Checkerboard() == Even);

  const CloverFieldType *Clover = GetCompatibleCloverField(in, inv);
  assert(Clover != nullptr);

  if(dag == DaggerYes)
    this->MultClovDagInternal(*Clover,1,Clover->oSites(),in,out);
  else
    this->MultClovInternal(*Clover,1,Clover->oSites(),in,out);
}

template <class Impl>
const typename WilsonCloverFermion<Impl>::CloverFieldType* WilsonCloverFermion<Impl>::GetCompatibleCloverField(const FermionField &in, int inv)
{
  if (in.Grid()->_isCheckerBoarded) {
    if (in.Checkerboard() == Odd) {
      return (inv == InverseYes) ? &CloverTermInvOdd : &CloverTermOdd;
    }
    else {
      return (inv == InverseYes) ? &CloverTermInvEven : &CloverTermEven;
    }
  }
  else {
    return (inv == InverseYes) ? &CloverTermInv : &CloverTerm;
  }
}

// the following 2 functions should move to a clover kernels class
template <class Impl>
void WilsonCloverFermion<Impl>::MultClovInternal(const CloverFieldType &clov, int Ls, int Nsite, const FermionField &in, FermionField &out)
{
  auto clov_v = clov.View();
  auto in_v   = in.View();
  auto out_v  = out.View();

#ifdef GRID_NVCC
  const uint64_t NN = Nsite*Ls;
  accelerator_forNB( ss, NN, Simd::Nsimd(), {
    int sF = ss;
    int sU = ss/Ls;
    coalescedWrite(out_v[sF], clov_v(sU) * in_v(sF));
  });
#else
  thread_for( ss, Nsite, {
    int sU = ss;
    for(int s=0;s<Ls;s++) {
      int sF = ss*Ls+s;
      vstream(out_v[sF], clov_v(sU) * in_v(sF));
    }
  });
#endif
}

template <class Impl>
void WilsonCloverFermion<Impl>::MultClovDagInternal(const CloverFieldType &clov, int Ls, int Nsite, const FermionField &in, FermionField &out)
{
  auto clov_v = clov.View();
  auto in_v   = in.View();
  auto out_v  = out.View();

#ifdef GRID_NVCC
  const uint64_t NN = Nsite*Ls;
  accelerator_forNB( ss, NN, Simd::Nsimd(), {
    int sF = ss;
    int sU = ss/Ls;
    coalescedWrite(out_v[sF], adj(clov_v(sU)) * in_v(sF));
  });
#else
  thread_for( ss, Nsite, {
    int sU = ss;
    for(int s=0;s<Ls;s++) {
      int sF = ss*Ls+s;
      vstream(out_v[sF], adj(clov_v(sU)) * in_v(sF));
    }
  });
#endif
}

// Derivative parts
template <class Impl>
void WilsonCloverFermion<Impl>::MooDeriv(GaugeField &mat, const FermionField &X, const FermionField &Y, int dag)
{
  assert(0);
}

// Derivative parts
template <class Impl>
void WilsonCloverFermion<Impl>::MeeDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag)
{
  assert(0); // not implemented yet
}

NAMESPACE_END(Grid);
