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
#include <Grid/Eigen/Dense>
#include <Grid/qcd/spin/Dirac.h>

namespace Grid
{
namespace QCD
{

// *NOT* EO
template <class Impl>
RealD WilsonCloverFermion<Impl>::M(const FermionField &in, FermionField &out)
{
  // Wilson term
  out.checkerboard = in.checkerboard;
  this->Dhop(in, out, DaggerNo);
  // Clover term
  // apply the sigma and Fmunu
  FermionField temp(out._grid);
  Mooee(in, temp);
  // overall factor
  out += temp;
  return axpy_norm(out, 4 + this->mass, in, out);
}

template <class Impl>
RealD WilsonCloverFermion<Impl>::Mdag(const FermionField &in, FermionField &out)
{
  // Wilson term
  out.checkerboard = in.checkerboard;
  this->Dhop(in, out, DaggerYes);
  // Clover term
  // apply the sigma and Fmunu
  FermionField temp(out._grid);
  MooeeDag(in, temp);
  out+=temp;
  return axpy_norm(out, 4 + this->mass, in, out);
}

template <class Impl>
void WilsonCloverFermion<Impl>::ImportGauge(const GaugeField &_Umu)
{
  WilsonFermion<Impl>::ImportGauge(_Umu);
  GridBase *grid = _Umu._grid;
  typename Impl::GaugeLinkField Bx(grid), By(grid), Bz(grid), Ex(grid), Ey(grid), Ez(grid);

  // Compute the field strength terms
  WilsonLoops<Impl>::FieldStrength(Bx, _Umu, Ydir, Zdir);
  WilsonLoops<Impl>::FieldStrength(By, _Umu, Zdir, Xdir);
  WilsonLoops<Impl>::FieldStrength(Bz, _Umu, Xdir, Ydir);
  WilsonLoops<Impl>::FieldStrength(Ex, _Umu, Tdir, Xdir);
  WilsonLoops<Impl>::FieldStrength(Ey, _Umu, Tdir, Ydir);
  WilsonLoops<Impl>::FieldStrength(Ez, _Umu, Tdir, Zdir);

  // Compute the Clover Operator acting on Colour and Spin
  CloverTerm  = fillCloverYZ(Bx);
  CloverTerm += fillCloverXZ(By);
  CloverTerm += fillCloverXY(Bz);
  CloverTerm += fillCloverXT(Ex);
  CloverTerm += fillCloverYT(Ey);
  CloverTerm += fillCloverZT(Ez) ;
  CloverTerm *= csw;

  int lvol = _Umu._grid->lSites();
  int DimRep = Impl::Dimension;

  Eigen::MatrixXcd EigenCloverOp = Eigen::MatrixXcd::Zero(Ns * DimRep, Ns * DimRep);
  Eigen::MatrixXcd EigenInvCloverOp = Eigen::MatrixXcd::Zero(Ns * DimRep, Ns * DimRep);

  std::vector<int> lcoor;
  typename SiteCloverType::scalar_object Qx = zero, Qxinv = zero;

  for (int site = 0; site < lvol; site++)
  {
    grid->LocalIndexToLocalCoor(site, lcoor);
    EigenCloverOp = Eigen::MatrixXcd::Zero(Ns * DimRep, Ns * DimRep);
    peekLocalSite(Qx, CloverTerm, lcoor);
    Qxinv = zero;
    for (int j = 0; j < Ns; j++)
      for (int k = 0; k < Ns; k++)
        for (int a = 0; a < DimRep; a++)
          for (int b = 0; b < DimRep; b++)
            EigenCloverOp(a + j * DimRep, b + k * DimRep) = Qx()(j, k)(a, b);
    //std::cout << EigenCloverOp << std::endl;


    EigenInvCloverOp = EigenCloverOp.inverse();
    //std::cout << EigenInvCloverOp << std::endl;
    for (int j = 0; j < Ns; j++)
      for (int k = 0; k < Ns; k++)
        for (int a = 0; a < DimRep; a++)
          for (int b = 0; b < DimRep; b++)
            Qxinv()(j, k)(a, b) = EigenInvCloverOp(a + j * DimRep, b + k * DimRep);

    pokeLocalSite(Qxinv, CloverTermInv, lcoor);
  }

  // Separate the even and odd parts.
  pickCheckerboard(Even, CloverTermEven, CloverTerm);
  pickCheckerboard( Odd,  CloverTermOdd, CloverTerm);

  pickCheckerboard(Even, CloverTermInvEven, CloverTermInv);
  pickCheckerboard( Odd,  CloverTermInvOdd, CloverTermInv);

}

template <class Impl>
void WilsonCloverFermion<Impl>::Mooee(const FermionField &in, FermionField &out)
{
  conformable(in,out);
  this->MooeeInternal(in, out, DaggerNo, InverseNo);
}

template <class Impl>
void WilsonCloverFermion<Impl>::MooeeDag(const FermionField &in, FermionField &out)
{
  this->MooeeInternal(in, out, DaggerNo, InverseYes);
}

template <class Impl>
void WilsonCloverFermion<Impl>::MooeeInv(const FermionField &in, FermionField &out)
{
  this->MooeeInternal(in, out, DaggerNo, InverseYes);
}

template <class Impl>
void WilsonCloverFermion<Impl>::MooeeInvDag(const FermionField &in, FermionField &out)
{
  this->MooeeInternal(in, out, DaggerNo, InverseYes);
}

template <class Impl>
void WilsonCloverFermion<Impl>::MooeeInternal(const FermionField &in, FermionField &out, int dag, int inv)
{
  out.checkerboard = in.checkerboard;
  CloverFieldType *Clover;
  assert(in.checkerboard == Odd || in.checkerboard == Even);

  if (in._grid->_isCheckerBoarded)
  {
    if (in.checkerboard == Odd)
    {
      std::cout << "Calling clover term Odd" << std::endl;
      Clover = (inv) ? &CloverTermInvOdd : &CloverTermOdd;
    }
    else
    {
      std::cout << "Calling clover term Even" << std::endl;
      Clover = (inv) ? &CloverTermInvEven : &CloverTermEven;
    }
  }
  else
  {
    Clover = (inv) ? &CloverTermInv : &CloverTerm;
  }

  std::cout << GridLogMessage << "*Clover.checkerboard "  << (*Clover).checkerboard << std::endl;
  if (dag){ out = adj(*Clover) * in;} else { out = *Clover * in;}
} // MooeeInternal

// Derivative parts
template <class Impl>
void WilsonCloverFermion<Impl>::MDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag)
{
  GaugeField tmp(mat._grid);

  conformable(U._grid, V._grid);
  conformable(U._grid, mat._grid);

  mat.checkerboard = U.checkerboard;
  tmp.checkerboard = U.checkerboard;

  this->DhopDeriv(mat, U, V, dag);
  MooDeriv(tmp, U, V, dag);
  mat += tmp;
}

// Derivative parts
template <class Impl>
void WilsonCloverFermion<Impl>::MooDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag)
{
  // Compute the 8 terms of the derivative
  assert(0); // not implemented yet
}

// Derivative parts
template <class Impl>
void WilsonCloverFermion<Impl>::MeeDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag)
{
  assert(0); // not implemented yet
}

FermOpTemplateInstantiate(WilsonCloverFermion); // now only for the fundamental representation
//AdjointFermOpTemplateInstantiate(WilsonCloverFermion);
//TwoIndexFermOpTemplateInstantiate(WilsonCloverFermion);
//GparityFermOpTemplateInstantiate(WilsonCloverFermion);
}
}
