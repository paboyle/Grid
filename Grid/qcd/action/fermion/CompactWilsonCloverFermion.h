/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/action/fermion/CompactWilsonCloverFermion.h

    Copyright (C) 2020 - 2022

    Author: Daniel Richtmann <daniel.richtmann@gmail.com>
    Author: Nils Meyer <nils.meyer@ur.de>

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

#include <Grid/qcd/action/fermion/WilsonCloverTypes.h>
#include <Grid/qcd/action/fermion/WilsonCloverHelpers.h>

NAMESPACE_BEGIN(Grid);

// see Grid/qcd/action/fermion/WilsonCloverFermion.h for description
//
// Modifications done here:
//
// Original: clover term = 12x12 matrix per site
//
// But: Only two diagonal 6x6 hermitian blocks are non-zero (also true for original, verified by running)
// Sufficient to store/transfer only the real parts of the diagonal and one triangular part
// 2 * (6 + 15 * 2) = 72 real or 36 complex words to be stored/transfered
//
// Here: Above but diagonal as complex numbers, i.e., need to store/transfer
// 2 * (6 * 2 + 15 * 2) = 84 real or 42 complex words
//
// Words per site and improvement compared to original (combined with the input and output spinors):
//
// - Original: 2*12 + 12*12 = 168 words -> 1.00 x less
// - Minimal:  2*12 + 36    =  60 words -> 2.80 x less
// - Here:     2*12 + 42    =  66 words -> 2.55 x less
//
// These improvements directly translate to wall-clock time
//
// Data layout:
//
// - diagonal and triangle part as separate lattice fields,
//   this was faster than as 1 combined field on all tested machines
// - diagonal: as expected
// - triangle: store upper right triangle in row major order
// - graphical:
//        0  1  2  3  4
//           5  6  7  8
//              9 10 11 = upper right triangle indices
//                12 13
//                   14
//     0
//        1
//           2
//              3       = diagonal indices
//                 4
//                    5
//     0
//     1  5
//     2  6  9          = lower left triangle indices
//     3  7 10 12
//     4  8 11 13 14
//
// Impact on total memory consumption:
// - Original: (2 * 1 + 8 * 1/2) 12x12 matrices = 6 12x12 matrices = 864 complex words per site
// - Here:     (2 * 1 + 4 * 1/2) diagonal parts = 4 diagonal parts =  24 complex words per site
//           + (2 * 1 + 4 * 1/2) triangle parts = 4 triangle parts =  60 complex words per site
//                                                                 =  84 complex words per site

template<class Impl>
class CompactWilsonCloverFermion : public WilsonFermion<Impl>,
                                   public WilsonCloverHelpers<Impl>,
                                   public CompactWilsonCloverHelpers<Impl> {
  /////////////////////////////////////////////
  // Sizes
  /////////////////////////////////////////////

public:

  INHERIT_COMPACT_CLOVER_SIZES(Impl);

  /////////////////////////////////////////////
  // Type definitions
  /////////////////////////////////////////////

public:

  INHERIT_IMPL_TYPES(Impl);
  INHERIT_CLOVER_TYPES(Impl);
  INHERIT_COMPACT_CLOVER_TYPES(Impl);

  typedef WilsonFermion<Impl>              WilsonBase;
  typedef WilsonCloverHelpers<Impl>        Helpers;
  typedef CompactWilsonCloverHelpers<Impl> CompactHelpers;

  /////////////////////////////////////////////
  // Constructors
  /////////////////////////////////////////////

public:

  CompactWilsonCloverFermion(GaugeField& _Umu,
			    GridCartesian& Fgrid,
			    GridRedBlackCartesian& Hgrid,
			    const RealD _mass,
			    const RealD _csw_r = 0.0,
			    const RealD _csw_t = 0.0,
			    const RealD _cF = 1.0,
			    const WilsonAnisotropyCoefficients& clover_anisotropy = WilsonAnisotropyCoefficients(),
			    const ImplParams& impl_p = ImplParams())
    : WilsonBase(_Umu, Fgrid, Hgrid, _mass, impl_p, clover_anisotropy)
    , csw_r(_csw_r)
    , csw_t(_csw_t)
    , cF(_cF)
    , open_boundaries(impl_p.boundary_phases[Nd-1] == 0.0)
    , Diagonal(&Fgrid),        Triangle(&Fgrid)
    , DiagonalEven(&Hgrid),    TriangleEven(&Hgrid)
    , DiagonalOdd(&Hgrid),     TriangleOdd(&Hgrid)
    , DiagonalInv(&Fgrid),     TriangleInv(&Fgrid)
    , DiagonalInvEven(&Hgrid), TriangleInvEven(&Hgrid)
    , DiagonalInvOdd(&Hgrid),  TriangleInvOdd(&Hgrid)
    , Tmp(&Fgrid)
    , BoundaryMask(&Fgrid)
    , BoundaryMaskEven(&Hgrid), BoundaryMaskOdd(&Hgrid)
  {
    csw_r *= 0.5;
    csw_t *= 0.5;
    if (clover_anisotropy.isAnisotropic)
      csw_r /= clover_anisotropy.xi_0;

    ImportGauge(_Umu);
    if (open_boundaries)
      CompactHelpers::SetupMasks(this->BoundaryMask, this->BoundaryMaskEven, this->BoundaryMaskOdd);
  }

  /////////////////////////////////////////////
  // Member functions (implementing interface)
  /////////////////////////////////////////////

public:

  virtual void Instantiatable() {};
  int          ConstEE()     override { return 0; };
  int          isTrivialEE() override { return 0; };


  void Dhop(const FermionField& in, FermionField& out, int dag) override {
    WilsonBase::Dhop(in, out, dag);
    if(open_boundaries) ApplyBoundaryMask(out);
  }

  void DhopOE(const FermionField& in, FermionField& out, int dag) override {
    WilsonBase::DhopOE(in, out, dag);
    if(open_boundaries) ApplyBoundaryMask(out);
  }

  void DhopEO(const FermionField& in, FermionField& out, int dag) override {
    WilsonBase::DhopEO(in, out, dag);
    if(open_boundaries) ApplyBoundaryMask(out);
  }

  void DhopDir(const FermionField& in, FermionField& out, int dir, int disp) override {
    WilsonBase::DhopDir(in, out, dir, disp);
    if(this->open_boundaries) ApplyBoundaryMask(out);
  }

  void DhopDirAll(const FermionField& in, std::vector<FermionField>& out) /* override */ {
    WilsonBase::DhopDirAll(in, out);
    if(this->open_boundaries) {
      for(auto& o : out) ApplyBoundaryMask(o);
    }
  }

  void M(const FermionField& in, FermionField& out) override {
    out.Checkerboard() = in.Checkerboard();
    WilsonBase::Dhop(in, out, DaggerNo); // call base to save applying bc
    Mooee(in, Tmp);
    axpy(out, 1.0, out, Tmp);
    if(open_boundaries) ApplyBoundaryMask(out);
  }

  void Mdag(const FermionField& in, FermionField& out) override {
    out.Checkerboard() = in.Checkerboard();
    WilsonBase::Dhop(in, out, DaggerYes);  // call base to save applying bc
    MooeeDag(in, Tmp);
    axpy(out, 1.0, out, Tmp);
    if(open_boundaries) ApplyBoundaryMask(out);
  }

  void Meooe(const FermionField& in, FermionField& out) override {
    WilsonBase::Meooe(in, out);
    if(open_boundaries) ApplyBoundaryMask(out);
  }

  void MeooeDag(const FermionField& in, FermionField& out) override {
    WilsonBase::MeooeDag(in, out);
    if(open_boundaries) ApplyBoundaryMask(out);
  }

  void Mooee(const FermionField& in, FermionField& out) override {
    if(in.Grid()->_isCheckerBoarded) {
      if(in.Checkerboard() == Odd) {
        MooeeInternal(in, out, DiagonalOdd, TriangleOdd);
      } else {
        MooeeInternal(in, out, DiagonalEven, TriangleEven);
      }
    } else {
      MooeeInternal(in, out, Diagonal, Triangle);
    }
    if(open_boundaries) ApplyBoundaryMask(out);
  }

  void MooeeDag(const FermionField& in, FermionField& out) override {
    Mooee(in, out); // blocks are hermitian
  }

  void MooeeInv(const FermionField& in, FermionField& out) override {
    if(in.Grid()->_isCheckerBoarded) {
      if(in.Checkerboard() == Odd) {
        MooeeInternal(in, out, DiagonalInvOdd, TriangleInvOdd);
      } else {
        MooeeInternal(in, out, DiagonalInvEven, TriangleInvEven);
      }
    } else {
      MooeeInternal(in, out, DiagonalInv, TriangleInv);
    }
    if(open_boundaries) ApplyBoundaryMask(out);
  }

  void MooeeInvDag(const FermionField& in, FermionField& out) override {
    MooeeInv(in, out); // blocks are hermitian
  }

  void Mdir(const FermionField& in, FermionField& out, int dir, int disp) override {
    DhopDir(in, out, dir, disp);
  }

  void MdirAll(const FermionField& in, std::vector<FermionField>& out) override {
    DhopDirAll(in, out);
  }

  void MDeriv(GaugeField& force, const FermionField& X, const FermionField& Y, int dag) override {
    assert(!open_boundaries); // TODO check for changes required for open bc

    // NOTE: code copied from original clover term
    conformable(X.Grid(), Y.Grid());
    conformable(X.Grid(), force.Grid());
    GaugeLinkField force_mu(force.Grid()), lambda(force.Grid());
    GaugeField clover_force(force.Grid());
    PropagatorField Lambda(force.Grid());

    // Guido: Here we are hitting some performance issues:
    // need to extract the components of the DoubledGaugeField
    // for each call
    // Possible solution
    // Create a vector object to store them? (cons: wasting space)
    std::vector<GaugeLinkField> U(Nd, this->Umu.Grid());

    Impl::extractLinkField(U, this->Umu);

    force = Zero();
    // Derivative of the Wilson hopping term
    this->DhopDeriv(force, X, Y, dag);

    ///////////////////////////////////////////////////////////
    // Clover term derivative
    ///////////////////////////////////////////////////////////
    Impl::outerProductImpl(Lambda, X, Y);
    //std::cout << "Lambda:" << Lambda << std::endl;

    Gamma::Algebra sigma[] = {
        Gamma::Algebra::SigmaXY,
        Gamma::Algebra::SigmaXZ,
        Gamma::Algebra::SigmaXT,
        Gamma::Algebra::MinusSigmaXY,
        Gamma::Algebra::SigmaYZ,
        Gamma::Algebra::SigmaYT,
        Gamma::Algebra::MinusSigmaXZ,
        Gamma::Algebra::MinusSigmaYZ,
        Gamma::Algebra::SigmaZT,
        Gamma::Algebra::MinusSigmaXT,
        Gamma::Algebra::MinusSigmaYT,
        Gamma::Algebra::MinusSigmaZT};

    /*
      sigma_{\mu \nu}=
      | 0         sigma[0]  sigma[1]  sigma[2] |
      | sigma[3]    0       sigma[4]  sigma[5] |
      | sigma[6]  sigma[7]     0      sigma[8] |
      | sigma[9]  sigma[10] sigma[11]   0      |
    */

    int count = 0;
    clover_force = Zero();
    for (int mu = 0; mu < 4; mu++)
    {
      force_mu = Zero();
      for (int nu = 0; nu < 4; nu++)
      {
        if (mu == nu)
          continue;

        RealD factor;
        if (nu == 4 || mu == 4)
        {
          factor = 2.0 * csw_t;
        }
        else
        {
          factor = 2.0 * csw_r;
        }
        PropagatorField Slambda = Gamma(sigma[count]) * Lambda; // sigma checked
        Impl::TraceSpinImpl(lambda, Slambda);                   // traceSpin ok
        force_mu -= factor*Helpers::Cmunu(U, lambda, mu, nu);   // checked
        count++;
      }

      pokeLorentz(clover_force, U[mu] * force_mu, mu);
    }
    //clover_force *= csw;
    force += clover_force;
  }

  void MooDeriv(GaugeField& mat, const FermionField& U, const FermionField& V, int dag) override {
    assert(0);
  }

  void MeeDeriv(GaugeField& mat, const FermionField& U, const FermionField& V, int dag) override {
    assert(0);
  }

  /////////////////////////////////////////////
  // Member functions (internals)
  /////////////////////////////////////////////

  void MooeeInternal(const FermionField&        in,
                     FermionField&              out,
                     const CloverDiagonalField& diagonal,
                     const CloverTriangleField& triangle) {
    assert(in.Checkerboard() == Odd || in.Checkerboard() == Even);
    out.Checkerboard() = in.Checkerboard();
    conformable(in, out);
    conformable(in, diagonal);
    conformable(in, triangle);

    CompactHelpers::MooeeKernel(diagonal.oSites(), 1, in, out, diagonal, triangle);
  }

  /////////////////////////////////////////////
  // Helpers
  /////////////////////////////////////////////

  void ImportGauge(const GaugeField& _Umu) override {
    // NOTE: parts copied from original implementation

    // Import gauge into base class
    double t0 = usecond();
    WilsonBase::ImportGauge(_Umu); // NOTE: called here and in wilson constructor -> performed twice, but can't avoid that

    // Initialize temporary variables
    double t1 = usecond();
    conformable(_Umu.Grid(), this->GaugeGrid());
    GridBase* grid = _Umu.Grid();
    typename Impl::GaugeLinkField Bx(grid), By(grid), Bz(grid), Ex(grid), Ey(grid), Ez(grid);
    CloverField TmpOriginal(grid);

    // Compute the field strength terms mu>nu
    double t2 = usecond();
    WilsonLoops<Impl>::FieldStrength(Bx, _Umu, Zdir, Ydir);
    WilsonLoops<Impl>::FieldStrength(By, _Umu, Zdir, Xdir);
    WilsonLoops<Impl>::FieldStrength(Bz, _Umu, Ydir, Xdir);
    WilsonLoops<Impl>::FieldStrength(Ex, _Umu, Tdir, Xdir);
    WilsonLoops<Impl>::FieldStrength(Ey, _Umu, Tdir, Ydir);
    WilsonLoops<Impl>::FieldStrength(Ez, _Umu, Tdir, Zdir);

    // Compute the Clover Operator acting on Colour and Spin
    // multiply here by the clover coefficients for the anisotropy
    double t3 = usecond();
    TmpOriginal  = Helpers::fillCloverYZ(Bx) * csw_r;
    TmpOriginal += Helpers::fillCloverXZ(By) * csw_r;
    TmpOriginal += Helpers::fillCloverXY(Bz) * csw_r;
    TmpOriginal += Helpers::fillCloverXT(Ex) * csw_t;
    TmpOriginal += Helpers::fillCloverYT(Ey) * csw_t;
    TmpOriginal += Helpers::fillCloverZT(Ez) * csw_t;
    TmpOriginal += this->diag_mass;

    // Convert the data layout of the clover term
    double t4 = usecond();
    CompactHelpers::ConvertLayout(TmpOriginal, Diagonal, Triangle);

    // Possible modify the boundary values
    double t5 = usecond();
    if(open_boundaries) CompactHelpers::ModifyBoundaries(Diagonal, Triangle, csw_t, cF, this->diag_mass);

    // Invert the clover term in the improved layout
    double t6 = usecond();
    CompactHelpers::Invert(Diagonal, Triangle, DiagonalInv, TriangleInv);

    // Fill the remaining clover fields
    double t7 = usecond();
    pickCheckerboard(Even, DiagonalEven,    Diagonal);
    pickCheckerboard(Even, TriangleEven,    Triangle);
    pickCheckerboard(Odd,  DiagonalOdd,     Diagonal);
    pickCheckerboard(Odd,  TriangleOdd,     Triangle);
    pickCheckerboard(Even, DiagonalInvEven, DiagonalInv);
    pickCheckerboard(Even, TriangleInvEven, TriangleInv);
    pickCheckerboard(Odd,  DiagonalInvOdd,  DiagonalInv);
    pickCheckerboard(Odd,  TriangleInvOdd,  TriangleInv);

    // Report timings
    double t8 = usecond();
#if 0
    std::cout << GridLogMessage << "CompactWilsonCloverFermion::ImportGauge timings:"
              << " WilsonFermion::Importgauge = " << (t1 - t0) / 1e6
              << ", allocations = "               << (t2 - t1) / 1e6
              << ", field strength = "            << (t3 - t2) / 1e6
              << ", fill clover = "               << (t4 - t3) / 1e6
              << ", convert = "                   << (t5 - t4) / 1e6
              << ", boundaries = "                << (t6 - t5) / 1e6
              << ", inversions = "                << (t7 - t6) / 1e6
              << ", pick cbs = "                  << (t8 - t7) / 1e6
              << ", total = "                     << (t8 - t0) / 1e6
              << std::endl;
#endif
  }

  /////////////////////////////////////////////
  // Helpers
  /////////////////////////////////////////////

private:

  template<class Field>
  const MaskField* getCorrectMaskField(const Field &in) const {
    if(in.Grid()->_isCheckerBoarded) {
      if(in.Checkerboard() == Odd) {
        return &this->BoundaryMaskOdd;
      } else {
        return &this->BoundaryMaskEven;
      }
    } else {
      return &this->BoundaryMask;
    }
  }

  template<class Field>
  void ApplyBoundaryMask(Field& f) {
    const MaskField* m = getCorrectMaskField(f); assert(m != nullptr);
    assert(m != nullptr);
    CompactHelpers::ApplyBoundaryMask(f, *m);
  }

  /////////////////////////////////////////////
  // Member Data
  /////////////////////////////////////////////

public:

  RealD csw_r;
  RealD csw_t;
  RealD cF;

  bool open_boundaries;

  CloverDiagonalField Diagonal,    DiagonalEven,    DiagonalOdd;
  CloverDiagonalField DiagonalInv, DiagonalInvEven, DiagonalInvOdd;

  CloverTriangleField Triangle,    TriangleEven,    TriangleOdd;
  CloverTriangleField TriangleInv, TriangleInvEven, TriangleInvOdd;

  FermionField Tmp;

  MaskField BoundaryMask, BoundaryMaskEven, BoundaryMaskOdd;
};

NAMESPACE_END(Grid);
