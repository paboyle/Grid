/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/scalar/CovariantLaplacian.h

Copyright (C) 2016

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */

#ifndef COVARIANT_LAPLACIAN_H
#define COVARIANT_LAPLACIAN_H

namespace Grid
{
namespace QCD
{

////////////////////////////////////////////////////////////
// Laplacian operator L on fermion fields
//
// phi: fermion field
//
// L phi(x) = Sum_mu [ U_mu(x) phi(x+mu) + U_mu(x-mu) phi(x-mu) - 2phi(x)]
//
// Operator designed to be encapsulated by
// an HermitianLinearOperator<.. , ..>
////////////////////////////////////////////////////////////

// has to inherit from a fermion implementation
template <class Impl>
class Laplacian
{
public:
  INHERIT_IMPL_TYPES(Impl);

  // add a bool to smear only in the spatial directions
  Laplacian(GridBase *grid, bool spatial = false)
      : U(Nd, grid), spatial_laplacian(spatial){};

  void Mdir(const FermionField &, FermionField &, int, int) { assert(0); }
  void Mdiag(const FermionField &, FermionField &) { assert(0); }

  void ImportGauge(const GaugeField &_U)
  {
    for (int mu = 0; mu < Nd; mu++)
      U[mu] = PeekIndex<LorentzIndex>(_U, mu);
  }

  void M(const FermionField &in, FermionField &out)
  {
    int dims = spatial_laplacian ? (Nd - 1) : Nd;

    out = -2.0 * dims * in;
    // eventually speed up with the stencil operator, if necessary
    for (int mu = 0; mu < dims; mu++)
      out += Impl::CovShiftForward(U[mu], mu, in) + Impl::CovShiftBackward(U[mu], mu, in);
  }

private:
  bool spatial_laplacian;
  std::vector<GaugeLinkField> U;
}; // Laplacian

} // QCD
} // Grid
#endif
