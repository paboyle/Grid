/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/action/fermion/WilsonCloverTypes.h

    Copyright (C) 2021 - 2022

    Author: Daniel Richtmann <daniel.richtmann@gmail.com>

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

NAMESPACE_BEGIN(Grid);

template<class Impl>
class WilsonCloverTypes {
public:
  INHERIT_IMPL_TYPES(Impl);

  template <typename vtype> using iImplClover = iScalar<iMatrix<iMatrix<vtype, Impl::Dimension>, Ns>>;

  typedef iImplClover<Simd> SiteClover;

  typedef Lattice<SiteClover> CloverField;
};

template<class Impl>
class CompactWilsonCloverTypes {
public:
  INHERIT_IMPL_TYPES(Impl);

  static_assert(Nd == 4 && Nc == 3 && Ns == 4 && Impl::Dimension == 3, "Wrong dimensions");

  static constexpr int Nred      = Nc * Nhs;        // 6
  static constexpr int Nblock    = Nhs;             // 2
  static constexpr int Ndiagonal = Nred;            // 6
  static constexpr int Ntriangle = (Nred - 1) * Nc; // 15

  template<typename vtype> using iImplCloverDiagonal = iScalar<iVector<iVector<vtype, Ndiagonal>, Nblock>>;
  template<typename vtype> using iImplCloverTriangle = iScalar<iVector<iVector<vtype, Ntriangle>, Nblock>>;

  typedef iImplCloverDiagonal<Simd> SiteCloverDiagonal;
  typedef iImplCloverTriangle<Simd> SiteCloverTriangle;
  typedef iSinglet<Simd>            SiteMask;

  typedef Lattice<SiteCloverDiagonal> CloverDiagonalField;
  typedef Lattice<SiteCloverTriangle> CloverTriangleField;
  typedef Lattice<SiteMask>           MaskField;
};

#define INHERIT_CLOVER_TYPES(Impl)                                 \
  typedef typename WilsonCloverTypes<Impl>::SiteClover SiteClover; \
  typedef typename WilsonCloverTypes<Impl>::CloverField CloverField;

#define INHERIT_COMPACT_CLOVER_TYPES(Impl) \
  typedef typename CompactWilsonCloverTypes<Impl>::SiteCloverDiagonal  SiteCloverDiagonal; \
  typedef typename CompactWilsonCloverTypes<Impl>::SiteCloverTriangle  SiteCloverTriangle; \
  typedef typename CompactWilsonCloverTypes<Impl>::SiteMask            SiteMask; \
  typedef typename CompactWilsonCloverTypes<Impl>::CloverDiagonalField CloverDiagonalField; \
  typedef typename CompactWilsonCloverTypes<Impl>::CloverTriangleField CloverTriangleField; \
  typedef typename CompactWilsonCloverTypes<Impl>::MaskField           MaskField; \
  /* ugly duplication but needed inside functionality classes */ \
  template<typename vtype> using iImplCloverDiagonal = \
    iScalar<iVector<iVector<vtype, CompactWilsonCloverTypes<Impl>::Ndiagonal>, CompactWilsonCloverTypes<Impl>::Nblock>>; \
  template<typename vtype> using iImplCloverTriangle = \
    iScalar<iVector<iVector<vtype, CompactWilsonCloverTypes<Impl>::Ntriangle>, CompactWilsonCloverTypes<Impl>::Nblock>>;

#define INHERIT_COMPACT_CLOVER_SIZES(Impl)                                    \
  static constexpr int Nred      = CompactWilsonCloverTypes<Impl>::Nred;      \
  static constexpr int Nblock    = CompactWilsonCloverTypes<Impl>::Nblock;    \
  static constexpr int Ndiagonal = CompactWilsonCloverTypes<Impl>::Ndiagonal; \
  static constexpr int Ntriangle = CompactWilsonCloverTypes<Impl>::Ntriangle;

NAMESPACE_END(Grid);
