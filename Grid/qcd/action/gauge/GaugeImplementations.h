/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/gauge/GaugeImplementations.h

Copyright (C) 2015

Author: paboyle <paboyle@ph.ed.ac.uk>

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
#ifndef GRID_QCD_GAUGE_IMPLEMENTATIONS_H
#define GRID_QCD_GAUGE_IMPLEMENTATIONS_H

#include "GaugeImplTypes.h"

NAMESPACE_BEGIN(Grid);

// Composition with smeared link, bc's etc.. probably need multiple inheritance
// Variable precision "S" and variable Nc
template <class GimplTypes> class PeriodicGaugeImpl : public GimplTypes {
public:
  INHERIT_GIMPL_TYPES(GimplTypes);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Support needed for the assembly of loops including all boundary condition
  // effects such as Conjugate bcs
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////

  template <class covariant>
  static inline Lattice<covariant>
  CovShiftForward(const GaugeLinkField &Link, int mu,
                  const Lattice<covariant> &field) {
    return PeriodicBC::CovShiftForward(Link, mu, field);
  }

  template <class covariant>
  static inline Lattice<covariant>
  CovShiftBackward(const GaugeLinkField &Link, int mu,
                   const Lattice<covariant> &field) {
    return PeriodicBC::CovShiftBackward(Link, mu, field);
  }
  static inline GaugeLinkField
  CovShiftIdentityBackward(const GaugeLinkField &Link, int mu) {
    return PeriodicBC::CovShiftIdentityBackward(Link, mu);
  }
  static inline GaugeLinkField
  CovShiftIdentityForward(const GaugeLinkField &Link, int mu) {
    return PeriodicBC::CovShiftIdentityForward(Link,mu);
  }
  static inline GaugeLinkField ShiftStaple(const GaugeLinkField &Link, int mu) {
    return PeriodicBC::ShiftStaple(Link,mu);
  }

  static inline bool isPeriodicGaugeField(void) { return true; }
};

// Composition with smeared link, bc's etc.. probably need multiple inheritance
// Variable precision "S" and variable Nc
class ConjugateGaugeImplBase {
protected:
  static std::vector<int> _conjDirs;
};

  template <class GimplTypes> class ConjugateGaugeImpl : public GimplTypes, ConjugateGaugeImplBase {
private:
public:
  INHERIT_GIMPL_TYPES(GimplTypes);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Support needed for the assembly of loops including all boundary condition
  // effects such as Gparity.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template <class covariant>
  static Lattice<covariant> CovShiftForward(const GaugeLinkField &Link, int mu,
                                            const Lattice<covariant> &field)
  {
    assert(_conjDirs.size() == Nd);
    if(_conjDirs[mu]) 
      return ConjugateBC::CovShiftForward(Link, mu, field);
    else
      return PeriodicBC::CovShiftForward(Link, mu, field);
  }

  template <class covariant>
  static Lattice<covariant> CovShiftBackward(const GaugeLinkField &Link, int mu,
                                             const Lattice<covariant> &field)
  {
    assert(_conjDirs.size() == Nd);
    if(_conjDirs[mu]) 
      return ConjugateBC::CovShiftBackward(Link, mu, field);
    else 
      return PeriodicBC::CovShiftBackward(Link, mu, field);
  }

  static inline GaugeLinkField
  CovShiftIdentityBackward(const GaugeLinkField &Link, int mu)
  {
    assert(_conjDirs.size() == Nd);
    if(_conjDirs[mu]) 
      return ConjugateBC::CovShiftIdentityBackward(Link, mu);
    else 
      return PeriodicBC::CovShiftIdentityBackward(Link, mu);
  }
  static inline GaugeLinkField
  CovShiftIdentityForward(const GaugeLinkField &Link, int mu)
  {
    assert(_conjDirs.size() == Nd);
    if(_conjDirs[mu]) 
      return ConjugateBC::CovShiftIdentityForward(Link,mu);
    else
      return PeriodicBC::CovShiftIdentityForward(Link,mu);
  }

  static inline GaugeLinkField ShiftStaple(const GaugeLinkField &Link, int mu)
  {
    assert(_conjDirs.size() == Nd);
    if(_conjDirs[mu]) 
      return ConjugateBC::ShiftStaple(Link,mu);
    else     
      return PeriodicBC::ShiftStaple(Link,mu);
  }

  static inline void       setDirections(std::vector<int> &conjDirs) { _conjDirs=conjDirs; }
  static inline std::vector<int> getDirections(void) { return _conjDirs; }
  static inline bool isPeriodicGaugeField(void) { return false; }
};

typedef PeriodicGaugeImpl<GimplTypesR> PeriodicGimplR; // Real.. whichever prec
typedef PeriodicGaugeImpl<GimplTypesF> PeriodicGimplF; // Float
typedef PeriodicGaugeImpl<GimplTypesD> PeriodicGimplD; // Double

typedef PeriodicGaugeImpl<GimplAdjointTypesR> PeriodicGimplAdjR; // Real.. whichever prec
typedef PeriodicGaugeImpl<GimplAdjointTypesF> PeriodicGimplAdjF; // Float
typedef PeriodicGaugeImpl<GimplAdjointTypesD> PeriodicGimplAdjD; // Double

typedef ConjugateGaugeImpl<GimplTypesR> ConjugateGimplR; // Real.. whichever prec
typedef ConjugateGaugeImpl<GimplTypesF> ConjugateGimplF; // Float
typedef ConjugateGaugeImpl<GimplTypesD> ConjugateGimplD; // Double

NAMESPACE_END(Grid);

#endif
