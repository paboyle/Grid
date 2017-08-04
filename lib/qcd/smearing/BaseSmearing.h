/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/modules/plaquette.h

Copyright (C) 2017

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
/*
  @brief Declares base smearing class Smear
 */
#ifndef BASE_SMEAR_
#define BASE_SMEAR_

template <class Gimpl>
class Smear{
public:
  INHERIT_GIMPL_TYPES(Gimpl) // inherits the types for the gauge fields

  virtual ~Smear(){}
  virtual void smear     (GaugeField&,const GaugeField&)const = 0;
  virtual void derivative(GaugeField&, const GaugeField&,const GaugeField&) const = 0;
};
#endif
