/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Utilities/Contractor.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>

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
#ifndef  Hadrons_Contractor_hpp_
#define Hadrons_Contractor_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE



END_HADRONS_NAMESPACE

#define BEGIN_CONTRACTOR_NAMESPACE namespace Contractor{
BEGIN_CONTRACTOR_NAMESPACE

using Grid::Serializable;
using Grid::Reader;
using Grid::Writer;
using Grid::ComplexD;

class TrajRange: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(TrajRange,
                                  unsigned int, start,
                                  unsigned int, end,
                                  unsigned int, step);
};

class GlobalPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(GlobalPar,
                                  TrajRange, trajCounter,
                                  unsigned int, nt,
                                  std::string, diskVectorDir,
                                  std::string, output);
};

class A2AMatrixPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMatrixPar,
                                  std::string, file,
                                  std::string, dataset,
                                  unsigned int, cacheSize,
                                  std::string, name);
};

class ProductPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(ProductPar,
                                  std::string, terms,
                                  std::vector<std::string>, times,
                                  std::string, translations,
                                  bool, translationAverage);
};

class CorrelatorResult: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(CorrelatorResult,
                                  std::vector<Contractor::A2AMatrixPar>,  a2aMatrix,
                                  ProductPar, contraction,
                                  std::vector<unsigned int>, times,
                                  std::vector<ComplexD>, correlator);
};

struct ContractorPar
{
  Contractor::GlobalPar                  global;
  std::vector<Contractor::A2AMatrixPar>  a2aMatrix;
  std::vector<Contractor::ProductPar>    product;
};

// Useful ... so long as there's a ContractorPar named par in scope
#define TIME_MOD(t) (((t) + par.global.nt) % par.global.nt)

#define END_CONTRACTOR_NAMESPACE }
END_CONTRACTOR_NAMESPACE

#endif // Hadrons_Contractor_hpp_
