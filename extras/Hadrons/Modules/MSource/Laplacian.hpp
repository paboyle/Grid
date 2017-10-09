/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MSource/Laplacian.hpp

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

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

#ifndef Hadrons_MSource_Laplacian_hpp_
#define Hadrons_MSource_Laplacian_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Laplacian smearing source
 -----------------------------
 
 * options:
 - source: name of source object to be smeared (string)
 - N: number of steps (integer)
 - alpha: smearing parameter (real)
 
 */

/******************************************************************************
 *                          Laplace smearing operator                         *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class LaplacianPar : Serializable
{
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LaplacianPar,
                                    std::string, source,
                                    std::string, gauge,
                                    unsigned int, N,
                                    double, alpha);
};

template <typename FImpl>
class TLaplacian : public Module<LaplacianPar>
{
  public:
    FERM_TYPE_ALIASES(FImpl, );

  public:
    // constructor
    TLaplacian(const std::string name);
    // destructor
    virtual ~TLaplacian(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(LaplaceSmearing, TLaplacian<FIMPL>, MSource);

/******************************************************************************
 *                       TLaplacian template implementation                   *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLaplacian<FImpl>::TLaplacian(const std::string name)
    : Module<LaplacianPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLaplacian<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().source, par().gauge};

    return in;
}

template <typename FImpl>
std::vector<std::string> TLaplacian<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLaplacian<FImpl>::setup(void)
{
    env().template registerLattice<PropagatorField>(getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLaplacian<FImpl>::execute(void)
{

    FermionField source(env().getGrid()), tmp(env().getGrid());
    PropagatorField &SmrSrc = *env().template createLattice<PropagatorField>(getName());
    PropagatorField &fullSrc = *env().template getObject<PropagatorField>(par().source);
    auto &U      = *env().template getObject<LatticeGaugeField>(par().gauge);
    Laplacian<FImpl> LaplaceOperator(env().getGrid());
    LaplaceOperator.ImportGauge(U);
    double prefactor = par().alpha / (double)(par().N);

    for (unsigned int s = 0; s < Ns; ++s)
    {
        for (unsigned int c = 0; c < Nc; ++c)
        {
            PropToFerm(source, fullSrc, s, c);
            for (int smr = 0; smr < par().N; ++smr)
            {
                LaplaceOperator.M(source, tmp);
                source += prefactor * tmp;
            }
            FermToProp(SmrSrc, source, s, c);
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_Z2_hpp_
