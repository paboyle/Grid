/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MSource/Z2.hpp

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

#ifndef Hadrons_MSource_Z2_hpp_
#define Hadrons_MSource_Z2_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Z_2 stochastic source
 -----------------------------
 * src_x = eta_x * theta(x_3 - tA) * theta(tB - x_3)
 
 the eta_x are independent uniform random numbers in {+/- 1 +/- i}
 
 * options:
 - tA: begin timeslice (integer)
 - tB: end timesilce (integer)
 
 */
 
/******************************************************************************
 *                          Z2 stochastic source                              *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class Z2Par: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(Z2Par,
                                    unsigned int, tA,
                                    unsigned int, tB);
};

template <typename FImpl>
class TZ2: public Module<Z2Par>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TZ2(const std::string name);
    // destructor
    virtual ~TZ2(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool        hasT_{false};
    std::string tName_;
};

MODULE_REGISTER_TMP(Z2,       TZ2<FIMPL>,        MSource);
MODULE_REGISTER_TMP(ScalarZ2, TZ2<ScalarImplCR>, MSource);

/******************************************************************************
 *                       TZ2 template implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TZ2<FImpl>::TZ2(const std::string name)
: Module<Z2Par>(name)
, tName_ (name + "_t")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TZ2<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TZ2<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TZ2<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
    envCache(Lattice<iScalar<vInteger>>, tName_, 1, envGetGrid(LatticeComplex));
    envTmpLat(LatticeComplex, "eta");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TZ2<FImpl>::execute(void)
{
    if (par().tA == par().tB)
    {
        LOG(Message) << "Generating Z_2 wall source at t= " << par().tA
                     << std::endl;
    }
    else
    {
        LOG(Message) << "Generating Z_2 band for " << par().tA << " <= t <= "
                     << par().tB << std::endl;
    }
    
    auto    &src = envGet(PropagatorField, getName());
    auto    &t   = envGet(Lattice<iScalar<vInteger>>, tName_);
    Complex shift(1., 1.);

    if (!hasT_)
    {
        LatticeCoordinate(t, Tp);
        hasT_ = true;
    }
    envGetTmp(LatticeComplex, eta);
    bernoulli(rng4d(), eta);
    eta = (2.*eta - shift)*(1./::sqrt(2.));
    eta = where((t >= par().tA) and (t <= par().tB), eta, 0.*eta);
    src = 1.;
    src = src*eta;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_Z2_hpp_
