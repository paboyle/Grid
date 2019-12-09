/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MSource/MomentumPhase.hpp

Copyright (C) 2015-2019

Author: Felix Erben <ferben@ed.ac.uk>

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

#ifndef Hadrons_MSource_MomentumPhase_hpp_
#define Hadrons_MSource_MomentumPhase_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Multiply source field by momentum phase
 -----------------------------
 * can be used to multiply a Z2-wall source by a phase e^{ipx}, which is needed for 
 * 2pt-correlation functions of mesons with nonzero momenta.
 
 */

/******************************************************************************
 *                         MomentumPhase                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class MomentumPhasePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MomentumPhasePar,
                                    std::string,    src,
                                    std::string,    mom);
};

template <typename FImpl>
class TMomentumPhase: public Module<MomentumPhasePar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TMomentumPhase(const std::string name);
    // destructor
    virtual ~TMomentumPhase(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool        hasPhase_{false};
    std::string momphName_, tName_;
};

MODULE_REGISTER_TMP(MomentumPhase, TMomentumPhase<FIMPL>, MSource);
MODULE_REGISTER_TMP(ZMomentumPhase, TMomentumPhase<ZFIMPL>, MSource);

/******************************************************************************
 *                         TMomentumPhase implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TMomentumPhase<FImpl>::TMomentumPhase(const std::string name)
: Module<MomentumPhasePar>(name)
, momphName_ (name + "_momph")
, tName_ (name + "_t")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TMomentumPhase<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().src};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TMomentumPhase<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMomentumPhase<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
    envCache(Lattice<iScalar<vInteger>>, tName_, 1, envGetGrid(LatticeComplex));
    envCacheLat(LatticeComplex, momphName_);
    envTmpLat(LatticeComplex, "coor");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMomentumPhase<FImpl>::execute(void)
{
    LOG(Message) << "Multiplying source field " << par().src
                 << " by momentum phase "
                 << par().mom << std::endl;
    auto  &out = envGet(PropagatorField, getName());
    auto  &src   = envGet(PropagatorField, par().src);
    auto  &ph  = envGet(LatticeComplex, momphName_);
    auto  &t   = envGet(Lattice<iScalar<vInteger>>, tName_);
    
    if (!hasPhase_)
    {
        Complex           i(0.0,1.0);
        std::vector<Real> p;

        envGetTmp(LatticeComplex, coor);
        p  = strToVec<Real>(par().mom);
        ph = Zero();
        for(unsigned int mu = 0; mu < env().getNd(); mu++)
        {
            LatticeCoordinate(coor, mu);
            ph = ph + (p[mu]/env().getDim(mu))*coor;
        }
        ph = exp((Real)(2*M_PI)*i*ph);
        LatticeCoordinate(t, Tp);
        hasPhase_ = true;
    }
    out = ph*src;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_MomentumPhase_hpp_
