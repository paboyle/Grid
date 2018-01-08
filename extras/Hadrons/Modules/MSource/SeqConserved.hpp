/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MContraction/SeqConserved.hpp

Copyright (C) 2017

Author: Andrew Lawson    <andrew.lawson1991@gmail.com>

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

#ifndef Hadrons_MSource_SeqConserved_hpp_
#define Hadrons_MSource_SeqConserved_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Sequential source
 -----------------------------
 * src_x = q_x * theta(x_3 - tA) * theta(tB - x_3) * J_mu * exp(i x.mom)
 
 * options:
 - q: input propagator (string)
 - action: fermion action used for propagator q (string)
 - tA: begin timeslice (integer)
 - tB: end timesilce (integer)
 - curr_type: type of conserved current to insert (Current)
 - mu: Lorentz index of current to insert (integer)
 - mom: momentum insertion, space-separated float sequence (e.g ".1 .2 1. 0.")
 
 */

/******************************************************************************
 *                              SeqConserved                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class SeqConservedPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SeqConservedPar,
                                    std::string,  q,
                                    std::string,  action,
                                    unsigned int, tA,
                                    unsigned int, tB,
                                    Current,      curr_type,
                                    unsigned int, mu,
                                    std::string,  mom);
};

template <typename FImpl>
class TSeqConserved: public Module<SeqConservedPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TSeqConserved(const std::string name);
    // destructor
    virtual ~TSeqConserved(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(SeqConserved, TSeqConserved<FIMPL>, MSource);

/******************************************************************************
 *                      TSeqConserved implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSeqConserved<FImpl>::TSeqConserved(const std::string name)
: Module<SeqConservedPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSeqConserved<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().action};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSeqConserved<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqConserved<FImpl>::setup(void)
{
    auto Ls_ = env().getObjectLs(par().action);
    envCreateLat(PropagatorField, getName(), Ls_);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqConserved<FImpl>::execute(void)
{
    if (par().tA == par().tB)
    {
        LOG(Message) << "Generating sequential source with conserved "
                     << par().curr_type << " current insertion (mu = " 
                     << par().mu << ") at " << "t = " << par().tA << std::endl;
    }
    else
    {
        LOG(Message) << "Generating sequential source with conserved "
                     << par().curr_type << " current insertion (mu = " 
                     << par().mu << ") for " << par().tA << " <= t <= " 
                     << par().tB << std::endl;
    }
    auto &src = envGet(PropagatorField, getName());
    auto &q   = envGet(PropagatorField, par().q);
    auto &mat = envGet(FMat, par().action);

    std::vector<Real> mom = strToVec<Real>(par().mom);
    mat.SeqConservedCurrent(q, src, par().curr_type, par().mu, 
                            mom, par().tA, par().tB);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_SeqConserved_hpp_
