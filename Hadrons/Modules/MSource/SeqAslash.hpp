/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MSource/SeqAslash.hpp

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Lanny91 <andrew.lawson@gmail.com>
Author: Vera Guelpers <Vera.Guelpers@ed.ac.uk>

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

#ifndef Hadrons_MSource_SeqAslash_hpp_
#define Hadrons_MSource_SeqAslash_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Sequential source
 -----------------------------
 * src_x = q_x * theta(x_3 - tA) * theta(tB - x_3) * i * A_mu g_mu * exp(i x.mom)
 
 * options:
 - q: input propagator (string)
 - tA: begin timeslice (integer)
 - tB: end timesilce (integer)
 - emField: input photon field (string)
 - mom: momentum insertion, space-separated float sequence (e.g ".1 .2 1. 0.")
 
 */

/******************************************************************************
 *                         SeqAslash                             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class SeqAslashPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SeqAslashPar,
                                    std::string,    q,
                                    unsigned int,   tA,
                                    unsigned int,   tB,
                                    std::string,    emField,
                                    std::string,    mom);
};

template <typename FImpl>
class TSeqAslash: public Module<SeqAslashPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    typedef PhotonR::GaugeField     EmField;
public:
    // constructor
    TSeqAslash(const std::string name);
    // destructor
    virtual ~TSeqAslash(void) {};
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

MODULE_REGISTER_TMP(SeqAslash, TSeqAslash<FIMPL>, MSource);

/******************************************************************************
 *                         TSeqAslash implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSeqAslash<FImpl>::TSeqAslash(const std::string name)
: Module<SeqAslashPar>(name)
, momphName_ (name + "_momph")
, tName_ (name + "_t")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSeqAslash<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q,par().emField};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSeqAslash<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqAslash<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
    envCache(Lattice<iScalar<vInteger>>, tName_, 1, envGetGrid(LatticeComplex));
    envCacheLat(LatticeComplex, momphName_);
    envTmpLat(LatticeComplex, "coor");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqAslash<FImpl>::execute(void)
{
    if (par().tA == par().tB)
    {
        LOG(Message) << "Generating Aslash sequential source at t= " << par().tA 
		     << " using the photon field " << par().emField << std::endl; 
    }
    else
    {
        LOG(Message) << "Generating Aslash sequential source for "
                     << par().tA << " <= t <= " << par().tB 
		     << " using the photon field " << par().emField << std::endl;
    }
    auto  &src = envGet(PropagatorField, getName()); src=zero;
    auto  &q   = envGet(PropagatorField, par().q);
    auto  &ph  = envGet(LatticeComplex, momphName_);
    auto  &t   = envGet(Lattice<iScalar<vInteger>>, tName_);

    if (!hasPhase_)
    {
        Complex           i(0.0,1.0);
        std::vector<Real> p;

        envGetTmp(LatticeComplex, coor);
        p  = strToVec<Real>(par().mom);
        ph = zero;
        for(unsigned int mu = 0; mu < env().getNd(); mu++)
        {
            LatticeCoordinate(coor, mu);
            ph = ph + (p[mu]/env().getDim(mu))*coor;
        }
        ph = exp((Real)(2*M_PI)*i*ph);
        LatticeCoordinate(t, Tp);
        hasPhase_ = true;
    }
    auto &stoch_photon = envGet(EmField,  par().emField);
    Complex ci(0.0,1.0);
    for(unsigned int mu=0;mu<=3;mu++)
    {
	Gamma gmu(Gamma::gmu[mu]);
	src = src + where((t >= par().tA) and (t <= par().tB), ci * PeekIndex<LorentzIndex>(stoch_photon, mu) *ph*(gmu*q), 0.*q);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_SeqAslash_hpp_
