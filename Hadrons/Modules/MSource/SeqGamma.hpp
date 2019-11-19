/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MSource/SeqGamma.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Lanny91 <andrew.lawson@gmail.com>

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

#ifndef Hadrons_MSource_SeqGamma_hpp_
#define Hadrons_MSource_SeqGamma_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Sequential source
 -----------------------------
 * src_x = q_x * theta(x_3 - tA) * theta(tB - x_3) * gamma * exp(i x.mom)
 
 * options:
 - q: input propagator (string)
 - tA: begin timeslice (integer)
 - tB: end timesilce (integer)
 - gamma: gamma product to insert (integer)
 - mom: momentum insertion, space-separated float sequence (e.g ".1 .2 1. 0.")
 
 */

/******************************************************************************
 *                         SeqGamma                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class SeqGammaPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SeqGammaPar,
                                    std::string,    q,
                                    unsigned int,   tA,
                                    unsigned int,   tB,
                                    Gamma::Algebra, gamma,
                                    std::string,    mom);
};

template <typename FImpl>
class TSeqGamma: public Module<SeqGammaPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TSeqGamma(const std::string name);
    // destructor
    virtual ~TSeqGamma(void) {};
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

template <typename FImpl>
class TStagSeqGamma: public Module<SeqGammaPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TStagSeqGamma(const std::string name);
    // destructor
    virtual ~TStagSeqGamma(void) {};
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
    std::string momphName_, tName_, stagPhaseName_;
};

MODULE_REGISTER_TMP(SeqGamma, TSeqGamma<FIMPL>, MSource);
MODULE_REGISTER_TMP(ZSeqGamma, TSeqGamma<ZFIMPL>, MSource);
MODULE_REGISTER_TMP(StagSeqGamma, TStagSeqGamma<STAGIMPL>, MSource);

/******************************************************************************
 *                         TSeqGamma implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSeqGamma<FImpl>::TSeqGamma(const std::string name)
: Module<SeqGammaPar>(name)
, momphName_ (name + "_momph")
, tName_ (name + "_t")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSeqGamma<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSeqGamma<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqGamma<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
    envCache(Lattice<iScalar<vInteger>>, tName_, 1, envGetGrid(LatticeComplex));
    envCacheLat(LatticeComplex, momphName_);
    envTmpLat(LatticeComplex, "coor");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqGamma<FImpl>::execute(void)
{
    if (par().tA == par().tB)
    {
        LOG(Message) << "Generating gamma_" << par().gamma
                     << " sequential source at t= " << par().tA << std::endl;
    }
    else
    {
        LOG(Message) << "Generating gamma_" << par().gamma
                     << " sequential source for "
                     << par().tA << " <= t <= " << par().tB << std::endl;
    }
    auto  &src = envGet(PropagatorField, getName());
    auto  &q   = envGet(PropagatorField, par().q);
    auto  &ph  = envGet(LatticeComplex, momphName_);
    auto  &t   = envGet(Lattice<iScalar<vInteger>>, tName_);
    Gamma g(par().gamma);
    
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
    src = where((t >= par().tA) and (t <= par().tB), ph*(g*q), 0.*q);
}

/******************************************************************************
 *                         TStagSeqGamma implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TStagSeqGamma<FImpl>::TStagSeqGamma(const std::string name)
: Module<SeqGammaPar>(name)
, momphName_ (name + "_momph")
, tName_ (name + "_t")
, stagPhaseName_ (name + "_stagphase")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TStagSeqGamma<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TStagSeqGamma<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStagSeqGamma<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
    envCache(Lattice<iScalar<vInteger>>, tName_, 1, envGetGrid(LatticeComplex));
    envCacheLat(LatticeComplex, momphName_);
    envCacheLat(LatticeComplex, stagPhaseName_);
    envTmpLat(LatticeComplex, "coor");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStagSeqGamma<FImpl>::execute(void)
{
    if (par().tA == par().tB)
    {
        LOG(Message) << "Generating gamma_" << par().gamma
        << " sequential source at t= " << par().tA << std::endl;
        LOG(Message) << "Using propagator " << par().q << std::endl;
    }
    else
    {
        LOG(Message) << "Generating gamma_" << par().gamma
        << " sequential source for "
        << par().tA << " <= t <= " << par().tB << std::endl;
        LOG(Message) << "Using propagator " << par().q << std::endl;
    }
    auto  &src = envGet(PropagatorField, getName());
    auto  &q   = envGet(PropagatorField, par().q);
    auto  &ph  = envGet(LatticeComplex, momphName_);
    auto  &stag_ph  = envGet(LatticeComplex, stagPhaseName_);
    Lattice<iScalar<vInteger> > t(env().getGrid()); LatticeCoordinate(t,3);
    if (!hasPhase_)
    {
        Complex           i(0.0,1.0);
        std::vector<Real> p;
        
        envGetTmp(LatticeComplex, coor);
        p  = strToVec<Real>(par().mom);
        //ph = zero;
        ph = Zero();
        for(unsigned int mu = 0; mu < env().getNd(); mu++)
        {
            LatticeCoordinate(coor, mu);
            ph = ph + (p[mu]/env().getDim(mu))*coor;
        }
        ph = exp((Real)(2*M_PI)*i*ph);
        LOG(Message) << "Computing staggered phase." << std::endl;
        // based on phases in Grid/qcd/action/fermion/FermionOperatorImpl.h
        // Staggered Phase does not include (-1)^x from "hermiticity" trans.
        stag_ph = 1.0;
        Lattice<iScalar<vInteger> > x(env().getGrid()); LatticeCoordinate(x,0);
        Lattice<iScalar<vInteger> > y(env().getGrid()); LatticeCoordinate(y,1);
        Lattice<iScalar<vInteger> > z(env().getGrid()); LatticeCoordinate(z,2);
        Lattice<iScalar<vInteger> > t(env().getGrid()); LatticeCoordinate(t,3);
        Lattice<iScalar<vInteger> > lin_x(env().getGrid()); lin_x=y+z+t;
        Lattice<iScalar<vInteger> > lin_y(env().getGrid()); lin_y=x+z+t;
        Lattice<iScalar<vInteger> > lin_z(env().getGrid()); lin_z=x+y+t;
        Lattice<iScalar<vInteger> > lin_5(env().getGrid()); lin_5=x+y+z+t;
        
        // local taste non-singlet ops from Degrand and Detar, Tab. 11.2
        if ( par().gamma == Gamma::Algebra::Gamma5 )      stag_ph = where( mod(lin_5,2)==(Integer)0, stag_ph,-stag_ph);
        else if ( par().gamma == Gamma::Algebra::GammaX ) stag_ph = where( mod(lin_x,2)==(Integer)0, stag_ph,-stag_ph);
        else if ( par().gamma == Gamma::Algebra::GammaY ) stag_ph = where( mod(lin_y,2)==(Integer)0, stag_ph,-stag_ph);
        else if ( par().gamma == Gamma::Algebra::GammaZ ) stag_ph = where( mod(lin_z,2)==(Integer)0, stag_ph,-stag_ph);
        else {
            std::cout << par().gamma << " not implemented for staggered fermon seq. source" << std::endl;
            assert(0);
        }
        
        LatticeCoordinate(t, Tp);
        hasPhase_ = true;
    }
    src = where((t >= par().tA) and (t <= par().tB), ph*stag_ph*q, 0.*q);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_SeqGamma_hpp_
