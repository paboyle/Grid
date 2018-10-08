/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MSink/Point.hpp

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Guido Cossu <guido.cossu@ed.ac.uk>
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

#ifndef Hadrons_MSink_Point_hpp_
#define Hadrons_MSink_Point_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                   Point                                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSink)

class PointPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PointPar,
                                    std::string, mom);
};

template <typename FImpl>
class TPoint: public Module<PointPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
    SINK_TYPE_ALIASES();
public:
    // constructor
    TPoint(const std::string name);
    // destructor
    virtual ~TPoint(void) {};
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
    std::string momphName_;
};

MODULE_REGISTER_TMP(Point,       TPoint<FIMPL>,        MSink);
MODULE_REGISTER_TMP(ScalarPoint, TPoint<ScalarImplCR>, MSink);

/******************************************************************************
 *                          TPoint implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TPoint<FImpl>::TPoint(const std::string name)
: Module<PointPar>(name)
, momphName_ (name + "_momph")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TPoint<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TPoint<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPoint<FImpl>::setup(void)
{
    envTmpLat(LatticeComplex, "coor");
    envCacheLat(LatticeComplex, momphName_);
    envCreate(SinkFn, getName(), 1, nullptr);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPoint<FImpl>::execute(void)
{   
    LOG(Message) << "Setting up point sink function for momentum ["
                 << par().mom << "]" << std::endl;

    auto &ph = envGet(LatticeComplex, momphName_);
    
    if (!hasPhase_)
    {
        Complex           i(0.0,1.0);
        std::vector<Real> p;

        envGetTmp(LatticeComplex, coor);
        p  = strToVec<Real>(par().mom);
        ph = zero;
        for(unsigned int mu = 0; mu < p.size(); mu++)
        {
            LatticeCoordinate(coor, mu);
            ph = ph + (p[mu]/env().getDim(mu))*coor;
        }
        ph = exp((Real)(2*M_PI)*i*ph);
        hasPhase_ = true;
    }
    auto sink = [&ph](const PropagatorField &field)
    {
        SlicedPropagator res;
        PropagatorField  tmp = ph*field;
        
        sliceSum(tmp, res, Tp);
        
        return res;
    };
    envGet(SinkFn, getName()) = sink;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSink_Point_hpp_
