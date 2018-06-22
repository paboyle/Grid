/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MFermion/GaugeProp.hpp

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Guido Cossu <guido.cossu@ed.ac.uk>
Author: Lanny91 <andrew.lawson@gmail.com>
Author: pretidav <david.preti@csic.es>

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

#ifndef Hadrons_MFermion_GaugeProp_hpp_
#define Hadrons_MFermion_GaugeProp_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                GaugeProp                                   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MFermion)

class GaugePropPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GaugePropPar,
                                    std::string, source,
                                    std::string, solver);
};

template <typename FImpl>
class TGaugeProp: public Module<GaugePropPar>
{
public:
    FG_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TGaugeProp(const std::string name);
    // destructor
    virtual ~TGaugeProp(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Ls_;
    Solver       *solver_{nullptr};
};

MODULE_REGISTER_TMP(GaugeProp, TGaugeProp<FIMPL>, MFermion);
MODULE_REGISTER_TMP(ZGaugeProp, TGaugeProp<ZFIMPL>, MFermion);

/******************************************************************************
 *                      TGaugeProp implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TGaugeProp<FImpl>::TGaugeProp(const std::string name)
: Module<GaugePropPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TGaugeProp<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().source, par().solver};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TGaugeProp<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_5d"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TGaugeProp<FImpl>::setup(void)
{
    Ls_ = env().getObjectLs(par().solver);
    envCreateLat(PropagatorField, getName());
    envTmpLat(FermionField, "source", Ls_);
    envTmpLat(FermionField, "sol", Ls_);
    envTmpLat(FermionField, "tmp");
    if (Ls_ > 1)
    {
        envCreateLat(PropagatorField, getName() + "_5d", Ls_);
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TGaugeProp<FImpl>::execute(void)
{
    LOG(Message) << "Computing quark propagator '" << getName() << "'"
                 << std::endl;
    
    std::string propName = (Ls_ == 1) ? getName() : (getName() + "_5d");
    auto        &prop    = envGet(PropagatorField, propName);
    auto        &fullSrc = envGet(PropagatorField, par().source);
    auto        &solver  = envGet(Solver, par().solver);
    auto        &mat     = solver.getFMat();
    
    envGetTmp(FermionField, source);
    envGetTmp(FermionField, sol);
    envGetTmp(FermionField, tmp);
    LOG(Message) << "Inverting using solver '" << par().solver
                 << "' on source '" << par().source << "'" << std::endl;
    for (unsigned int s = 0; s < Ns; ++s)
    for (unsigned int c = 0; c < FImpl::Dimension; ++c)
    {
        LOG(Message) << "Inversion for spin= " << s << ", color= " << c
                     << std::endl;
        // source conversion for 4D sources
        LOG(Message) << "Import source" << std::endl;
        if (!env().isObject5d(par().source))
        {
            if (Ls_ == 1)
            {
               PropToFerm<FImpl>(source, fullSrc, s, c);
            }
            else
            {
                PropToFerm<FImpl>(tmp, fullSrc, s, c);
                mat.ImportPhysicalFermionSource(tmp, source);
            }
        }
        // source conversion for 5D sources
        else
        {
            if (Ls_ != env().getObjectLs(par().source))
            {
                HADRONS_ERROR(Size, "Ls mismatch between quark action and source");
            }
            else
            {
                PropToFerm<FImpl>(source, fullSrc, s, c);
            }
        }
        LOG(Message) << "Solve" << std::endl;
        sol = zero;
        solver(sol, source);
        LOG(Message) << "Export solution" << std::endl;
        FermToProp<FImpl>(prop, sol, s, c);
        // create 4D propagators from 5D one if necessary
        if (Ls_ > 1)
        {
            PropagatorField &p4d = envGet(PropagatorField, getName());
            mat.ExportPhysicalFermionSolution(sol, tmp);
            FermToProp<FImpl>(p4d, tmp, s, c);
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_GaugeProp_hpp_
