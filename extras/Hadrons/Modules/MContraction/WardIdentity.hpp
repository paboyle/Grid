/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MContraction/WardIdentity.hpp

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

#ifndef Hadrons_WardIdentity_hpp_
#define Hadrons_WardIdentity_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
  Ward Identity contractions
 -----------------------------
 
 * options:
 - q:          propagator, 5D if available (string)
 - action:     action module used for propagator solution (string)
 - mass:       mass of quark (double)
 - test_axial: whether or not to test PCAC relation.
*/

/******************************************************************************
 *                              WardIdentity                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class WardIdentityPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WardIdentityPar,
                                    std::string, q,
                                    std::string, action,
                                    double,      mass,
                                    bool,        test_axial);
};

template <typename FImpl>
class TWardIdentity: public Module<WardIdentityPar>
{
public:
    TYPE_ALIASES(FImpl,);
public:
    // constructor
    TWardIdentity(const std::string name);
    // destructor
    virtual ~TWardIdentity(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Ls_;
};

MODULE_REGISTER_NS(WardIdentity, TWardIdentity<FIMPL>, MContraction);

/******************************************************************************
 *                     TWardIdentity implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TWardIdentity<FImpl>::TWardIdentity(const std::string name)
: Module<WardIdentityPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TWardIdentity<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().action};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TWardIdentity<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWardIdentity<FImpl>::setup(void)
{
    Ls_ = env().getObjectLs(par().q);
    if (Ls_ != env().getObjectLs(par().action))
    {
        HADRON_ERROR("Ls mismatch between quark action and propagator");
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWardIdentity<FImpl>::execute(void)
{
    LOG(Message) << "Performing Ward Identity checks for quark '" << par().q
                 << "'." << std::endl;

    PropagatorField psi(env().getGrid()), tmp(env().getGrid()),
                    vector_WI(env().getGrid());
    PropagatorField &q    = *env().template getObject<PropagatorField>(par().q);
    FMat            &act = *(env().template getObject<FMat>(par().action));
    Gamma           g5(Gamma::Algebra::Gamma5);
    
    // Compute D_mu V_mu, D here is backward derivative.
    vector_WI    = zero;
    for (unsigned int mu = 0; mu < Nd; ++mu)
    {
        act.ContractConservedCurrent(q, q, tmp, Current::Vector, mu);
        tmp -= Cshift(tmp, mu, -1);
        vector_WI += tmp;
    }

    LOG(Message) << "Vector Ward Identity check Delta_mu V_mu = " 
                 << norm2(vector_WI) << std::endl;

    if (par().test_axial)
    {
        LatticeComplex PP(env().getGrid()), axial_defect(env().getGrid()),
                       PJ5q(env().getGrid());

        // Compute D_mu A_mu, D is backwards derivative.
        axial_defect = zero;
        for (unsigned int mu = 0; mu < Nd; ++mu)
        {
            act.ContractConservedCurrent(q, q, tmp, Current::Axial, mu);
            tmp -= Cshift(tmp, mu, -1);
            axial_defect += trace(g5*tmp);
        }

        // Get PJ5q for 5D (zero for 4D) and PP.
        PJ5q = zero;
        if (Ls_ > 1)
        {
            // PP
            ExtractSlice(tmp, q, 0, 0);
            psi  = (tmp - g5*tmp);
            ExtractSlice(tmp, q, Ls_ - 1, 0);
            psi += (tmp + g5*tmp);
            PP = trace(adj(psi)*psi);

            // P5Jq
            ExtractSlice(tmp, q, Ls_/2 - 1, 0);
            psi  = 0.5 * (tmp + g5*tmp);
            ExtractSlice(tmp, q, Ls_/2, 0);
            psi += 0.5 * (tmp - g5*tmp);
            PJ5q = trace(adj(psi)*psi);
        }
        else
        {
            PP = trace(adj(q)*q);
        }

        // Test ward identities, D_mu V_mu = 0; D_mu A_mu = 2m<PP> + 2 PJ5q
        axial_defect -= 2.*PJ5q;
        axial_defect -= 2.*(par().mass)*PP;
        LOG(Message) << "Axial Ward Identity defect Delta_mu A_mu = "
                     << norm2(axial_defect) << std::endl;
        LOG(Message) << "norm2(PP) = " << norm2(PP) << std::endl; 
        LOG(Message) << "norm2(PJ5q) = " << norm2(PJ5q) << std::endl; 
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_WardIdentity_hpp_
