/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MGauge/GaugeFix.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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

#ifndef Hadrons_MGaugeFix_hpp_
#define Hadrons_MGaugeFix_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Grid/qcd/utils/GaugeFix.h>
BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                              Fix gauge                                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class GaugeFixPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GaugeFixPar,
                                    std::string, gauge,
                                    Real,  alpha,
                                    int, maxiter, 
                                    Real, Omega_tol, 
                                    Real, Phi_tol,
                                    bool, Fourier);
};

template <typename GImpl>
class TGaugeFix: public Module<GaugeFixPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TGaugeFix(const std::string name);
    // destructor
    virtual ~TGaugeFix(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(GaugeFix, TGaugeFix<GIMPL>, MGauge);

/******************************************************************************
*                            TGaugeFix implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TGaugeFix<GImpl>::TGaugeFix(const std::string name)
: Module<GaugeFixPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TGaugeFix<GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    return in;
}

template <typename GImpl>
std::vector<std::string> TGaugeFix<GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TGaugeFix<GImpl>::setup(void)
{
    envCreateLat(GaugeField, getName());
}


// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TGaugeFix<GImpl>::execute(void)
//Loads the gauge and fixes it
{
    std::cout << "executing" << std::endl;
    LOG(Message) << "Fixing the Gauge" << std::endl;
    LOG(Message) << par().gauge << std::endl;
    auto &U     = envGet(GaugeField, par().gauge);
    auto &Umu   = envGet(GaugeField, getName());
    LOG(Message) << "Gauge Field fetched" << std::endl;
    //do we allow maxiter etc to be user set?
    Real alpha     = par().alpha;
    int  maxiter   = par().maxiter;
    Real Omega_tol = par().Omega_tol;
    Real Phi_tol   = par().Phi_tol;
    bool Fourier   = par().Fourier;
    FourierAcceleratedGaugeFixer<PeriodicGimplR>::SteepestDescentGaugeFix(U,alpha,maxiter,Omega_tol,Phi_tol,Fourier);
    Umu = U;
    LOG(Message) << "Gauge Fixed" << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGaugeFix_hpp_
