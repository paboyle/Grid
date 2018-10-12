/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MGauge/Fix.hpp

Copyright (C) 2015
Copyright (C) 2016

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
    virtual ~TGaugeFix(void) = default;
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(GaugeFix, TGaugeFix<GIMPL>, MGauge);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGaugeFix_hpp_
