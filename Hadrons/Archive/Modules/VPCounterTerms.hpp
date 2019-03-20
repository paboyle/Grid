/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Archive/Modules/VPCounterTerms.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: James Harrison <jch1g10@soton.ac.uk>

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
#ifndef Hadrons_MScalar_VPCounterTerms_hpp_
#define Hadrons_MScalar_VPCounterTerms_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         VPCounterTerms                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalar)

class VPCounterTermsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(VPCounterTermsPar,
                                    std::string, source,
                                    double,      mass,
                                    std::string, output,
                                    std::vector<std::string>, outputMom);
};

class TVPCounterTerms: public Module<VPCounterTermsPar>
{
public:
    BASIC_TYPE_ALIASES(SIMPL,);
    class Result: Serializable
    {
    public:
        class Projection: Serializable
        {
        public:
            GRID_SERIALIZABLE_CLASS_MEMBERS(Projection,
                                            std::vector<int>,     momentum,
                                            std::vector<std::vector<std::vector<Complex>>>, twoScalar,
                                            std::vector<std::vector<std::vector<Complex>>>, threeScalar,
                                            std::vector<std::vector<std::vector<Complex>>>, pSquaredInsertion);
        };
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::vector<int>,        lattice_size,
                                        double,                  mass,
                                        std::vector<Projection>, projection);
    };
public:
    // constructor
    TVPCounterTerms(const std::string name);
    // destructor
    virtual ~TVPCounterTerms(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    void project(std::vector<Complex> &projection, const ScalarField &vp, int i_p);
private:
    std::string                freeMomPropName_, GFSrcName_, phatsqName_, prop0Name_,
                               twoscalarName_, twoscalarVertexName_,
                               psquaredName_, psquaredVertexName_;
    std::vector<std::string>   phaseName_, momPhaseName_;
    std::vector<ScalarField *> phase_, momPhase_;
};

MODULE_REGISTER(VPCounterTerms, TVPCounterTerms, MScalar);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalar_VPCounterTerms_hpp_
