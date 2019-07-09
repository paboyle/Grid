/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/Nucleon.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Felix Erben <felix.erben@ed.ac.uk>

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

#ifndef Hadrons_MContraction_Nucleon_hpp_
#define Hadrons_MContraction_Nucleon_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                               Nucleon                                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class NucleonPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(NucleonPar,
                                    std::string, q1,
                                    std::string, q2,
                                    std::string, q3,
                                    std::string, output);
};

template <typename FImpl1, typename FImpl2, typename FImpl3>
class TNucleon: public Module<NucleonPar>
{
public:
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
    FERM_TYPE_ALIASES(FImpl3, 3);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::vector<std::vector<std::vector<Complex>>>, corr);
    };
public:
    // constructor
    TNucleon(const std::string name);
    // destructor
    virtual ~TNucleon(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Nucleon, ARG(TNucleon<FIMPL, FIMPL, FIMPL>), MContraction);

/******************************************************************************
 *                         TNucleon implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
TNucleon<FImpl1, FImpl2, FImpl3>::TNucleon(const std::string name)
: Module<NucleonPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TNucleon<FImpl1, FImpl2, FImpl3>::getInput(void)
{
    std::vector<std::string> input = {par().q1, par().q2, par().q3};
    
    return input;
}

template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TNucleon<FImpl1, FImpl2, FImpl3>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TNucleon<FImpl1, FImpl2, FImpl3>::setup(void)
{
    envTmpLat(LatticeComplex, "c");
    envTmpLat(LatticeComplex, "diquark");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TNucleon<FImpl1, FImpl2, FImpl3>::execute(void)
{
    LOG(Message) << "Computing nucleon contractions '" << getName() << "' using"
                 << " quarks '" << par().q1 << "', '" << par().q2 << "', and '"
                 << par().q3 << "'" << std::endl;
    
    auto       &q1 = envGet(PropagatorField1, par().q1);
    auto       &q2 = envGet(PropagatorField2, par().q2);
    auto       &q3 = envGet(PropagatorField3, par().q2);
    envGetTmp(LatticeComplex, c);
    //envGetTmp(LatticeComplex, quark2);
    //envGetTmp(LatticeComplex, quark3);
    envGetTmp(LatticeComplex, diquark);
    Result     result;
   
    // C = i gamma_2 gamma_4 => C gamma_5 = - i gamma_1 gamma_3  
    Gamma Cg5(Gamma::Algebra::SigmaXZ);
    Gamma g4(Gamma::Algebra::GammaT); //needed for parity P_\pm = 0.5*(1 \pm \gamma_4)

    std::vector<std::vector<int>> epsilon = {{0,1,2},{1,2,0},{2,0,1},{0,2,1},{2,1,0},{1,0,2}};
    std::vector<int> epsilon_sgn = {1,1,1,-1,-1,-1};
    // This is the \delta_{123}^{123} part
    for (int ie_src=0; ie_src < 6 ; ie_src++){
       int c1_src = epsilon[ie_src][0];
       int c2_src = epsilon[ie_src][1];
       int c3_src = epsilon[ie_src][2];
      for (int ie_snk=0; ie_snk < 6 ; ie_snk++){
         int c1_snk = epsilon[ie_snk][0];
         int c2_snk = epsilon[ie_snk][1];
         int c3_snk = epsilon[ie_snk][2];
         auto Dcc = peekColour(q1,c1_snk,c1_src); //D_{gamma' gamma}
         auto Daa = peekColour(q2,c2_snk,c2_src); //D_{alpha' alpha}
         auto Dbb = peekColour(q3,c3_snk,c3_src); //D_{beta' beta}
         diquark = trace(Cg5 * Daa * Cg5 * Dbb); //Daa transposed????
           //diquark = q2()()(c2,1) * Cg5 * q3()()(c3,2); //Why does this not work??
         auto temp = Dcc * diquark;
         auto g4_temp = g4 * temp; 
         int parity = 1;
         c += epsilon_sgn[ie_src] * epsilon_sgn[ie_snk] * 0.5 * (double)parity * trace(temp + g4_temp);
      }
    }


    // This is the \delta_{123}^{213} part
    for (int ie_src=0; ie_src < 6 ; ie_src++){
       int c1_src = epsilon[ie_src][0];
       int c2_src = epsilon[ie_src][1];
       int c3_src = epsilon[ie_src][2];
      for (int ie_snk=0; ie_snk < 6 ; ie_snk++){
         int c1_snk = epsilon[ie_snk][0];
         int c2_snk = epsilon[ie_snk][1];
         int c3_snk = epsilon[ie_snk][2];
         auto Dca = peekColour(q1,c1_snk,c2_src); //D_{gamma' alpha}
         auto Dac = peekColour(q2,c2_snk,c1_src); //D_{alpha' gamma}
         auto Dbb = peekColour(q3,c3_snk,c3_src); //D_{beta' beta}
         auto temp = Dca * Cg5 * Dbb * Cg5 * Dac; //(Dbb*Cg5) transposed???
         auto g4_temp = g4 * temp; 
         int parity = 1;
         c -= epsilon_sgn[ie_src] * epsilon_sgn[ie_snk] * 0.5 * (double)parity * trace(temp + g4_temp);
      }
    }

    // saveResult(par().output, "meson", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_Nucleon_hpp_
