/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/Baryon.hpp

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

#ifndef Hadrons_MContraction_Baryon_hpp_
#define Hadrons_MContraction_Baryon_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Grid/qcd/utils/BaryonUtils.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                               Baryon                                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class BaryonPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(BaryonPar,
                                    std::string, q1,
                                    std::string, q2,
                                    std::string, q3,
                                    std::string, gamma,
                                    std::string, output);
};

template <typename FImpl1, typename FImpl2, typename FImpl3>
class TBaryon: public Module<BaryonPar>
{
public:
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
    FERM_TYPE_ALIASES(FImpl3, 3);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::vector<Complex>, corr);
    };
public:
    // constructor
    TBaryon(const std::string name);
    // destructor
    virtual ~TBaryon(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Baryon, ARG(TBaryon<FIMPL, FIMPL, FIMPL>), MContraction);

/******************************************************************************
 *                         TBaryon implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
TBaryon<FImpl1, FImpl2, FImpl3>::TBaryon(const std::string name)
: Module<BaryonPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TBaryon<FImpl1, FImpl2, FImpl3>::getInput(void)
{
    std::vector<std::string> input = {par().q1, par().q2, par().q3};
    
    return input;
}

template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TBaryon<FImpl1, FImpl2, FImpl3>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TBaryon<FImpl1, FImpl2, FImpl3>::setup(void)
{
    envTmpLat(LatticeComplex, "c");
    envTmpLat(LatticeComplex, "diquark");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TBaryon<FImpl1, FImpl2, FImpl3>::execute(void)
{
    LOG(Message) << "Computing nucleon contractions '" << getName() << "' using"
                 << " quarks '" << par().q1 << "', '" << par().q2 << "', and '"
                 << par().q3 << "'" << std::endl;
    
    auto       &q1 = envGet(PropagatorField1, par().q1);
    auto       &q2 = envGet(PropagatorField2, par().q2);
    auto       &q3 = envGet(PropagatorField3, par().q2);
    envGetTmp(LatticeComplex, c);
    envGetTmp(LatticeComplex, diquark);
    Result     result;
    int nt = env().getDim(Tp);
    result.corr.resize(nt);
    const std::string gamma{ par().gamma };
    std::vector<TComplex> buf;
    // C = i gamma_2 gamma_4 => C gamma_5 = - i gamma_1 gamma_3  
/*    Gamma GammaA(Gamma::Algebra::Identity); //Still hardcoded 1
    Gamma GammaB(Gamma::Algebra::SigmaXZ); //Still hardcoded Cg5
    Gamma g4(Gamma::Algebra::GammaT); //needed for parity P_\pm = 0.5*(1 \pm \gamma_4)

    std::vector<std::vector<int>> epsilon = {{0,1,2},{1,2,0},{2,0,1},{0,2,1},{2,1,0},{1,0,2}};
    std::vector<int> epsilon_sgn = {1,1,1,-1,-1,-1};

    char left[] = "uud";
    char right[] = "uud";
    std::vector<int> wick_contraction = {0,0,0,0,0,0};

    for (int ie=0; ie < 6 ; ie++)
      if (left[0] == right[epsilon[ie][0]] && left[1] == right[epsilon[ie][1]] && left[2] == right[epsilon[ie][2]])
        wick_contraction[ie]=1;


    int parity = 1;


    for (int ie_src=0; ie_src < 6 ; ie_src++){
       int a_src = epsilon[ie_src][0]; //a
       int b_src = epsilon[ie_src][1]; //b
       int c_src = epsilon[ie_src][2]; //c
      for (int ie_snk=0; ie_snk < 6 ; ie_snk++){
         int a_snk = epsilon[ie_snk][0]; //a'
         int b_snk = epsilon[ie_snk][1]; //b'
         int c_snk = epsilon[ie_snk][2]; //c'
         auto Daa = peekColour(q2,a_snk,a_src); //D_{alpha' alpha}
         auto Dbb = peekColour(q3,b_snk,b_src); //D_{beta' beta}
         auto Dcc = peekColour(q1,c_snk,c_src); //D_{gamma' gamma}
         auto Dab = peekColour(q2,a_snk,b_src); //D_{alpha' beta}
         auto Dac = peekColour(q2,a_snk,c_src); //D_{alpha' gamma}
         auto Dba = peekColour(q3,b_snk,a_src); //D_{beta' alpha}
         auto Dbc = peekColour(q3,b_snk,c_src); //D_{beta' gamma}
         auto Dca = peekColour(q1,c_snk,a_src); //D_{gamma' alpha}
         auto Dcb = peekColour(q1,c_snk,b_src); //D_{gamma' beta}
         // This is the \delta_{123}^{123} part
         if (wick_contraction[0]){
           diquark = trace(GammaB * Daa * GammaB * Dbb); //1st GammaB and Daa transposed????
           auto temp = GammaA * Dcc * diquark;
           auto g4_temp = GammaA * g4 * temp; 
           c += epsilon_sgn[ie_src] * epsilon_sgn[ie_snk] * 0.5 * trace(GammaA * temp + (double)parity * g4_temp);
         }
         // This is the \delta_{123}^{231} part
         if (wick_contraction[1]){
           auto temp = GammaA * Dca * GammaB * Dab * GammaB * Dbc; //Dab transposed???
           auto g4_temp = GammaA * g4 * temp; 
           c += epsilon_sgn[ie_src] * epsilon_sgn[ie_snk] * 0.5 * trace(GammaA * temp + (double)parity * g4_temp);
         }
         // This is the \delta_{123}^{312} part
         if (wick_contraction[2]){
           auto temp = GammaA * Dcb * GammaB * Dba * GammaB * Dac; //both GammaB and Dba transposed???
           auto g4_temp = GammaA * g4 * temp; 
           c += epsilon_sgn[ie_src] * epsilon_sgn[ie_snk] * 0.5 * trace(GammaA * temp + (double)parity * g4_temp);
         }
         // This is the \delta_{123}^{132} part
         if (wick_contraction[3]){
           diquark = trace(GammaB * Dba * GammaB * Dab); //2nd GammaB and Dab transposed????
           auto temp = GammaA * Dcc * diquark;
           auto g4_temp = GammaA * g4 * temp; 
           c -= epsilon_sgn[ie_src] * epsilon_sgn[ie_snk] * 0.5 * trace(GammaA * temp + (double)parity * g4_temp);
         }
         // This is the \delta_{123}^{321} part
         if (wick_contraction[4]){
           auto temp = GammaA * Dcb * GammaB * Daa * GammaB * Dbc; //1st GammaB and Daa transposed???
           auto g4_temp = GammaA * g4 * temp; 
           c -= epsilon_sgn[ie_src] * epsilon_sgn[ie_snk] * 0.5 * trace(GammaA * temp + (double)parity * g4_temp);
         }
         // This is the \delta_{123}^{213} part
         if (wick_contraction[5]){
           auto temp = GammaA * Dca * GammaB * Dbb * GammaB * Dac; //(Dbb*GammaB) transposed???
           auto g4_temp = GammaA * g4 * temp; 
           c -= epsilon_sgn[ie_src] * epsilon_sgn[ie_snk] * 0.5 * trace(GammaA * temp + (double)parity * g4_temp);
         }
      }
    }
*/

    Gamma GammaA(Gamma::Algebra::Identity);
    Gamma GammaB(Gamma::Algebra::SigmaXZ); //Still hardcoded Cg5
    if (gamma.compare("X") ==0){
      std::cout << "using interpolator C gamma_X";
      Gamma GammaB(Gamma::Algebra::GammaZGamma5); //Still hardcoded CgX = i gamma_3 gamma_5
    } 
    if (gamma.compare("Y") ==0){
      std::cout << "using interpolator C gamma_Y";
      Gamma GammaB(Gamma::Algebra::GammaT); //Still hardcoded CgX = - gamma_4
    } 
    if (gamma.compare("Z")==0){
      std::cout << "using interpolator C gamma_Z";
      Gamma GammaB(Gamma::Algebra::GammaXGamma5); //Still hardcoded CgX = i gamma_1 gamma_5
    } 

    BaryonUtils<FIMPL>::ContractBaryons(q1,q2,q3,GammaA,GammaB,c);

    sliceSum(c,buf,Tp);

    for (unsigned int t = 0; t < buf.size(); ++t)
    {
        result.corr[t] = TensorRemove(buf[t]);
    }

    saveResult(par().output, "baryon", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_Baryon_hpp_
