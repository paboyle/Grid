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

typedef std::pair<Gamma::Algebra, Gamma::Algebra> GammaAB;
typedef std::pair<GammaAB, GammaAB> GammaABPair;

class BaryonPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(BaryonPar,
                                    std::string, q1,
                                    std::string, q2,
                                    std::string, q3,
                                    std::string, gammas,
                                    std::string, quarks,
                                    std::string, prefactors,
                                    std::string, parity,
                                    std::string, sink,
                                    std::string, output);
};

template <typename FImpl1, typename FImpl2, typename FImpl3>
class TBaryon: public Module<BaryonPar>
{
public:
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
    FERM_TYPE_ALIASES(FImpl3, 3);
    BASIC_TYPE_ALIASES(ScalarImplCR, Scalar);
    SINK_TYPE_ALIASES(Scalar);
    class Metadata: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        Gamma::Algebra, gammaA_left,
                                        Gamma::Algebra, gammaB_left,
                                        Gamma::Algebra, gammaA_right,
                                        Gamma::Algebra, gammaB_right,
                                        std::string, quarks,
                                        std::string, prefactors,
                                        int, parity);
    };
    typedef Correlator<Metadata> Result;
public:
    // constructor
    TBaryon(const std::string name);
    // destructor
    virtual ~TBaryon(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void parseGammaString(std::vector<GammaABPair> &gammaList);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
    // Which gamma algebra was specified
    Gamma::Algebra  al;
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
    std::vector<std::string> input = {par().q1, par().q2, par().q3, par().sink};
    
    return input;
}

template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TBaryon<FImpl1, FImpl2, FImpl3>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

template <typename FImpl1, typename FImpl2, typename FImpl3>
void TBaryon<FImpl1, FImpl2,FImpl3>::parseGammaString(std::vector<GammaABPair> &gammaList)
{
    gammaList.clear();
    
    std::string gammaString = par().gammas;
    //Shorthands for standard baryon operators
    gammaString = regex_replace(gammaString, std::regex("j12"),"(Identity SigmaXZ)");
    gammaString = regex_replace(gammaString, std::regex("j32X"),"(Identity MinusGammaZGamma5)");
    gammaString = regex_replace(gammaString, std::regex("j32Y"),"(Identity GammaT)");
    gammaString = regex_replace(gammaString, std::regex("j32Z"),"(Identity GammaXGamma5)");
    //Shorthands for less common baryon operators
    gammaString = regex_replace(gammaString, std::regex("j12_alt1"),"(Gamma5 MinusSigmaYT)");
    gammaString = regex_replace(gammaString, std::regex("j12_alt2"),"(Identity GammaYGamma5)");
    
    //A single gamma matrix 
    std::regex rex_g("([0-9a-zA-Z]+)");
    //The full string we expect
    std::regex rex("( *\\(( *\\(([0-9a-zA-Z]+) +([0-9a-zA-Z]+) *\\)){2} *\\) *)+");
    std::smatch sm;
    std::regex_match(gammaString, sm, rex);
    assert(sm[0].matched && "invalid gamma structure.");

    auto gamma_begin = std::sregex_iterator(gammaString.begin(), gammaString.end(), rex_g);
    auto gamma_end = std::sregex_iterator();

    int nGamma = std::distance(gamma_begin, gamma_end); 
    //couldn't find out how to count the size in the iterator, other than looping through it...
  /*  int nGamma=0;
    for (std::sregex_iterator i = gamma_begin; i != gamma_end; ++i) {
	nGamma++;
    }
*/   
    gammaList.resize(nGamma/4);
    std::vector<std::string> gS;
    gS.resize(nGamma);
    //even more ugly workarounds here...
    int iG=0;
    for (std::sregex_iterator i = gamma_begin; i != gamma_end; ++i) {
        std::smatch match = *i;                                                 
        gS[iG] = match.str(); 
	iG++;
    }
    for (int i = 0; i < gammaList.size(); i++){
	std::vector<Gamma::Algebra> gS1 = strToVec<Gamma::Algebra>(gS[4*i]);
	std::vector<Gamma::Algebra> gS2 = strToVec<Gamma::Algebra>(gS[4*i+1]);
	std::vector<Gamma::Algebra> gS3 = strToVec<Gamma::Algebra>(gS[4*i+2]);
	std::vector<Gamma::Algebra> gS4 = strToVec<Gamma::Algebra>(gS[4*i+3]);
        gammaList[i].first.first=gS1[0];
        gammaList[i].first.second=gS2[0];
        gammaList[i].second.first=gS3[0];
        gammaList[i].second.second=gS4[0];
    }
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TBaryon<FImpl1, FImpl2, FImpl3>::setup(void)
{
    envTmpLat(LatticeComplex, "c");
    envTmpLat(LatticeComplex, "c2");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TBaryon<FImpl1, FImpl2, FImpl3>::execute(void)
{

    std::vector<std::string> quarks = strToVec<std::string>(par().quarks);    
    std::vector<double> prefactors = strToVec<double>(par().prefactors);    
    int nQ=quarks.size();
    const int  parity {par().parity.size()>0 ? std::stoi(par().parity) : 1};

    std::vector<GammaABPair> gammaList;
    parseGammaString(gammaList);

    assert(prefactors.size()==nQ && "number of prefactors needs to match number of quark-structures.");
    for (int iQ = 0; iQ < nQ; iQ++)
        assert(quarks[iQ].size()==3 && "quark-structures must consist of 3 quarks each.");

    LOG(Message) << "Computing baryon contractions '" << getName() << "'" << std::endl;
    for (int iQ1 = 0; iQ1 < nQ; iQ1++)
        for (int iQ2 = 0; iQ2 < nQ; iQ2++)
            LOG(Message) << prefactors[iQ1]*prefactors[iQ2] << "*<" << quarks[iQ1] << "|" << quarks[iQ2] << ">" << std::endl;
    LOG(Message) << " using quarks " << par().q1 << "', " << par().q2 << "', and '" << par().q3 << std::endl;
     for (int iG = 0; iG < gammaList.size(); iG++)
         LOG(Message) << "' with (Gamma^A,Gamma^B)_left = ( " << gammaList[iG].first.first << " , " << gammaList[iG].first.second << "') and (Gamma^A,Gamma^B)_right = ( " << gammaList[iG].second.first << " , " << gammaList[iG].second.second << ")" << std::endl; 
      LOG(Message) << "and parity " << parity << " using sink " << par().sink << "." << std::endl;
        
    envGetTmp(LatticeComplex, c);
    envGetTmp(LatticeComplex, c2);
    int nt = env().getDim(Tp);
    std::vector<TComplex> buf;
    TComplex cs;
    TComplex ch;

    std::vector<Result> result;
    Result              r;
    r.info.parity = parity;
    r.info.quarks = par().quarks;
    r.info.prefactors = par().prefactors;

    if (envHasType(SlicedPropagator1, par().q1) and
        envHasType(SlicedPropagator2, par().q2) and
        envHasType(SlicedPropagator3, par().q3))
    {
        auto &q1 = envGet(SlicedPropagator1, par().q1);
        auto &q2 = envGet(SlicedPropagator2, par().q2);
        auto &q3 = envGet(SlicedPropagator3, par().q3);
        for (unsigned int i = 0; i < gammaList.size(); ++i)
        {
            r.info.gammaA_left = gammaList[i].first.first;
            r.info.gammaB_left = gammaList[i].first.second;
            r.info.gammaA_right = gammaList[i].second.first;
            r.info.gammaB_right = gammaList[i].second.second;

            Gamma gAl(gammaList[i].first.first);
            Gamma gBl(gammaList[i].first.second);
            Gamma gAr(gammaList[i].second.first);
            Gamma gBr(gammaList[i].second.second);
        
            LOG(Message) << "(propagator already sinked)" << std::endl;
            r.corr.clear();
            for (unsigned int t = 0; t < buf.size(); ++t)
            {
                cs = Zero();
                for (int iQ1 = 0; iQ1 < nQ; iQ1++){
                    for (int iQ2 = 0; iQ2 < nQ; iQ2++){
                        BaryonUtils<FIMPL>::ContractBaryons_Sliced(q1[t],q2[t],q3[t],gAl,gBl,gAr,gBr,quarks[iQ1].c_str(),quarks[iQ2].c_str(),parity,ch);
                        cs += prefactors[iQ1]*prefactors[iQ2]*ch;
                    }
                }
                r.corr.push_back(TensorRemove(cs));
            }
            result.push_back(r);
        }
    }
    else
    {
        auto       &q1 = envGet(PropagatorField1, par().q1);
        auto       &q2 = envGet(PropagatorField2, par().q2);
        auto       &q3 = envGet(PropagatorField3, par().q3);
        for (unsigned int i = 0; i < gammaList.size(); ++i)
        {
            r.info.gammaA_left = gammaList[i].first.first;
            r.info.gammaB_left = gammaList[i].first.second;
            r.info.gammaA_right = gammaList[i].second.first;
            r.info.gammaB_right = gammaList[i].second.second;

            Gamma gAl(gammaList[i].first.first);
            Gamma gBl(gammaList[i].first.second);
            Gamma gAr(gammaList[i].second.first);
            Gamma gBr(gammaList[i].second.second);
        
            std::string ns;
                
            ns = vm().getModuleNamespace(env().getObjectModule(par().sink));
            if (ns == "MSource")
            {
                c=Zero();
                for (int iQ1 = 0; iQ1 < nQ; iQ1++){
                    for (int iQ2 = 0; iQ2 < nQ; iQ2++){
                        BaryonUtils<FIMPL>::ContractBaryons(q1,q2,q3,gAl,gBl,gAr,gBr,quarks[iQ1].c_str(),quarks[iQ2].c_str(),parity,c2);
                        c+=prefactors[iQ1]*prefactors[iQ2]*c2;
                    }
                }
                PropagatorField1 &sink = envGet(PropagatorField1, par().sink);
                auto test = closure(trace(sink*c));     
                sliceSum(test, buf, Tp); 
            }
            else if (ns == "MSink")
            {
                c=Zero();
                for (int iQ1 = 0; iQ1 < nQ; iQ1++){
                    for (int iQ2 = 0; iQ2 < nQ; iQ2++){
                        BaryonUtils<FIMPL>::ContractBaryons(q1,q2,q3,gAl,gBl,gAr,gBr,quarks[iQ1].c_str(),quarks[iQ2].c_str(),parity,c2);
                        c+=prefactors[iQ1]*prefactors[iQ2]*c2;
                    }
                }
                SinkFnScalar &sink = envGet(SinkFnScalar, par().sink);
                buf = sink(c);
            } 
            r.corr.clear();
            for (unsigned int t = 0; t < buf.size(); ++t)
            {
                r.corr.push_back(TensorRemove(buf[t]));
            }
            result.push_back(r);
        }
    }

    saveResult(par().output, "baryon", result);

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_Baryon_hpp_
