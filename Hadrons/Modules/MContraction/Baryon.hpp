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
                                    std::string, q1_src,
                                    std::string, q2_src,
                                    std::string, q3_src,
                                    std::string, GammaA,
                                    std::string, GammaB,
                                    std::string, quarks_snk,
                                    std::string, quarks_src,
                                    int, parity,
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
    std::vector<std::string> input = {par().q1_src, par().q2_src, par().q3_src, par().sink};
    
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
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TBaryon<FImpl1, FImpl2, FImpl3>::execute(void)
{
    LOG(Message) << "Computing baryon contraction '" << getName() << "' < " << par().quarks_snk << " | " << par().quarks_src << " > using"
                 << " quarks '" << par().q1_src << "', and a diquark formed of ('" << par().q2_src << "', and '"
                 << par().q3_src << "') at the source and (Gamma^A,Gamma^B) = ( " << par().GammaA << " , " << par().GammaB 
                 << " ) and parity " << par().parity << "." << std::endl;
    
    auto       &q1_src = envGet(PropagatorField1, par().q1_src);
    auto       &q2_src = envGet(PropagatorField2, par().q2_src);
    auto       &q3_src = envGet(PropagatorField3, par().q3_src);
    envGetTmp(LatticeComplex, c);
    Result     result;
    int nt = env().getDim(Tp);
    result.corr.resize(nt);
    std::vector<Gamma::Algebra> ggA = strToVec<Gamma::Algebra>(par().GammaA);
    Gamma GammaA(ggA[0]);
    std::vector<Gamma::Algebra> ggB = strToVec<Gamma::Algebra>(par().GammaB);
    Gamma GammaB(ggB[0]);
    std::vector<TComplex> buf;
    vTComplex cs;
    const int  parity {par().parity};
    const char * quarks_snk{par().quarks_snk.c_str()};
    const char * quarks_src{par().quarks_src.c_str()};

    if (envHasType(SlicedPropagator1, par().q1_src) and
        envHasType(SlicedPropagator2, par().q2_src) and
        envHasType(SlicedPropagator3, par().q3_src))
    {
        auto &q1_src = envGet(SlicedPropagator1, par().q1_src);
        auto &q2_src = envGet(SlicedPropagator2, par().q2_src);
        auto &q3_src = envGet(SlicedPropagator3, par().q3_src);
        
        LOG(Message) << "(propagator already sinked)" << std::endl;
        for (unsigned int t = 0; t < buf.size(); ++t)
        {
            //TODO: Get this to compile without the casts. Templates? 
            //BaryonUtils<FIMPL>::ContractBaryons_Sliced(*reinterpret_cast<Grid::iScalar<Grid::iMatrix<Grid::iMatrix<Grid::vComplex, 3>, 4>>*>(&q1_src[t]),*reinterpret_cast<Grid::iScalar<Grid::iMatrix<Grid::iMatrix<Grid::vComplex, 3>, 4>>*>(&q2_src[t]),*reinterpret_cast<Grid::iScalar<Grid::iMatrix<Grid::iMatrix<Grid::vComplex, 3>, 4>>*>(&q3_src[t]),GammaA,GammaB,quarks_snk,quarks_src,parity,cs);
            //result.corr[t] = TensorRemove(*reinterpret_cast<Grid::TComplex*>(&cs));
      //      BaryonUtils<FIMPL>::ContractBaryons_Sliced(q1_src[t],q2_src[t],q3_src[t],GammaA,GammaB,quarks_snk,quarks_src,parity,cs);
    //        result.corr[t] = TensorRemove(cs);
        }
    }
    else
    {
        std::string ns;
                
        ns = vm().getModuleNamespace(env().getObjectModule(par().sink));
        if (ns == "MSource")
        {
         //TODO: Understand what this is and then get it to compile. Hopefully no new function needed. The following lines are from the Meson.hpp module.
        /*    PropagatorField1 &sink = envGet(PropagatorField1, par().sink);
                
            c = trace(mesonConnected(q1, q2, gSnk, gSrc)*sink);
            sliceSum(c, buf, Tp); */
// My attempt at some code, which doesn't work. I also don't know whether anything like this is what we want here.
           /* BaryonUtils<FIMPL>::ContractBaryons(q1_src,q2_src,q3_src,GammaA,GammaB,quarks_snk,quarks_src,parity,c);
            PropagatorField1 &sink = envGet(PropagatorField1, par().sink);
            auto test = trace(c*sink);     
            sliceSum(test, buf, Tp); */
        }
        else if (ns == "MSink")
        {
            BaryonUtils<FIMPL>::ContractBaryons(q1_src,q2_src,q3_src,GammaA,GammaB,quarks_snk,quarks_src,parity,c);

            SinkFnScalar &sink = envGet(SinkFnScalar, par().sink);
            buf = sink(c);
        } 
        for (unsigned int t = 0; t < buf.size(); ++t)
        {
            result.corr[t] = TensorRemove(buf[t]);
        }
    }

    saveResult(par().output, "baryon", result);

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_Baryon_hpp_
