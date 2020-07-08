/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/Meson.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Lanny91 <andrew.lawson@gmail.com>
Author: Vera Guelpers <vmg1n14@soton.ac.uk>
Author: Tom Blum (conserved currents)

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

#ifndef Hadrons_MContraction_MesonCC_hpp_
#define Hadrons_MContraction_MesonCC_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MSource/Point.hpp>
//#include <Hadrons/utils_memory.h>

BEGIN_HADRONS_NAMESPACE

/*
 
 2pt conserved-conserved staggered correlation function
 -----------------------------
 
 * options:
 - q1: input propagator 1 (string) src at y
 - q2: input propagator 2 (string) src at y+hat mu
 
*/

/******************************************************************************
 *                                TMesonCC                                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)


class MesonCCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MesonCCPar,
                                    std::string, gauge,
                                    std::string, q1,
                                    std::string, q2,
                                    std::string, source,
                                    std::string, gammas,
                                    std::string, output);
};

template <typename FImpl1, typename FImpl2>
class TStagMesonCC: public Module<MesonCCPar>
{
public:
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
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
    TStagMesonCC(const std::string name);
    // destructor
    virtual ~TStagMesonCC(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    protected:
    // execution
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
};

MODULE_REGISTER_TMP(StagMesonCC, ARG(TStagMesonCC<STAGIMPL, STAGIMPL>), MContraction);

/******************************************************************************
 *                           TStagMesonCC implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
TStagMesonCC<FImpl1, FImpl2>::TStagMesonCC(const std::string name)
: Module<MesonCCPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
std::vector<std::string> TStagMesonCC<FImpl1, FImpl2>::getInput(void)
{
    std::vector<std::string> input = {par().gauge, par().q1, par().q2, par().source};
    
    return input;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TStagMesonCC<FImpl1, FImpl2>::getOutput(void)
{
    std::vector<std::string> output = {};
    
    return output;
}

// setup ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TStagMesonCC<FImpl1, FImpl2>::setup(void)
{
    envTmpLat(LatticeComplex, "corr");
    envTmpLat(LatticePropagator, "qshift");

    // grid can't handle real * prop, so use complex
    envTmpLat(LatticeComplex,  "herm_phase");
    envGetTmp(LatticeComplex, herm_phase);
    
    // sink
    Lattice<iScalar<vInteger> > x(env().getGrid()); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(env().getGrid()); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > z(env().getGrid()); LatticeCoordinate(z,2);
    Lattice<iScalar<vInteger> > t(env().getGrid()); LatticeCoordinate(t,3);
    Lattice<iScalar<vInteger> > s(env().getGrid());
    
    // coordinate of source
    std::vector<int> src_coor = strToVec<int>(static_cast<MSource::StagPoint *>(vm().getModule(par().source))->par().position);
    
    //``Hermiticity" phase, (-1)^(x+y)
    herm_phase = 1.0;
    s=x+y+z+t;
    herm_phase = where( mod(s,2)==(Integer)0, herm_phase, -herm_phase);
    if((src_coor[0]+src_coor[1]+src_coor[2]+src_coor[3])%2)
        herm_phase = -herm_phase;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TStagMesonCC<FImpl1, FImpl2>::execute(void)
{
    LOG(Message) << "Computing meson contractions '" << getName() << "' using"
    << " quarks '" << par().q1 << "' and '" << par().q2 << "'"
    << std::endl;
    
    //printMem("MesonCC execute() ", env().getGrid()->ThisRank());
    
    std::vector<TComplex>  buf;
    Result    result;
    int                    nt = env().getDim(Tp);
    
    result.corr.resize(nt);
    
    auto &U = envGet(LatticeGaugeField, par().gauge);
    auto &q1 = envGet(PropagatorField1, par().q1);
    auto &q2 = envGet(PropagatorField2, par().q2);
        
    envGetTmp(LatticeComplex, corr);
    envGetTmp(LatticeComplex, herm_phase);
    PropagatorField1 qshift(U.Grid());
    FermionField1 tmp(U.Grid());
    
    // Do spatial gamma's only
    // Staggered Phases, put into links
    //Lattice<iScalar<vInteger> > coor(U.Grid());
    Lattice<iScalar<vInteger> > x(U.Grid()); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(U.Grid()); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > lin_z(U.Grid()); lin_z=x+y;
    LatticeComplex phases(U.Grid());
    phases=1.0;
    ComplexD src_phase=1.0;
    
    Coordinate srcSite = strToVec<int>(static_cast<MSource::StagPoint *>(vm().getModule(par().source))->par().position);
    
    int mu;
    std::vector<Gamma::Algebra>        gamma_;
    gamma_ = strToVec<Gamma::Algebra>(par().gammas);
    if(gamma_[0]==Gamma::Algebra::GammaX){
        mu=0;
    }else if(gamma_[0]==Gamma::Algebra::GammaY){
        mu=1;
        phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
        if(srcSite[0]%2) src_phase = -src_phase;
    } else if(gamma_[0]==Gamma::Algebra::GammaZ){
        mu=2;
        phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
        if((srcSite[0]+srcSite[1])%2) src_phase = -src_phase;
    } else assert(0);

    LatticeColourMatrix Umu(U.Grid());
    Umu = PeekIndex<LorentzIndex>(U,mu);
    Umu *= phases;
    
    ColourMatrix UmuSrc;
    peekSite(UmuSrc, Umu, srcSite);
    
    LOG(Message) << "StagMesonCC src_xyzt " << srcSite << " mu " << mu << std::endl;
    
    qshift = Cshift(q2, mu, 1);
    corr = trace(adj(qshift) * adj(Umu) * q1 * UmuSrc);
    corr += trace(adj(q1) * Umu * qshift * adj(UmuSrc));
 
    qshift = Cshift(q1, mu, 1);
    corr -= trace(adj(q2) * Umu * qshift * UmuSrc); // -1^muhat
    corr -= trace(adj(qshift) * adj(Umu) * q2 * adj(UmuSrc)); //-1^muhat
    
    corr *= herm_phase;
    
    sliceSum(corr, buf, Tp);
    for (unsigned int t = 0; t < buf.size(); ++t){
        result.corr[t] = TensorRemove(buf[t]);
    }
    saveResult(par().output, "mesonCC", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_MesonCC_hpp_
