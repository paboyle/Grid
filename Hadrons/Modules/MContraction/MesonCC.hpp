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


class MesonPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MesonPar,
                                    std::string, gauge,
                                    std::string, q1,
                                    std::string, q2,
                                    std::string, source,
                                    std::string, output);
};

template <typename FImpl1, typename FImpl2>
class TStagMesonCC: public Module<MesonPar>
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
    std::vector<RealD>               stag_phase_source;
};

MODULE_REGISTER_TMP(StagMesonCC, ARG(TStagMesonCC<STAGIMPL, STAGIMPL>), MContraction);

/******************************************************************************
 *                           TStagMesonCC implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
TStagMesonCC<FImpl1, FImpl2>::TStagMesonCC(const std::string name)
: Module<MesonPar>(name)
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
    envTmpLat(LatticeComplex, "c");
    parseGammaString();
    int Ngam=gammaList.size();

    // grid can't handle real * prop, so use complex
    envTmp(std::vector<LatticeComplex>,  "stag_phase_sink", 1, Ngam,
           LatticeComplex(env().getGrid()));
    envGetTmp(std::vector<LatticeComplex>,stag_phase_sink);
    stag_phase_source.resize(Ngam);
    
    Lattice<iScalar<vInteger> > x(env().getGrid()); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(env().getGrid()); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > z(env().getGrid()); LatticeCoordinate(z,2);
    
    // coordinate of source
    std::vector<int> src_coor = strToVec<int>(static_cast<MSource::StagPoint *>(vm().getModule(par().source))->par().position);
    // local taste non-singlet ops, including ``Hermiticity" phase,
    // see Tab. 11.2 in Degrand and Detar
    for(int i=0; i < gammaList.size(); i++){

        stag_phase_sink[i] = 1.0;
        stag_phase_source[i] = 1.0;
        
        LOG(Message) << "Using gamma: " << gammaList[i] << std::endl;
        switch(gammaList[i]) {
                
            case Gamma::Algebra::GammaX  :
                stag_phase_sink[i] = where( mod(x,2)==(Integer)0, stag_phase_sink[i], -stag_phase_sink[i]);
                if((src_coor[0])%2) stag_phase_source[i]= -stag_phase_source[i];
                break;
                
            case Gamma::Algebra::GammaY  :
                stag_phase_sink[i] = where( mod(y,2)==(Integer)0, stag_phase_sink[i], -stag_phase_sink[i]);
                if((src_coor[1])%2) stag_phase_source[i]= -stag_phase_source[i];
                break;
                
            case Gamma::Algebra::GammaZ  :
                stag_phase_sink[i] = where( mod(z,2)==(Integer)0, stag_phase_sink[i], -stag_phase_sink[i]);
                if((src_coor[2])%2) stag_phase_source[i] = -stag_phase_source[i];
                break;

            case Gamma::Algebra::Gamma5  :
                break;

            default :
                std::cout << "your gamma is not supported for stag meson" << std::endl;
                assert(0);
        }
    }
}

// execution ///////////////////////////////////////////////////////////////////
#define StagMesonCCConnected(q1, q2, gSnk, gSrc) \
(gSnk)*(q1)*adj(q2)*(gSrc)

template <typename FImpl1, typename FImpl2>
void TStagMesonCC<FImpl1, FImpl2>::execute(void)
{
    LOG(Message) << "Computing meson contractions '" << getName() << "' using"
    << " quarks '" << par().q1 << "' and '" << par().q2 << "'"
    << std::endl;
    
    
    std::vector<TComplex>  buf;
    Result    result;
    int                    nt = env().getDim(Tp);
    // staggered gammas
    envGetTmp(std::vector<LatticeComplex>,stag_phase_sink);
    
    result.corr.resize(nt);
    
    auto &U = envGet(LatticeGaugeField, par().gauge);
    auto &q1 = envGet(PropagatorField1, par().q1);
    auto &q2 = envGet(PropagatorField2, par().q2);
        
    envGetTmp(LatticeComplex, c);
    
    // Do spatial gamma's only
    // Staggered Phases.
    //Lattice<iScalar<vInteger> > coor(U.Grid());
    Lattice<iScalar<vInteger> > x(U.Grid()); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(U.Grid()); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > lin_z(U.Grid()); lin_z=x+y;
    ComplexField phases(U.Grid());
    phases=1.0;
    int mu;
    if(gamma_[0]==Gamma::Algebra::GammaX)mu=0;
    else if(gamma_[0]==Gamma::Algebra::GammaY){
        mu=1;
        phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
    } else if(gamma_[0]==Gamma::Algebra::GammaZ){
        mu=2;
        phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
    } else assert(0);
    //if ( mu == 3 ) phases = where( mod(lin_t,2)==(Integer)0, phases,-phases);
    // U_mu(x) right(x+mu)
    LatticeColourMatrix Umu(U.Grid());
    Umu = PeekIndex<LorentzIndex>(U,mu);
    Umu *= phases;
    
    for(int j=0;j<N_j;j++)
        right[j] = Umu*Cshift(right[j], mu, 1);
    
    PropagatorField1 &sink = envGet(PropagatorField1, par().sink);
    c = trace(StagMesonCCConnected(q1, q2, stag_phase_sink[i], stag_phase_source[i]));
    sliceSum(c, buf, Tp);

    for (unsigned int t = 0; t < buf.size(); ++t){
        result.corr[t] = TensorRemove(buf[t]);
    }

    saveResult(par().output, "mesonCC", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_MesonCC_hpp_
