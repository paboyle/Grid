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

#ifndef Hadrons_MContraction_MesonCCLoop_hpp_
#define Hadrons_MContraction_MesonCCLoop_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MSource/Point.hpp>
#include <Hadrons/Solver.hpp>
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
 *                                TMesonCCLoop                                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)


class MesonCCLoopPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MesonCCLoopPar,
                                    std::string, gauge,
                                    std::string, output,
                                    std::string, solver,
                                    int, inc,
                                    int, tinc);
};

template <typename FImpl1, typename FImpl2>
class TStagMesonCCLoop: public Module<MesonCCLoopPar>
{
public:
    typedef typename FImpl1::FermionField FermionField;
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
    SOLVER_TYPE_ALIASES(FImpl1,);
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
    TStagMesonCCLoop(const std::string name);
    // destructor
    virtual ~TStagMesonCCLoop(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    protected:
    // execution
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    Solver       *solver_{nullptr};
};

MODULE_REGISTER_TMP(StagMesonCCLoop, ARG(TStagMesonCCLoop<STAGIMPL, STAGIMPL>), MContraction);

/******************************************************************************
 *                           TStagMesonCCLoop implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
TStagMesonCCLoop<FImpl1, FImpl2>::TStagMesonCCLoop(const std::string name)
: Module<MesonCCLoopPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
std::vector<std::string> TStagMesonCCLoop<FImpl1, FImpl2>::getInput(void)
{
    std::vector<std::string> input = {par().gauge, par().solver};
    
    return input;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TStagMesonCCLoop<FImpl1, FImpl2>::getOutput(void)
{
    std::vector<std::string> output = {};
    
    return output;
}

// setup ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TStagMesonCCLoop<FImpl1, FImpl2>::setup(void)
{
    envTmpLat(LatticeComplex, "corr");
    envTmpLat(LatticePropagator, "q1");
    envTmpLat(LatticePropagator, "q2");
    envTmpLat(LatticePropagator, "qshift");
    envTmpLat(FermionField, "source");
    envTmpLat(FermionField, "sol");
    
    // grid can't handle real * prop, so use complex
    envTmpLat(LatticeComplex,  "herm_phase");
    envGetTmp(LatticeComplex, herm_phase);
    
    // sink
    Lattice<iScalar<vInteger> > x(env().getGrid()); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(env().getGrid()); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > z(env().getGrid()); LatticeCoordinate(z,2);
    Lattice<iScalar<vInteger> > t(env().getGrid()); LatticeCoordinate(t,3);
    Lattice<iScalar<vInteger> > s(env().getGrid());
    
    //``Hermiticity" phase, (-1)^(x+y)
    // sink only for now
    herm_phase = 1.0;
    s=x+y+z+t;
    herm_phase = where( mod(s,2)==(Integer)0, herm_phase, -herm_phase);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TStagMesonCCLoop<FImpl1, FImpl2>::execute(void)
{
    LOG(Message) << "Computing Conserved Current Stag meson contractions " << std::endl;
    
    //printMem("MesonCCLoop execute() ", env().getGrid()->ThisRank());
    
    std::vector<TComplex>  buf;
    Result    result;
    int nt = env().getDim(Tp);
    int ns = env().getDim(Xp);
    
    result.corr.resize(nt);
    
    auto &U       = envGet(LatticeGaugeField, par().gauge);
    auto &solver  = envGet(Solver, par().solver);
    auto &mat     = solver.getFMat();
    
    envGetTmp(LatticeComplex, corr);
    envGetTmp(LatticeComplex, herm_phase);
    PropagatorField1 q1(U.Grid());
    PropagatorField2 q2(U.Grid());
    PropagatorField1 qshift(U.Grid());
    
    // Do spatial gamma's only
    Lattice<iScalar<vInteger> > x(U.Grid()); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(U.Grid()); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > lin_z(U.Grid()); lin_z=x+y;
    LatticeComplex phases(U.Grid());
    
    LatticeColourMatrix Umu(U.Grid());
    envGetTmp(FermionField, source);
    envGetTmp(FermionField, sol);
    Coordinate srcSite;
    ColourMatrix UmuSrc;
    ColourVector Csrc;
    std::string outFileName;
    
    // loop over source position
    for(int mu=0;mu<3;mu++){
        
        //staggered phases go into links
        Umu = PeekIndex<LorentzIndex>(U,mu);
        phases=1.0;
        if(mu==0){
        }else if(mu==1){
            phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
        }else if(mu==2){
            phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
        }else assert(0);
        Umu *= phases;
        ComplexD src_phase;
        
        for(int t=0; t<nt;t+=par().tinc){
            for(int z=0; z<ns;z+=par().inc){
                for(int y=0; y<ns;y+=par().inc){
                    for(int x=0; x<ns;x+=par().inc){
                        
                        srcSite[0]=x;
                        srcSite[1]=y;
                        srcSite[2]=z;
                        srcSite[3]=t;
                        
                        LOG(Message) << "StagMesonCCLoop src_xyzt " << srcSite[0] <<" "<< srcSite[1]<<" "<<srcSite[2] <<" "<< srcSite[3] <<" mu "<< mu << std::endl;
                        
                        peekSite(UmuSrc, Umu, srcSite);
    
                        for (unsigned int c = 0; c < FImpl1::Dimension; ++c){
                            source = Zero();
                            Csrc=Zero();
                            Csrc()()(c)=Complex(1.0,0.0);
                            pokeSite(Csrc,source,srcSite);
                            sol = Zero();
                            solver(sol, source);
                            FermToProp<FImpl1>(q1, sol, c);
                        }
                        
                        srcSite[mu]=(srcSite[mu]+1)%ns;
                        
                        for (unsigned int c = 0; c < FImpl1::Dimension; ++c){
                            source = Zero();
                            Csrc=Zero();
                            Csrc()()(c)=Complex(1.0,0.0);
                            pokeSite(Csrc,source,srcSite);
                            sol = Zero();
                            solver(sol, source);
                            FermToProp<FImpl1>(q2, sol, c);
                        }
                        
                        qshift = Cshift(q2, mu, 1);
                        corr = trace(adj(qshift) * adj(Umu) * q1 * UmuSrc);
                        corr += trace(adj(q1) * Umu * qshift * adj(UmuSrc));

                        qshift = Cshift(q1, mu, 1);
                        // minus signs from herm transform of prop
                        corr -= trace(adj(q2) * Umu * qshift * UmuSrc); // -1^muhat
                        corr -= trace(adj(qshift) * adj(Umu) * q2 * adj(UmuSrc)); //-1^muhat

                        corr *= herm_phase; // sink
                        if((x+y+z+t)%2)corr *= -1.0; // source
                        
                        sliceSum(corr, buf, Tp);
                        for (unsigned int tsnk = 0; tsnk < buf.size(); ++tsnk){
                            result.corr[tsnk] = TensorRemove(buf[tsnk]);
                        }
                        outFileName = par().output+"_"+
                            std::to_string(x)+"_"+
                            std::to_string(y)+"_"+
                            std::to_string(z)+"_"+
                            std::to_string(t)+"_mu_"+
                            std::to_string(mu);
                        saveResult(outFileName, "mesonCCLoop", result);
                    }
                }
            }
        }
    }
}
END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_MesonCCLoop_hpp_
