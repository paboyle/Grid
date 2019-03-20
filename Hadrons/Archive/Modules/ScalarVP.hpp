/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Archive/Modules/ScalarVP.hpp

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
#ifndef Hadrons_MScalar_ScalarVP_hpp_
#define Hadrons_MScalar_ScalarVP_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Scalar vacuum polarisation                         *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalar)

class ScalarVPPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ScalarVPPar,
                                    std::string, emField,
                                    std::string, scalarProp,
                                    std::string, output,
                                    std::vector<std::string>, outputMom);
};

class TScalarVP: public Module<ScalarVPPar>
{
public:
    BASIC_TYPE_ALIASES(SIMPL,);
    typedef PhotonR::GaugeField     EmField;
    typedef PhotonR::GaugeLinkField EmComp;
    class Result: Serializable
    {
    public:
        class Projection: Serializable
        {
        public:
            GRID_SERIALIZABLE_CLASS_MEMBERS(Projection,
                                            std::vector<int>,     momentum,
                                            std::vector<std::vector<std::vector<Complex>>>, pi,
                                            std::vector<std::vector<std::vector<Complex>>>, pi_free,
                                            std::vector<std::vector<std::vector<Complex>>>, pi_2E,
                                            std::vector<std::vector<std::vector<Complex>>>, pi_2T,
                                            std::vector<std::vector<std::vector<Complex>>>, pi_S,
                                            std::vector<std::vector<std::vector<Complex>>>, pi_4C,
                                            std::vector<std::vector<std::vector<Complex>>>, pi_X,
                                            std::vector<std::vector<std::vector<Complex>>>, pi_srcT,
                                            std::vector<std::vector<std::vector<Complex>>>, pi_snkT);
        };
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::vector<int>,        lattice_size,
                                        double,                  mass,
                                        double,                  charge,
                                        std::vector<Projection>, projection);
    };
public:
    // constructor
    TScalarVP(const std::string name);
    // destructor
    virtual ~TScalarVP(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    void makeCaches(void);
    // conserved vector two-point contraction
    void vpContraction(ScalarField &vp,
                       ScalarField &prop_0_x, ScalarField &prop_nu_x,
                       TComplex u_src, ScalarField &u_snk, int mu);
    // conserved vector two-point contraction with unit gauge link at sink
    void vpContraction(ScalarField &vp,
                       ScalarField &prop_0_x, ScalarField &prop_nu_x,
                       TComplex u_src, int mu);
    // write momentum-projected vacuum polarisation to file(s)
    void project(std::vector<Complex> &projection, const ScalarField &vp,
                 int i_p);
    // momentum-space Delta_1 insertion
    void momD1(ScalarField &s, FFT &fft);
private:
    bool                                        momPhasesDone_;
    std::string                                 freeMomPropName_, GFSrcName_,
                                                prop0Name_, propQName_,
                                                propSunName_, propTadName_,
                                                fftName_;
    std::vector<std::string>                    phaseName_, muPropQName_,
                                                momPhaseName_;
    std::vector<std::vector<std::string> >      vpTensorName_;
    std::vector<ScalarField *>                  phase_, momPhase_;
};

MODULE_REGISTER(ScalarVP, TScalarVP, MScalar);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalar_ScalarVP_hpp_
