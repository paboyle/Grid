/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MScalarSUN/StochFreeField.hpp

Copyright (C) 2015-2019

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
#ifndef Hadrons_MScalarSUN_StochFreeField_hpp_
#define Hadrons_MScalarSUN_StochFreeField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                      stochastic free SU(N) scalar field                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class StochFreeFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StochFreeFieldPar,
                                    double, m2,
                                    double, g,
                                    double, smearing);
};

template <typename SImpl>
class TStochFreeField: public Module<StochFreeFieldPar>
{
public:
    typedef typename SImpl::Field                    Field;
    typedef typename SImpl::ComplexField             ComplexField;
    typedef typename SImpl::Group                    Group;
    typedef typename SImpl::SiteField::scalar_object Site;
public:
    // constructor
    TStochFreeField(const std::string name);
    // destructor
    virtual ~TStochFreeField(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool create_weight;
};

MODULE_REGISTER_TMP(StochFreeFieldSU2, TStochFreeField<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_TMP(StochFreeFieldSU3, TStochFreeField<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_TMP(StochFreeFieldSU4, TStochFreeField<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_TMP(StochFreeFieldSU5, TStochFreeField<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_TMP(StochFreeFieldSU6, TStochFreeField<ScalarNxNAdjImplR<6>>, MScalarSUN);

/******************************************************************************
 *                 TStochFreeField implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TStochFreeField<SImpl>::TStochFreeField(const std::string name)
: Module<StochFreeFieldPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TStochFreeField<SImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename SImpl>
std::vector<std::string> TStochFreeField<SImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TStochFreeField<SImpl>::setup(void)
{
    create_weight = false; 
    if (!env().hasCreatedObject("_" + getName() + "_weight"))
    {
        envCacheLat(ComplexField, "_" + getName() + "_weight");
        envTmpLat(ComplexField, "smear");
        create_weight = true;
    }
    envTmpLat(Field, "phift");
    envTmpLat(ComplexField, "ca");
    envCreateLat(Field, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TStochFreeField<SImpl>::execute(void)
{
    LOG(Message) << "Generating stochastic scalar field" << std::endl;
    
    const   unsigned int N    = Group::Dimension;
    const   unsigned int Nadj = Group::AdjointDimension;
    auto    &phi              = envGet(Field, getName());
    auto    &w                = envGet(ComplexField, "_" + getName() + "_weight");
    auto    &rng              = rng4d();
    double  trphi2;
    FFT     fft(envGetGrid(Field));
    Integer vol;

    vol = 1;
    for(int d = 0; d < env().getNd(); d++)
    {
        vol = vol*env().getDim(d);
    }
    if (create_weight)
    {
        LOG(Message) << "Caching momentum-space scalar action" << std::endl;
        
        envGetTmp(ComplexField, smear);
        SImpl::MomentaSquare(smear);
        smear = exp(-par().smearing*smear);
        SImpl::MomentumSpacePropagator(w, sqrt(par().m2));
        w *= par().g/N*smear;
        w  = sqrt(vol)*sqrt(w);
    }
    LOG(Message) << "Generating random momentum-space field" << std::endl;
    envGetTmp(Field, phift);
    envGetTmp(ComplexField, ca);
    phift = zero;
    for (int a = 0; a < Nadj; ++a) 
    {
        Site ta;

        gaussian(rng, ca);
        Group::generator(a, ta);
        phift += ca*ta;
    }
    phift *= w;
    LOG(Message) << "Field Fourier transform" << std::endl;
    fft.FFT_all_dim(phi, phift, FFT::backward);
    phi = 0.5*(phi - adj(phi));
    trphi2 = -TensorRemove(sum(trace(phi*phi))).real()/vol;
    LOG(Message) << "tr(phi^2)= " << trphi2 << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_StochFreeField_hpp_
