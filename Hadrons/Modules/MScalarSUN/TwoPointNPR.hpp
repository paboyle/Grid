/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MScalarSUN/TwoPointNPR.hpp

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
#ifndef Hadrons_MScalarSUN_TwoPointNPR_hpp_
#define Hadrons_MScalarSUN_TwoPointNPR_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MScalarSUN/Utils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         TwoPointNPR                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class TwoPointNPRPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TwoPointNPRPar,
                                    std::vector<std::string>, op,
                                    std::string,              field,
                                    std::string,              output);
};

class TwoPointNPRResult: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TwoPointNPRResult,
                                    std::string,          op,
                                    std::vector<Complex>, data);
};

template <typename SImpl>
class TTwoPointNPR: public Module<TwoPointNPRPar>
{
public:
    typedef typename SImpl::Field                    Field;
    typedef typename SImpl::SiteField::scalar_object Site;
    typedef typename SImpl::ComplexField             ComplexField;
public:
    // constructor
    TTwoPointNPR(const std::string name);
    // destructor
    virtual ~TTwoPointNPR(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(TwoPointNPRSU2, TTwoPointNPR<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_TMP(TwoPointNPRSU3, TTwoPointNPR<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_TMP(TwoPointNPRSU4, TTwoPointNPR<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_TMP(TwoPointNPRSU5, TTwoPointNPR<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_TMP(TwoPointNPRSU6, TTwoPointNPR<ScalarNxNAdjImplR<6>>, MScalarSUN);

/******************************************************************************
 *                 TTwoPointNPR implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TTwoPointNPR<SImpl>::TTwoPointNPR(const std::string name)
: Module<TwoPointNPRPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TTwoPointNPR<SImpl>::getInput(void)
{
    std::vector<std::string> in = par().op;

    in.push_back(par().field);

    return in;
}

template <typename SImpl>
std::vector<std::string> TTwoPointNPR<SImpl>::getOutput(void)
{
    std::vector<std::string> out;
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TTwoPointNPR<SImpl>::setup(void)
{
    const unsigned int nl = env().getDim(0);

    for (unsigned int mu = 1; mu < env().getNd(); ++mu)
    {
        if (nl != env().getDim(mu))
        {
            HADRONS_ERROR(Size, "non-cubic grid");
        }
    }
    envTmpLat(ComplexField, "ftBuf");
    envTmpLat(Field, "ftMatBuf");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TTwoPointNPR<SImpl>::execute(void)
{
    const unsigned int             nd   = env().getNd();
    const unsigned int             nl   = env().getDim(0);
    const Real                     invV = 1./env().getVolume();
    FFT                            fft(envGetGrid(Field));
    std::vector<TwoPointNPRResult> result;
    TwoPointNPRResult              twoPtp1, twoPtp2, twoPtDisc;
    auto                           &phi    = envGet(Field, par().field);
    bool                           doAux = true;

    envGetTmp(ComplexField, ftBuf);
    envGetTmp(Field, ftMatBuf);
    LOG(Message) << "FFT: field '" << par().field << "'" << std::endl;
    fft.FFT_all_dim(ftMatBuf, phi, FFT::forward);
    for (auto &opName: par().op)
    {
        auto              &op = envGet(ComplexField, opName);
        std::vector<int>  p1, p2, p;
        Site              phip1, phip2;
        TComplex          opp;
        TwoPointNPRResult r, rDisc;

        LOG(Message) << "FFT: operator '" << opName << "'" << std::endl;
        fft.FFT_all_dim(ftBuf, op, FFT::forward);
        LOG(Message) << "Generating vertex function" << std::endl;
        r.op = opName;
        r.data.resize(nl);
        rDisc.op = opName + "_disc";
        rDisc.data.resize(nl);
        if (doAux)
        {
            twoPtp1.op = "phi_prop_p1";
            twoPtp1.data.resize(nl);
            twoPtp2.op = "phi_prop_p2";
            twoPtp2.data.resize(nl);
            twoPtDisc.op = "phi_prop_disc";
            twoPtDisc.data.resize(nl);
        }
        for (unsigned int n = 0; n < nl; ++n)
        {
            p1.assign(nd, 0);
            p2.assign(nd, 0);
            p.assign(nd, 0);
            // non-exceptional RI/SMOM kinematic
            // p1 = mu*(1,1,0): in mom
            // p2 = mu*(0,1,1): out mom
            // p  = p1 - p2 = mu*(1,0,-1)
            // mu = 2*n*pi/L
            p1[0] = n;
            p1[1] = n;
            p2[1] = n;
            p2[2] = n;
            p[0]  = n;
            p[2]  = (nl - n) % nl;
            peekSite(phip1, ftMatBuf, p1);
            peekSite(phip2, ftMatBuf, p2);
            peekSite(opp, ftBuf, p);
            if (doAux)
            {
                twoPtp1.data[n]   = invV*TensorRemove(trace(phip1*adj(phip1)));
                twoPtp2.data[n]   = invV*TensorRemove(trace(phip2*adj(phip2)));
                twoPtDisc.data[n] = invV*TensorRemove(trace(phip2*adj(phip1)));
            }
            r.data[n]     = invV*TensorRemove(trace(phip2*adj(phip1))*opp);
            rDisc.data[n] = invV*TensorRemove(trace(phip1*adj(phip1))*opp);
        }
        if (doAux)
        {
            result.push_back(twoPtp1);
            result.push_back(twoPtp2);
            result.push_back(twoPtDisc);
        }
        result.push_back(r);
        result.push_back(rDisc);
        doAux = false;
    }
    saveResult(par().output, "twoptnpr", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_TwoPointNPR_hpp_
