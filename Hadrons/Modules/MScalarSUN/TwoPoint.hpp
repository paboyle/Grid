/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MScalarSUN/TwoPoint.hpp

Copyright (C) 2015-2018

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
#ifndef Hadrons_MScalarSUN_TwoPoint_hpp_
#define Hadrons_MScalarSUN_TwoPoint_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MScalarSUN/Utils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                 2-pt functions for a given set of operators                *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class TwoPointPar: Serializable
{
public:
    typedef std::pair<std::string, std::string> OpPair;
    GRID_SERIALIZABLE_CLASS_MEMBERS(TwoPointPar,
                                    std::vector<OpPair>,      op,
                                    std::vector<std::string>, mom,
                                    std::string,              output);
};

class TwoPointResult: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TwoPointResult,
                                    std::string, sink,
                                    std::string, source,
                                    std::vector<int>, mom,
                                    std::vector<Complex>, data);
};

template <typename SImpl>
class TTwoPoint: public Module<TwoPointPar>
{
public:
    typedef typename SImpl::Field         Field;
    typedef typename SImpl::ComplexField  ComplexField;
    typedef          std::vector<Complex> SlicedOp;
public:
    // constructor
    TTwoPoint(const std::string name);
    // destructor
    virtual ~TTwoPoint(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::vector<std::vector<int>> mom_;
};

MODULE_REGISTER_TMP(TwoPointSU2, TTwoPoint<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_TMP(TwoPointSU3, TTwoPoint<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_TMP(TwoPointSU4, TTwoPoint<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_TMP(TwoPointSU5, TTwoPoint<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_TMP(TwoPointSU6, TTwoPoint<ScalarNxNAdjImplR<6>>, MScalarSUN);

/******************************************************************************
 *                       TTwoPoint implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TTwoPoint<SImpl>::TTwoPoint(const std::string name)
: Module<TwoPointPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TTwoPoint<SImpl>::getInput(void)
{   
    std::vector<std::string> in;
    std::set<std::string>    ops;

    for (auto &p: par().op)
    {
        ops.insert(p.first);
        ops.insert(p.second);
    }
    for (auto &o: ops)
    {
        in.push_back(o);
    }

    return in;
}

template <typename SImpl>
std::vector<std::string> TTwoPoint<SImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TTwoPoint<SImpl>::setup(void)
{
    const unsigned int nd = env().getDim().size();

    mom_.resize(par().mom.size());
    for (unsigned int i = 0; i < mom_.size(); ++i)
    {
        mom_[i] = strToVec<int>(par().mom[i]);
        if (mom_[i].size() != nd - 1)
        {
            HADRONS_ERROR(Size, "momentum number of components different from " 
                               + std::to_string(nd-1));
        }
        for (unsigned int j = 0; j < nd - 1; ++j)
        {
            mom_[i][j] = (mom_[i][j] + env().getDim(j)) % env().getDim(j);
        }
    }
    envTmpLat(ComplexField, "ftBuf");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TTwoPoint<SImpl>::execute(void)
{
    LOG(Message) << "Computing 2-point functions" << std::endl;
    for (auto &p: par().op)
    {
        LOG(Message) << "  <" << p.first << " " << p.second << ">" << std::endl;
    }

    const unsigned int                           nd      = env().getNd();
    const unsigned int                           nt      = env().getDim().back();
    const unsigned int                           nop     = par().op.size();
    const unsigned int                           nmom    = mom_.size();
    double                                       partVol = 1.;
    std::vector<int>                             dMask(nd, 1);
    std::set<std::string>                        ops;
    std::vector<TwoPointResult>                  result;
    std::map<std::string, std::vector<SlicedOp>> slicedOp;
    FFT                                          fft(envGetGrid(Field));
    TComplex                                     buf;

    envGetTmp(ComplexField, ftBuf);
    dMask[nd - 1] = 0;
    for (unsigned int mu = 0; mu < nd - 1; ++mu)
    {
        partVol *= env().getDim()[mu];
    }
    for (auto &p: par().op)
    {
        ops.insert(p.first);
        ops.insert(p.second);
    }
    for (auto &o: ops)
    {
        auto &op = envGet(ComplexField, o);

        slicedOp[o].resize(nmom);
        LOG(Message) << "Operator '" << o << "' FFT" << std::endl;
        fft.FFT_dim_mask(ftBuf, op, dMask, FFT::forward);
        for (unsigned int m = 0; m < nmom; ++m)
        {
            auto qt = mom_[m];

            qt.resize(nd);
            slicedOp[o][m].resize(nt);
            for (unsigned int t = 0; t < nt; ++t)
            {
                qt[nd - 1] = t;
                peekSite(buf, ftBuf, qt);
                slicedOp[o][m][t] = TensorRemove(buf);
            }
        }
    }
    LOG(Message) << "Making contractions" << std::endl;
    for (unsigned int m = 0; m < nmom; ++m)
    for (auto &p: par().op)
    {
        TwoPointResult r;

        r.sink   = p.first;
        r.source = p.second;
        r.mom    = mom_[m];
        r.data   = makeTwoPoint(slicedOp[p.first][m], slicedOp[p.second][m], 
                                1./partVol);
        result.push_back(r);
    }
    saveResult(par().output, "twopt", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_TwoPoint_hpp_
