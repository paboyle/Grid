/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MContraction/Meson.hpp

Copyright (C) 2015
Copyright (C) 2016

Author: Antonin Portelli <antonin.portelli@me.com>
        Andrew Lawson    <andrew.lawson1991@gmail.com>

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

#ifndef Hadrons_Meson_hpp_
#define Hadrons_Meson_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Meson contractions
 -----------------------------
 
 * options:
 - q1: input propagator 1 (string)
 - q2: input propagator 2 (string)
 - gammas: gamma products to insert at source & sink, pairs of gamma matrices 
           (space-separated integers) in square brackets, in a sequence 
           (e.g. "[15 7][7 15][7 7]").

           Special values: "all" - perform all possible contractions.
 
 */

/******************************************************************************
 *                                TMeson                                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class MesonPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MesonPar,
                                    std::string, q1,
                                    std::string, q2,
                                    std::string, gammas,
                                    std::string, output);
};

template <typename FImpl1, typename FImpl2>
class TMeson: public Module<MesonPar>
{
public:
    TYPE_ALIASES(FImpl1, 1);
    TYPE_ALIASES(FImpl2, 2);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::vector<std::vector<std::vector<Complex>>>, corr);
    };
public:
    // constructor
    TMeson(const std::string name);
    // destructor
    virtual ~TMeson(void) = default;
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void parseGammaString(SpinMatrix*,
                                  unsigned int &, 
                                  unsigned int &,
                                  std::vector<std::vector<bool>> &);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(Meson, ARG(TMeson<FIMPL, FIMPL>), MContraction);

/******************************************************************************
 *                           TMeson implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
TMeson<FImpl1, FImpl2>::TMeson(const std::string name)
: Module<MesonPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
std::vector<std::string> TMeson<FImpl1, FImpl2>::getInput(void)
{
    std::vector<std::string> input = {par().q1, par().q2};
    
    return input;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TMeson<FImpl1, FImpl2>::getOutput(void)
{
    std::vector<std::string> output = {getName()};
    
    return output;
}

template <typename T>
std::vector<std::pair<T, T>> strToVecPair(const std::string s)
{
    std::vector<std::pair<T, T>> v;
    return v;
}

template <typename FImpl1, typename FImpl2>
void TMeson<FImpl1, FImpl2>::parseGammaString(SpinMatrix *g,
                                              unsigned int &n_snk, 
                                              unsigned int &n_src,
                                              std::vector<std::vector<bool>> &toDo)
{
    // Initialise counters for parsing gamma insertions at source & sink.
    int empty = -1;
    std::vector<int> gamma_inds(Ns*Ns, empty);
    unsigned int n_gam = 0;

    // Determine gamma matrices to insert at source/sink.
    if (par().gammas.compare("all") == 0)
    {
        // Do all contractions.
        toDo.resize(Ns*Ns);
        for (int i = 0; i < Ns*Ns; ++i)
        {
            g[i] = makeGammaProd(i);
            toDo[i].assign(Ns*Ns, true);
        }
    }
    else
    {
        // Parse individual contractions from input string.
        std::vector<std::pair<int, int>> gamma_pairs;
        gamma_pairs = strToVecPair<int>(par().gammas);

        // Function for gamma matrix counting & indexing at source/sink.
        auto index_gamma = [&empty, &g, &gamma_inds, &n_gam](int i, 
                                                             unsigned int &n)
        {
            if (i >= gamma_inds.size())
            {
                HADRON_ERROR("Invalid gamma matrix index " << i);
            }
            if (gamma_inds[i] == empty)
            {
                g[n_gam] = makeGammaProd(i);
                gamma_inds[i] = n_gam;
                ++n_gam;
                ++n;
            }
        };

        // Count no. of unique gamma matrices, then construct matrix of
        // contractions to do.
        for (unsigned int i = 0; i < gamma_inds.size(); ++i)
        {
            index_gamma(gamma_pairs[i].first, n_snk);
            index_gamma(gamma_pairs[i].second, n_src);
        }
        toDo.resize(n_gam);
        for (int i = 0; i < n_gam; ++i)
        {
            toDo[i].assign(n_gam, false);
        }
        for (int i = 0; i < gamma_inds.size(); ++i)
        {
            toDo[gamma_inds[gamma_pairs[i].first]]
                [gamma_inds[gamma_pairs[i].second]] = true;
        }
    }
}


// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TMeson<FImpl1, FImpl2>::execute(void)
{
    LOG(Message) << "Computing meson contractions '" << getName() << "' using"
                 << " quarks '" << par().q1 << "' and '" << par().q2 << "'"
                 << std::endl;
    
    XmlWriter             writer(par().output);
    PropagatorField1      &q1 = *env().template getObject<PropagatorField1>(par().q1);
    PropagatorField2      &q2 = *env().template getObject<PropagatorField2>(par().q2);
    LatticeComplex        c(env().getGrid());
    SpinMatrix            g[Ns*Ns], g5;
    std::vector<std::vector<bool>> toDo;
    std::vector<TComplex> buf;
    Result                result;
    unsigned int          n_snk, n_src;

    g5 = makeGammaProd(Ns*Ns - 1);
    parseGammaString(g, n_snk, n_src, toDo);

    result.corr.resize(n_snk);
    for (unsigned int iSink = 0; iSink < toDo.size(); ++iSink)
    {
        result.corr[iSink].resize(n_src);
        for (unsigned int iSrc = 0; iSrc < toDo.size(); ++iSrc)
        {
            if (toDo[iSink][iSrc])
            {
                c = trace(g[iSink]*q1*g[iSrc]*g5*adj(q2)*g5);
                sliceSum(c, buf, Tp);
                result.corr[iSink][iSrc].resize(buf.size());
                for (unsigned int t = 0; t < buf.size(); ++t)
                {
                    result.corr[iSink][iSrc][t] = TensorRemove(buf[t]);
                }
            }
        }
    }
    write(writer, "meson", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_Meson_hpp_
