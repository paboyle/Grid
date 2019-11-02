/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: Hadrons/Modules/MDistil/Noises.hpp
 
 Copyright (C) 2019
 
 Author: Felix Erben <ferben@ed.ac.uk>
 Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 
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

#ifndef Hadrons_MDistil_Noises_hpp_
#define Hadrons_MDistil_Noises_hpp_

#include <Hadrons/Modules/MDistil/DistilCommon.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MDistil)

/******************************************************************************
 *                         Noises                                 *
 ******************************************************************************/

class NoisesPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(NoisesPar,
                                    int, nnoise,
                                    int, nvec,
                                    std::string, UniqueIdentifier,
                                    std::string, TI,
                                    std::string, LI)
};

template <typename FImpl>
class TNoises: public Module<NoisesPar>
{
public:
    // constructor
    TNoises(const std::string name);
    // destructor
    virtual ~TNoises(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Noises, TNoises<FIMPL>, MDistil);

/******************************************************************************
 *                 TNoises implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TNoises<FImpl>::TNoises(const std::string name)
: Module<NoisesPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TNoises<FImpl>::getInput(void)
{
    return {};
}

template <typename FImpl>
std::vector<std::string> TNoises<FImpl>::getOutput(void)
{
    return {getName()};
}

// setup ///////////////////////////////////////////////////////////////////////

template <typename FImpl>
void TNoises<FImpl>::setup(void)
{
    const int Nt{env().getDim(Tdir)};
    const int nnoise{par().nnoise};
    const int nvec{par().nvec};
    envCreate(NoiseTensor, getName(), 1, nnoise, Nt, nvec, Ns);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TNoises<FImpl>::execute(void)
{
    const int Nt{env().getDim(Tdir)};
    const int nnoise{par().nnoise};
    const int nvec{par().nvec};
    const int TI{ Hadrons::MDistil::DistilParameters::ParameterDefault( par().TI, Nt, false) };
    const int LI{ Hadrons::MDistil::DistilParameters::ParameterDefault( par().LI, nvec, false) };
    const bool full_tdil{ TI == Nt }; \
    const bool exact_distillation{ full_tdil && LI == nvec }; \
    std::string UniqueIdentifier{par().UniqueIdentifier};
    if (UniqueIdentifier.empty())
        UniqueIdentifier = getName();
    UniqueIdentifier.append( std::to_string( vm().getTrajectory() ) );
    
    // We use our own seeds so we can specify different noises per quark
    GridSerialRNG sRNG;
    sRNG.SeedUniqueString(UniqueIdentifier);
    Real rn;
    auto &noise = envGet(NoiseTensor, getName());
    for (int inoise = 0; inoise < nnoise; inoise++) {
        for (int t = 0; t < Nt; t++) {
            for (int ivec = 0; ivec < nvec; ivec++) {
                for (int is = 0; is < Ns; is++) {
                    if (exact_distillation)
                        noise(inoise, t, ivec, is) = 1.;
                    else{
                        random(sRNG,rn);
                        // We could use a greater number of complex roots of unity
                        // ... but this seems to work well
                        noise(inoise, t, ivec, is) = (rn > 0.5) ? -1 : 1;
                    }
                }
            }
        }
    }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif // Hadrons_MDistil_Noises_hpp_
