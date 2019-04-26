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

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/DilutedNoise.hpp>

// These are members of Distillation
 #include <Hadrons/Distil.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Noises                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MDistil)

class NoisesPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(NoisesPar,
                                    std::string, UniqueIdentifier,
                                    int, nvec,
			            DistilParameters, Distil);
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
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TNoises<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TNoises<FImpl>::setup(void)
{
    const DistilParameters & Distil{par().Distil};
    const int nvec{par().nvec};
 
    envCreate(std::vector<Complex>, getName(), 1, nvec*Distil.Ns*Distil.Nt*Distil.nnoise);

}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TNoises<FImpl>::execute(void)
{
       const std::string &UniqueIdentifier{par().UniqueIdentifier};
       auto &noise   = envGet(std::vector<Complex>, getName());
       const int nvec{par().nvec};
       const DistilParameters & Distil{par().Distil};
       const int nnoise{Distil.nnoise};
       const int Nt{Distil.Nt};
       const int Ns{Distil.Ns};
       const int TI{Distil.TI};
       const int LI{Distil.LI};
       const bool full_tdil{TI==Nt};
       const bool exact_distillation{full_tdil && LI==nvec};

       GridSerialRNG sRNG; 
       sRNG.SeedUniqueString(UniqueIdentifier + std::to_string(vm().getTrajectory())); //maybe add more??
       Real rn;
		       
       for (int inoise=0;inoise<nnoise;inoise++) {
         for (int t=0;t<Nt;t++) {
           for (int ivec=0;ivec<nvec;ivec++) {
             for (int is=0;is<Ns;is++) {
               if (exact_distillation)
                 noise[inoise + nnoise*(t + Nt*(ivec+nvec*is))] = 1.;
               else{
                 random(sRNG,rn);
                 // We could use a greater number of complex roots of unity
                 // ... but this seems to work well
                 noise[inoise + nnoise*(t + Nt*(ivec+nvec*is))] = (rn > 0.5) ? -1 : 1;
               }
	     }
	   }
	 }
       }

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_Noises_hpp_
