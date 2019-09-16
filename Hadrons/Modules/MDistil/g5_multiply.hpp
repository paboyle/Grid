/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: Hadrons/Modules/MDistil/g5_multiply.hpp
 
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

#ifndef Hadrons_MDistil_g5_multiply_hpp_
#define Hadrons_MDistil_g5_multiply_hpp_

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
 *                         g5_multiply                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MDistil)

class g5_multiplyPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(g5_multiplyPar,
                                    std::string, input,
                                    int, nnoise,
                                    int, LI,
                                    int, Ns,
                                    int, Nt_inv);
};

template <typename FImpl>
class Tg5_multiply: public Module<g5_multiplyPar>
{
	public:
		    FERM_TYPE_ALIASES(FImpl,);

public:
    // constructor
    Tg5_multiply(const std::string name);
    // destructor
    virtual ~Tg5_multiply(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(g5_multiply, Tg5_multiply<FIMPL>, MDistil);

/******************************************************************************
 *                 Tg5_multiply implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
Tg5_multiply<FImpl>::Tg5_multiply(const std::string name)
: Module<g5_multiplyPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> Tg5_multiply<FImpl>::getInput(void)
{
    std::vector<std::string> in;

    in.push_back(par().input);

    return in;
}

template <typename FImpl>
std::vector<std::string> Tg5_multiply<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void Tg5_multiply<FImpl>::setup(void)
{
   int nnoise=par().nnoise;
   int LI=par().LI;
   int Ns=par().Ns;
   int Nt_inv=par().Nt_inv;

   envCreate(std::vector<FermionField>, getName(), 1, 
		                                    nnoise*LI*Ns*Nt_inv, envGetGrid(FermionField));

}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void Tg5_multiply<FImpl>::execute(void)
{
     auto        &input_vector       = envGet(std::vector<FermionField>, par().input);
     auto        &output_vector       = envGet(std::vector<FermionField>, getName());

     Gamma g5(Gamma::Algebra::Gamma5);

   int nnoise=par().nnoise;
   int LI=par().LI;
   int Ns=par().Ns;
   int Nt_inv=par().Nt_inv;
   int Ntot = nnoise*LI*Ns*Nt_inv;
   for (int i =0; i<Ntot;i++){
     output_vector[i] = g5*input_vector[i];
   } 

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_g5_multiply_hpp_
