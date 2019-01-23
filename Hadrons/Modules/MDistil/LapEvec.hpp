/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: Hadrons/Modules/MDistil/LapEvec.hpp

 Copyright (C) 2015-2019
 
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

#ifndef Hadrons_MDistil_LapEvec_hpp_
#define Hadrons_MDistil_LapEvec_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

// These are members of Distillation
#include <Hadrons/Modules/MDistil/Distil.hpp>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MDistil)

/******************************************************************************

 Laplacian eigenvectors - parameters

 ******************************************************************************/

struct StoutParameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(StoutParameters,
                                  int, steps,
                                  double, parm)
  StoutParameters() = default;
  template <class ReaderClass> StoutParameters(Reader<ReaderClass>& Reader){read(Reader,"StoutSmearing",*this);}
};

struct ChebyshevParameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(ChebyshevParameters,
                                  int, PolyOrder,
                                  double, alpha,
                                  double, beta)
  ChebyshevParameters() = default;
  template <class ReaderClass> ChebyshevParameters(Reader<ReaderClass>& Reader){read(Reader,"Chebyshev",*this);}
};

struct LanczosParameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(LanczosParameters,
                                  int, Nstart,
                                  int, Nvec,
                                  int, Nk,
                                  int, Nm,    // Not currently used
                                  int, Np,
                                  int, MaxIt,
                                  int, MinRes,
                                  double, resid)
  LanczosParameters() = default;
  template <class ReaderClass> LanczosParameters(Reader<ReaderClass>& Reader){read(Reader,"Lanczos",*this);}
};

struct DistilParameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(DistilParameters,
                                  int, TI,
                                  int, LI,
                                  int, Nnoise,
                                  int, Ls,       // For makeFiveDimGrid
                                  int, tSrc)
  DistilParameters() = default;
  template <class ReaderClass> DistilParameters(Reader<ReaderClass>& Reader){read(Reader,"Distil",*this);}
};

struct SolverParameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(SolverParameters,
                                  double, CGPrecision,
                                  int, MaxIterations,
                                  double, mass,
                                  double, M5)
  SolverParameters() = default;
  template <class ReaderClass> SolverParameters(Reader<ReaderClass>& Reader){read(Reader,"Solver",*this);}
};

// These are the actual parameters passed to the module during construction

class LapEvecPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(LapEvecPar,
                                  StoutParameters,     Stout,
                                  ChebyshevParameters, Cheby,
                                  LanczosParameters,   Lanczos,
                                  DistilParameters,    Distil,
                                  SolverParameters,    Solver);
};

/******************************************************************************
 
 Laplacian eigenvectors - Module (class) definition
 
 ******************************************************************************/

template <typename FImpl>
class TLapEvec: public Module<LapEvecPar>
{
public:
    // constructor
    TLapEvec();
    TLapEvec(const std::string name);
    // destructor
    virtual ~TLapEvec(void);
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
public:
  GridCartesian * gridLD = nullptr;
};

MODULE_REGISTER_TMP(LapEvec, TLapEvec<FIMPL>, MDistil);

/******************************************************************************
 *                 TLapEvec implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLapEvec<FImpl>::TLapEvec(const std::string name) : Module<LapEvecPar>(name)
{
  // NB: This constructor isn't used!!!
  LOG(Message) << "TLapEvec constructor : TLapEvec<FImpl>::TLapEvec(const std::string name)" << std::endl;
  LOG(Message) << "TLapEvec constructor : Setting gridLD=nullptr" << std::endl;
  gridLD = nullptr;
}

// destructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLapEvec<FImpl>::~TLapEvec()
{
  if( gridLD != nullptr ) {
    LOG(Message) << "Destroying lower dimensional grid" << std::endl;
    delete gridLD;
  }
  else
    LOG(Message) << "Lower dimensional grid was never created" << std::endl;
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLapEvec<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TLapEvec<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLapEvec<FImpl>::setup(void)
{
  LOG(Message) << "setup() : start" << std::endl;
  LOG(Message) << "Stout.steps=" << par().Stout.steps << ", Stout.parm=" << par().Stout.parm << std::endl;
  if( gridLD ) {
    LOG(Message) << "Didn't expect to need to destroy gridLD here!" << std::endl;
    delete gridLD;
  }
  LOG(Message) << "Creating lower dimensional grid" << std::endl;
  gridLD = MakeLowerDimGrid( env().getGrid() );
  LOG(Message) << "setup() : end" << std::endl;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLapEvec<FImpl>::execute(void)
{
  LOG(Message) << "execute() : start" << std::endl;
  LOG(Message) << "Stout.steps=" << par().Stout.steps << ", Stout.parm=" << par().Stout.parm << std::endl;
  LOG(Message) << "execute() : end" << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_LapEvec_hpp_
