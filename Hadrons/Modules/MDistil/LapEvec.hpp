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
#include <Hadrons/EigenPack.hpp>

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
                                  std::string,         gauge,
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
  GAUGE_TYPE_ALIASES(FImpl,);
  typedef std::vector<Grid::Hadrons::EigenPack<LatticeColourVector> > DistilEP;

public:
  // constructor
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
protected:
  // These variables are created in setup() and freed in Cleanup()
  GridCartesian * gridLD; // Owned by me, so I must delete it
  GridCartesian * gridHD; // Owned by environment (so I won't delete it)
  int Nx, Ny, Nz, Nt;

protected:
  void Cleanup(void);
};

MODULE_REGISTER_TMP(LapEvec, TLapEvec<FIMPL>, MDistil);

/******************************************************************************
 TLapEvec implementation
 ******************************************************************************/

// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLapEvec<FImpl>::TLapEvec(const std::string name) : gridLD{nullptr}, Module<LapEvecPar>(name)
{
  LOG(Message) << "TLapEvec constructor : TLapEvec<FImpl>::TLapEvec(const std::string name)" << std::endl;
  LOG(Message) << "this=" << this << ", gridLD=" << gridLD << std::endl;
}

// destructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLapEvec<FImpl>::~TLapEvec()
{
  Cleanup();
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLapEvec<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
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
  Cleanup();
  Environment & e{env()};
  gridHD = e.getGrid();
  gridLD = MakeLowerDimGrid( gridHD );
  Nx = gridHD->_gdimensions[Xdir];
  Ny = gridHD->_gdimensions[Ydir];
  Nz = gridHD->_gdimensions[Zdir];
  Nt = gridHD->_gdimensions[Tdir];
  // Temporaries
  envTmpLat(GaugeField, "Umu");
  envTmpLat(GaugeField, "Umu_stout");
  envTmpLat(GaugeField, "Umu_smear");
  // Output objects
  envCreate(DistilEP, getName(), 1, Nt);
}

// clean up any temporaries created by setup (that aren't stored in the environment)
template <typename FImpl>
void TLapEvec<FImpl>::Cleanup(void)
{
  if( gridLD != nullptr ) {
    delete gridLD;
    gridLD = nullptr;
  }
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
