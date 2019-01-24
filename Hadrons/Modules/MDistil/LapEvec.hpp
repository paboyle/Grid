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
                                  std::string,         EigenPackName,
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

//constexpr char szEigenPackSuffix[] = "_eigenPack";

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
    std::vector<std::string> out = {getName()}; // This is the higher dimensional eigenpack
    
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
  Nx = gridHD->_fdimensions[Xdir];
  Ny = gridHD->_fdimensions[Ydir];
  Nz = gridHD->_fdimensions[Zdir];
  Nt = gridHD->_fdimensions[Tdir];
  // Temporaries
  envTmpLat(GaugeField, "Umu");
  envTmpLat(GaugeField, "Umu_stout");
  envTmpLat(GaugeField, "Umu_smear");
  envTmp(LatticeGaugeField, "UmuNoTime",1,LatticeGaugeField(gridLD));
  envTmp(LatticeColourVector, "src",1,LatticeColourVector(gridLD));
  envTmp(std::vector<DistilEP>, "eig",1,std::vector<DistilEP>(Nt));
  // Output objects
  envCreate(DistilEP, getName(), 1, par().Lanczos.Nvec, gridHD );
}

// clean up any temporaries created by setup (that aren't stored in the environment)
template <typename FImpl>
void TLapEvec<FImpl>::Cleanup(void)
{
  if( gridLD != nullptr ) {
    delete gridLD;
    gridLD = nullptr;
  }
  gridHD = nullptr;
}

/******************************************************************************
 Calculate low-mode eigenvalues of the Laplacian
 ******************************************************************************/

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLapEvec<FImpl>::execute(void)
{
  LOG(Message) << "execute() : start" << std::endl;

  // Alii for parameters
  const int &TI{par().Distil.TI};
  const int &LI{par().Distil.LI};
  const int &nnoise{par().Distil.Nnoise};
  const int &tsrc{par().Distil.tSrc};
  const LanczosParameters &LPar{par().Lanczos};
  const int &nvec{LPar.Nvec};
  const bool exact_distillation{TI==Nt && LI==nvec};
  const bool full_tdil{TI==Nt};
  const int &Nt_inv{full_tdil ? 1 : TI};
  const ChebyshevParameters &ChebPar{par().Cheby};

  // Assertions on the parameters we read
  assert(TI>1);
  assert(LI>1);
  if(exact_distillation)
    assert(nnoise==1);
  else
    assert(nnoise>1);

  // Stout smearing
  envGetTmp(GaugeField, Umu);
  envGetTmp(GaugeField, Umu_smear);
  LOG(Message) << "Initial plaquette: " << WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu) << std::endl;
  {
    envGetTmp(GaugeField, Umu_stout);
    const int &Steps{par().Stout.steps};
    Smear_Stout<PeriodicGimplR> LS(par().Stout.parm);
    for (int i = 0; i < Steps; i++) {
      LS.smear(Umu_stout, Umu_smear);
      Umu_smear = Umu_stout;
    }
  }
  LOG(Message) << "Smeared plaquette: " << WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu_smear) << std::endl;

  // For debugging only, write logging output to a local file
  std::ofstream * ll = nullptr;
  const int rank{gridHD->ThisRank()};
  if((0)) { // debug to a local log file
    std::string filename{"Local_"};
    filename.append(std::to_string(rank));
    filename.append(".log");
    ll = new std::ofstream(filename);
  }

  ////////////////////////////////////////////////////////////////////////
  // Invert Peardon Nabla operator separately on each time-slice
  ////////////////////////////////////////////////////////////////////////
  
  std::string sEigenPackName(par().EigenPackName);
  bool bReturnValue = true;
  auto & eig4d = envGet(DistilEP, getName() );
  eig4d.resize(nvec,gridHD);
  envGetTmp(std::vector<DistilEP>, eig);   // Eigenpack for each timeslice
  envGetTmp(LatticeGaugeField, UmuNoTime); // Gauge field without time dimension
  envGetTmp(LatticeColourVector, src);
  const int Ntlocal{gridHD->LocalDimensions()[Tdir]};
  const int Ntfirst{gridHD->LocalStarts()[Tdir]};
  for(int t=Ntfirst;bReturnValue && t<Ntfirst+Ntlocal;t++){
    std::cout << GridLogMessage << "------------------------------------------------------------" << std::endl;
    std::cout << GridLogMessage << " Compute eigenpack, Timeslice  = " << t << std::endl;
    std::cout << GridLogMessage << "------------------------------------------------------------" << std::endl;
    
    LOG(Message) << "eig.size()=" << eig.size() << std::endl;
    eig[t].resize(LPar.Nk+LPar.Np,gridLD);
    LOG(Message) << "After eig[t].resize" << std::endl;
    
    // Construct smearing operator
    ExtractSliceLocal(UmuNoTime,Umu_smear,0,t-Ntfirst,Grid::QCD::Tdir); // switch to 3d/4d objects
    LinOpPeardonNabla<LatticeColourVector> PeardonNabla(UmuNoTime);
    std::cout << "Chebyshev preconditioning to order " << ChebPar.PolyOrder
    << " with parameters (alpha,beta) = (" << ChebPar.alpha << "," << ChebPar.beta << ")" << std::endl;
    Chebyshev<LatticeColourVector> Cheb(ChebPar.alpha,ChebPar.beta,ChebPar.PolyOrder);
    
    //from Test_Cheby.cc
    if ( ((0)) && Ntfirst == 0 && t==0) {
      std::ofstream of("cheby_" + std::to_string(ChebPar.alpha) + "_" + std::to_string(ChebPar.beta) + "_" + std::to_string(ChebPar.PolyOrder));
      Cheb.csv(of);
    }

    // Construct source vector according to Test_dwf_compressed_lanczos.cc
    src=11.0;
    RealD nn = norm2(src);
    nn = Grid::sqrt(nn);
    src = src * (1.0/nn);

    GridLogIRL.Active(1);
    LinOpPeardonNablaHerm<LatticeColourVector> PeardonNablaCheby(Cheb,PeardonNabla);
    ImplicitlyRestartedLanczos<LatticeColourVector> IRL(PeardonNablaCheby,PeardonNabla,LPar.Nvec,LPar.Nk,LPar.Nk+LPar.Np,LPar.resid,LPar.MaxIt);
    int Nconv = 0;
    
    if(ll) *ll << t << " : Before IRL.calc()" << std::endl;
    IRL.calc(eig[t].eval,eig[t].evec,src,Nconv);
    if(ll) *ll << t << " : After  IRL.calc()" << std::endl;
    if( Nconv < LPar.Nvec ) {
      bReturnValue = false;
      if(ll) *ll << t << " : Convergence error : Only " << Nconv << " converged!" << std::endl;
    } else {
      if( Nconv > LPar.Nvec )
        eig[t].resize( LPar.Nvec, gridLD );
      std::cout << GridLogMessage << "Timeslice " << t << " has " << eig[t].eval.size() << " eigenvalues and " << eig[t].evec.size() << " eigenvectors." << std::endl;
      
      // Now rotate the eigenvectors into our phase convention
      RotateEigen( eig[t].evec );
      
      // Write the eigenvectors and eigenvalues to disk
      //std::cout << GridLogMessage << "Writing eigenvalues/vectors to " << pszEigenPack << std::endl;
      eig[t].record.operatorXml="<OPERATOR>Michael</OPERATOR>";
      eig[t].record.solverXml="<SOLVER>Felix</SOLVER>";
      eig[t].write(sEigenPackName,false,t);
      //std::cout << GridLogMessage << "Written eigenvectors" << std::endl;
    }
    for (int i=0;i<LPar.Nvec;i++){
      std::cout << "Inserting Timeslice " << t << " into vector " << i << std::endl;
      InsertSliceLocal(eig[t].evec[i],eig4d.evec[i],0,t,3);
    }
  }

  // Close the local debugging log file
  if( ll ) {
    *ll << " Returning " << bReturnValue << std::endl;
    delete ll;
  }
  LOG(Message) << "execute() : end" << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_LapEvec_hpp_
