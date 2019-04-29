/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: Hadrons/Modules/MDistil/LapEvec.hpp

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

#ifndef Hadrons_MDistil_LapEvec_hpp_
#define Hadrons_MDistil_LapEvec_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

// These are members of Distillation
#include <Hadrons/Distil.hpp>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MDistil)

/******************************************************************************

 Laplacian eigenvectors - parameters

 ******************************************************************************/

struct StoutParameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(StoutParameters,
                                  int, steps,
                                  double, parm) // TODO: change name of this to rho
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
                                  //int, Nstart,
                                  int, Nvec,
                                  int, Nk,
                                  //int, Nm,    // Not currently used
                                  int, Np,
                                  int, MaxIt,
                                  //int, MinRes,
                                  double, resid,
                                  std::string, Log) // Any non-empty string will enable logging
  LanczosParameters() = default;
  template <class ReaderClass> LanczosParameters(Reader<ReaderClass>& Reader){read(Reader,"Lanczos",*this);}
};

// These are the actual parameters passed to the module during construction

class LapEvecPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(LapEvecPar,
                                  std::string,         gauge,
 //                                 std::string, ConfigFileDir,
  //                                std::string, ConfigFileName,
                                  //,std::string,         EigenPackName
                                  StoutParameters,     Stout
                                  ,ChebyshevParameters, Cheby
                                  ,LanczosParameters,   Lanczos
                                  //,DistilParameters,    Distil
                                  )//,SolverParameters,    Solver)
};

/******************************************************************************
 
 Laplacian eigenvectors - Module (class) definition
 
 ******************************************************************************/

template <typename GImpl>
class TLapEvec: public Module<LapEvecPar>
{
public:
  GAUGE_TYPE_ALIASES(GImpl,);
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
  virtual void Cleanup(void);
};

MODULE_REGISTER_TMP(LapEvec, TLapEvec<GIMPL>, MDistil);

/******************************************************************************
 TLapEvec implementation
 ******************************************************************************/

// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TLapEvec<GImpl>::TLapEvec(const std::string name) : gridLD{nullptr}, Module<LapEvecPar>(name)
{
}

// destructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TLapEvec<GImpl>::~TLapEvec()
{
  Cleanup();
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TLapEvec<GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    return in;
}

template <typename GImpl>
std::vector<std::string> TLapEvec<GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()}; // This is the higher dimensional eigenpack
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TLapEvec<GImpl>::setup(void)
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
  //envTmpLat(GaugeField, "Umu");
  envTmpLat(GaugeField, "Umu_stout");
  envTmpLat(GaugeField, "Umu_smear");
  envTmp(LatticeGaugeField, "UmuNoTime",1,LatticeGaugeField(gridLD));
  envTmp(LatticeColourVector, "src",1,LatticeColourVector(gridLD));
  envTmp(std::vector<DistilEP>, "eig",1,std::vector<DistilEP>(Nt));
  // Output objects
  envCreate(DistilEP, getName(), 1, par().Lanczos.Nvec, gridHD );
}

// clean up any temporaries created by setup (that aren't stored in the environment)
template <typename GImpl>
void TLapEvec<GImpl>::Cleanup(void)
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
template <typename GImpl>
void TLapEvec<GImpl>::execute(void)
{
  const ChebyshevParameters &ChebPar{par().Cheby};
  const LanczosParameters   &LPar{par().Lanczos};
  const int &nvec{LPar.Nvec};

  // Enable IRL logging if requested
  if( LPar.Log.size() > 0 )
    GridLogIRL.Active(1);

  //const bool exact_distillation{TI==Nt && LI==nvec};
  //const bool full_tdil{TI==Nt};
  //const int &Nt_inv{full_tdil ? 1 : TI};

  // Assertions on the parameters we read
  //assert(TI>1);
  //assert(LI>1);
  //if(exact_distillation)
    //assert(nnoise==1);
  //else
    //assert(nnoise>1);

  auto &Umu = envGet(GaugeField, par().gauge);
  envGetTmp(GaugeField, Umu_smear);

  // Stout smearing
  Umu_smear = Umu;
  LOG(Message) << "Initial plaquette: " << WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu) << std::endl;
  {
    const StoutParameters &Stout{par().Stout};
    envGetTmp(GaugeField, Umu_stout);
    Smear_Stout<PeriodicGimplR> LS(Stout.parm, Tdir); // spatial smearing only
    for (int i = 0; i < Stout.steps; i++) {
      LS.smear(Umu_stout, Umu_smear);
      Umu_smear = Umu_stout;
    }
  }
  LOG(Message) << "Smeared plaquette: " << WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu_smear) << std::endl;

  ////////////////////////////////////////////////////////////////////////
  // Invert Peardon Nabla operator separately on each time-slice
  ////////////////////////////////////////////////////////////////////////
  
  auto & eig4d = envGet(DistilEP, getName() );
  envGetTmp(std::vector<DistilEP>, eig);   // Eigenpack for each timeslice
  envGetTmp(LatticeGaugeField, UmuNoTime); // Gauge field without time dimension
  envGetTmp(LatticeColourVector, src);
  const int Ntlocal{gridHD->LocalDimensions()[Tdir]};
  const int Ntfirst{gridHD->LocalStarts()[Tdir]};
  const char DefaultOperatorXml[] = "<OPERATOR>Michael</OPERATOR>";
  const char DefaultsolverXml[]   = "<SOLVER>Felix</SOLVER>";
  for(int t = Ntfirst; t < Ntfirst + Ntlocal; t++ ) {
    LOG(Message) << "------------------------------------------------------------" << std::endl;
    LOG(Message) << " Compute eigenpack, Timeslice = " << t << " / " << Ntfirst + Ntlocal << std::endl;
    LOG(Message) << "------------------------------------------------------------" << std::endl;
    eig[t].resize(LPar.Nk+LPar.Np,gridLD);
    
    // Construct smearing operator
    ExtractSliceLocal(UmuNoTime,Umu_smear,0,t-Ntfirst,Grid::QCD::Tdir); // switch to 3d/4d objects
    LinOpPeardonNabla<LatticeColourVector> PeardonNabla(UmuNoTime);
    LOG(Debug) << "Chebyshev preconditioning to order " << ChebPar.PolyOrder
      << " with parameters (alpha,beta) = (" << ChebPar.alpha << "," << ChebPar.beta << ")" << std::endl;
    Chebyshev<LatticeColourVector> Cheb(ChebPar.alpha,ChebPar.beta,ChebPar.PolyOrder);
    
    //from Test_Cheby.cc
    //if( Ntfirst == 0 && t==0) {
      //std::ofstream of("cheby_" + std::to_string(ChebPar.alpha) + "_" + std::to_string(ChebPar.beta) + "_" + std::to_string(ChebPar.PolyOrder));
      //Cheb.csv(of);
    //}

    // Construct source vector according to Test_dwf_compressed_lanczos.cc
    src = 11.0;
    RealD nn = norm2(src);
    nn = Grid::sqrt(nn);
    src = src * (1.0/nn);

    LinOpPeardonNablaHerm<LatticeColourVector> PeardonNablaCheby(Cheb,PeardonNabla);
    ImplicitlyRestartedLanczos<LatticeColourVector>
      IRL(PeardonNablaCheby,PeardonNabla,LPar.Nvec,LPar.Nk,LPar.Nk+LPar.Np,LPar.resid,LPar.MaxIt);
    int Nconv = 0;
    IRL.calc(eig[t].eval,eig[t].evec,src,Nconv);
    assert( Nconv >= LPar.Nvec && "MDistil::LapEvec : Error - not enough eigenvectors converged" );
    if( Nconv > LPar.Nvec )
      eig[t].resize( LPar.Nvec, gridLD );
    RotateEigen( eig[t].evec ); // Rotate the eigenvectors into our phase convention

    for (int i=0;i<LPar.Nvec;i++){
      InsertSliceLocal(eig[t].evec[i],eig4d.evec[i],0,t,3);
      if(t==0)
        eig4d.eval[i] = eig[t].eval[i]; // TODO: Discuss: is this needed? Is there a better way?
    }
  }

  // Now write out the 4d eigenvectors
  eig4d.record.operatorXml = DefaultOperatorXml;
  eig4d.record.solverXml = DefaultsolverXml;
  std::string sEigenPackName(getName());
  sEigenPackName.append(".");
  sEigenPackName.append(std::to_string(vm().getTrajectory()));
  eig4d.write(sEigenPackName,false);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_LapEvec_hpp_
