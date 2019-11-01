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
                                  double, rho) 
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
                                  int, Nvec,
                                  int, Nk,
                                  int, Np,
                                  int, MaxIt,
                                  double, resid,
                                  int, IRLLog)
  LanczosParameters() = default;
  template <class ReaderClass> LanczosParameters(Reader<ReaderClass>& Reader){read(Reader,"Lanczos",*this);}
};

// These are the actual parameters passed to the module during construction

struct LapEvecPar: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(LapEvecPar
                                  ,std::string,         gauge
                                  ,StoutParameters,     Stout
                                  ,ChebyshevParameters, Cheby
                                  ,LanczosParameters,   Lanczos)
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
  std::string sGaugeName;
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
  sGaugeName = par().gauge;
  if( sGaugeName.size() == 0 ) {
    sGaugeName = getName();
    sGaugeName.append( "_gauge" );
  }
  return std::vector<std::string>{ sGaugeName };
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
  const int Ntlocal{gridHD->LocalDimensions()[Tdir]};
  // Temporaries
  envTmpLat(GaugeField, "Umu_stout");
  envTmpLat(GaugeField, "Umu_smear");
  envTmp(LatticeGaugeField, "UmuNoTime",1,LatticeGaugeField(gridLD));
  envTmp(LatticeColourVector, "src",1,LatticeColourVector(gridLD));
  envTmp(std::vector<LapEvecs>, "eig",1,std::vector<LapEvecs>(Ntlocal));
  // Output objects
  envCreate(LapEvecs, getName(), 1, par().Lanczos.Nvec, gridHD );
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

/*************************************************************************************

 -Grad^2 (Peardon, 2009, pg 2, equation 3, https://arxiv.org/abs/0905.2160)
 Field      Type of field the operator will be applied to
 GaugeField Gauge field the operator will smear using
 
 *************************************************************************************/

template<typename Field, typename GaugeField=LatticeGaugeField>
class Laplacian3D : public LinearOperatorBase<Field>, public LinearFunction<Field> {
  typedef typename GaugeField::vector_type vCoeff_t;
protected: // I don't really mind if _gf is messed with ... so make this public?
  //GaugeField & _gf;
  int          nd; // number of spatial dimensions
  std::vector<Lattice<iColourMatrix<vCoeff_t> > > U;
public:
  // Construct this operator given a gauge field and the number of dimensions it should act on
  Laplacian3D( GaugeField& gf, int dimSpatial = Tdir ) : /*_gf(gf),*/ nd{dimSpatial} {
    assert(dimSpatial>=1);
    for( int mu = 0 ; mu < nd ; mu++ )
      U.push_back(PeekIndex<LorentzIndex>(gf,mu));
      }
  
  // Apply this operator to "in", return result in "out"
  void operator()(const Field& in, Field& out) {
    assert( nd <= in.Grid()->Nd() );
    conformable( in, out );
    out = ( ( Real ) ( 2 * nd ) ) * in;
    Field _tmp(in.Grid());
    typedef typename GaugeField::vector_type vCoeff_t;
    //Lattice<iColourMatrix<vCoeff_t> > U(in.Grid());
    for( int mu = 0 ; mu < nd ; mu++ ) {
      //U = PeekIndex<LorentzIndex>(_gf,mu);
      out -= U[mu] * Cshift( in, mu, 1);
      _tmp = adj( U[mu] ) * in;
      out -= Cshift(_tmp,mu,-1);
    }
  }
  
  void OpDiag (const Field &in, Field &out) { assert(0); };
  void OpDir  (const Field &in, Field &out,int dir,int disp) { assert(0); };
  void Op     (const Field &in, Field &out) { assert(0); };
  void AdjOp  (const Field &in, Field &out) { assert(0); };
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2) { assert(0); };
  void HermOp(const Field &in, Field &out) { operator()(in,out); };
};

template<typename Field>
class Laplacian3DHerm : public LinearFunction<Field> {
public:
  OperatorFunction<Field>   & _poly;
  LinearOperatorBase<Field> &_Linop;
  
  Laplacian3DHerm(OperatorFunction<Field> & poly,LinearOperatorBase<Field>& linop)
  : _poly{poly}, _Linop{linop} {}
  
  void operator()(const Field& in, Field& out) {
    _poly(_Linop,in,out);
  }
};

/******************************************************************************
 Calculate low-mode eigenvalues of the Laplacian
 ******************************************************************************/

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TLapEvec<GImpl>::execute(void)
{
  const ChebyshevParameters &ChebPar{par().Cheby};
  const LanczosParameters   &LPar{par().Lanczos};

  // Disable IRL logging if requested
  LOG(Message) << "IRLLog=" << LPar.IRLLog << std::endl;
  const int PreviousIRLLogState{GridLogIRL.isActive()};
  GridLogIRL.Active( LPar.IRLLog == 0 ? 0 : 1 );

  // Stout smearing
  envGetTmp(GaugeField, Umu_smear);
  Umu_smear = envGet(GaugeField, sGaugeName); // The smeared field starts off as the Gauge field
  LOG(Message) << "Initial plaquette: " << WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu_smear) << std::endl;
  const StoutParameters &Stout{par().Stout};
  if( Stout.steps )
  {
    envGetTmp(GaugeField, Umu_stout);
    Smear_Stout<PeriodicGimplR> LS(Stout.rho, Tdir); // spatial smearing only
    for (int i = 0; i < Stout.steps; i++) {
      LS.smear(Umu_stout, Umu_smear);
      Umu_smear = Umu_stout;
    }
    LOG(Message) << "Smeared plaquette: " << WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu_smear) << std::endl;
  }

  ////////////////////////////////////////////////////////////////////////
  // Invert nabla operator separately on each time-slice
  ////////////////////////////////////////////////////////////////////////
  
  auto & eig4d = envGet(LapEvecs, getName() );
  envGetTmp(std::vector<LapEvecs>, eig);   // Eigenpack for each timeslice
  envGetTmp(LatticeGaugeField, UmuNoTime); // Gauge field without time dimension
  envGetTmp(LatticeColourVector, src);
  const int Ntlocal{gridHD->LocalDimensions()[Tdir]};
  const int Ntfirst{gridHD->LocalStarts()[Tdir]};
  uint32_t ConvergenceErrors{0};
  for(int t = 0; t < Ntlocal; t++ ) {
    LOG(Message) << "------------------------------------------------------------" << std::endl;
    LOG(Message) << " Compute eigenpack, local timeslice = " << t << " / " << Ntlocal << std::endl;
    LOG(Message) << "------------------------------------------------------------" << std::endl;
    eig[t].resize(LPar.Nk+LPar.Np,gridLD);
    
    // Construct smearing operator
    ExtractSliceLocal(UmuNoTime,Umu_smear,0,t,Tdir); // switch to 3d/4d objects
    Laplacian3D<LatticeColourVector> Nabla(UmuNoTime);
    LOG(Debug) << "Chebyshev preconditioning to order " << ChebPar.PolyOrder
      << " with parameters (alpha,beta) = (" << ChebPar.alpha << "," << ChebPar.beta << ")" << std::endl;
    Chebyshev<LatticeColourVector> Cheb(ChebPar.alpha,ChebPar.beta,ChebPar.PolyOrder);
    
    // Construct source vector according to Test_dwf_compressed_lanczos.cc
    src = 11.0;  //TODO: Why hard-coded 11?
    RealD nn = norm2(src);
    nn = Grid::sqrt(nn);
    src = src * (1.0/nn);

    Laplacian3DHerm<LatticeColourVector> NablaCheby(Cheb,Nabla);
    ImplicitlyRestartedLanczos<LatticeColourVector>
      IRL(NablaCheby,Nabla,LPar.Nvec,LPar.Nk,LPar.Nk+LPar.Np,LPar.resid,LPar.MaxIt);
    int Nconv = 0;
    IRL.calc(eig[t].eval,eig[t].evec,src,Nconv);
    if( Nconv < LPar.Nvec ) {
      // NB: Can't assert here since we are processing local slices - i.e. not all nodes would assert
      ConvergenceErrors = 1;
      LOG(Error) << "MDistil::LapEvec : Not enough eigenvectors converged. If this occurs in practice, we should modify the eigensolver to iterate once more to ensure the second convergence test does not take us below the requested number of eigenvectors" << std::endl;
    }
    if( Nconv != LPar.Nvec )
      eig[t].resize( LPar.Nvec, gridLD );
    RotateEigen( eig[t].evec ); // Rotate the eigenvectors into our phase convention

    for (int i=0;i<LPar.Nvec;i++){
      InsertSliceLocal(eig[t].evec[i],eig4d.evec[i],0,t,Tdir);
      if(t==0 && Ntfirst==0)
        eig4d.eval[i] = eig[t].eval[i]; // TODO: Discuss: is this needed? Is there a better way?
    }
  }
  GridLogIRL.Active( PreviousIRLLogState );
  gridHD->GlobalSum(ConvergenceErrors);
  assert(ConvergenceErrors==0 && "The eingensolver failed to find enough eigenvectors on at least one node");
#if DEBUG
  // Now write out the 4d eigenvectors
  eig4d.record.operatorXml = "<OPERATOR>Distillation</OPERATOR>";
  eig4d.record.solverXml = "<SOLVER>CG</SOLVER>";
  std::string sEigenPackName(getName());
  sEigenPackName.append(".");
  sEigenPackName.append(std::to_string(vm().getTrajectory()));
  eig4d.write(sEigenPackName,false);
#endif
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_LapEvec_hpp_
