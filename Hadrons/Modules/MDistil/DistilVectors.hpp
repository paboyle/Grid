/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: Hadrons/Modules/MDistil/DistilVectors.hpp
 
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

#ifndef Hadrons_MDistil_DistilVectors_hpp_
#define Hadrons_MDistil_DistilVectors_hpp_

#include <Hadrons/Modules/MDistil/DistilCommon.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MDistil)

/******************************************************************************
 *                         DistilVectors                                      *
 *                (Create rho and/or phi vectors)                             *
 ******************************************************************************/

class DistilVectorsPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(DistilVectorsPar,
                                  std::string, noise,
                                  std::string, perambulator,
                                  std::string, lapevec,
                                  std::string, source,
                                  std::string, sink,
                                  int, tsrc,
                                  std::string, nvec,
                                  std::string, TI)
};

template <typename FImpl>
class TDistilVectors: public Module<DistilVectorsPar>
{
public:
  FERM_TYPE_ALIASES(FImpl,);
  // constructor
  TDistilVectors(const std::string name);
  // destructor
  virtual ~TDistilVectors(void);
  // dependency relation
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);
  // setup
  virtual void setup(void);
  // execution
  virtual void execute(void);
protected:
  // These variables are created in setup() and freed in Cleanup()
  GridCartesian * grid3d; // Owned by me, so I must delete it
  GridCartesian * grid4d; // Owned by environment (so I won't delete it)
  virtual void Cleanup(void);
public:
  // These variables contain parameters
  std::string PerambulatorName;
  std::string NoiseVectorName;
  std::string LapEvecName;
  bool bMakeSource;
  bool bMakeSink;
  std::string SourceName;
  std::string SinkName;
};

MODULE_REGISTER_TMP(DistilVectors, TDistilVectors<FIMPL>, MDistil);

/******************************************************************************
 *                 TDistilVectors implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDistilVectors<FImpl>::TDistilVectors(const std::string name)
:  grid3d{nullptr}, grid4d{nullptr}, Module<DistilVectorsPar>(name)
{}
// destructor
template <typename FImpl>
TDistilVectors<FImpl>::~TDistilVectors(void)
{
  Cleanup();
};

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDistilVectors<FImpl>::getInput(void)
{
  PerambulatorName = par().perambulator;
  if( PerambulatorName.size() == 0 ) {
    PerambulatorName = getName();
    PerambulatorName.append( "_peramb" );
  }
  NoiseVectorName = par().noise;
  if( NoiseVectorName.size() == 0 ) {
    NoiseVectorName = PerambulatorName;
    NoiseVectorName.append( "_noise" );
  }
  LapEvecName = par().lapevec;
  if( LapEvecName.size() == 0 ) {
    LapEvecName = PerambulatorName;
    LapEvecName.append( "_lapevec" );
  }
  return { PerambulatorName, NoiseVectorName, LapEvecName };
}

template <typename FImpl>
std::vector<std::string> TDistilVectors<FImpl>::getOutput(void)
{
  SourceName  = par().source;
  SinkName    = par().sink;
  bMakeSource = ( SourceName.size() > 0 );
  bMakeSink   = (   SinkName.size() > 0 );
  if( !bMakeSource && !bMakeSink ) {
    SourceName = getName();
    SinkName   = SourceName;
    SourceName.append( "_rho" );
      SinkName.append( "_phi" );
    bMakeSource = true;
    bMakeSink   = true;
  }
  std::vector<std::string> out;
  if( bMakeSource )
    out.push_back( SourceName );
  if( bMakeSink )
    out.push_back( SinkName );
  return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilVectors<FImpl>::setup(void)
{
  Cleanup();
  auto &noise        = envGet(NoiseTensor, NoiseVectorName);
  auto &perambulator = envGet(PerambTensor, PerambulatorName);

  // We expect the perambulator to have been created with these indices
  assert( perambulator.ValidateIndexNames() && "Perambulator index names bad" );

  const int Nt{ env().getDim(Tdir) };
  assert( Nt == static_cast<int>( perambulator.tensor.dimension(0) ) && "PerambTensor time dimensionality bad" );
  const int TI{ Hadrons::MDistil::DistilParameters::ParameterDefault( par().TI, Nt, true) };
  const int LI{ static_cast<int>( perambulator.tensor.dimension(2) ) };
  const int SI{ static_cast<int>( perambulator.tensor.dimension(5) ) };
  const int Nt_inv{ static_cast<int>( perambulator.tensor.dimension(4) ) };
  const int nnoise{ static_cast<int>( perambulator.tensor.dimension(3) ) };
  assert( nnoise >= static_cast<int>( noise.dimension(0) ) && "Not enough noise vectors for perambulator" );
  // Nvec defaults to what's in the perambulator unless overriden
  const int nvec_per{ static_cast<int>( perambulator.tensor.dimension(1) ) };
  const int nvec{Hadrons::MDistil::DistilParameters::ParameterDefault(par().nvec, nvec_per, true) };
  assert( nvec <= nvec_per && "Not enough distillation sub-space vectors" );

  if( bMakeSource )
    envCreate(std::vector<FermionField>, SourceName, 1, nnoise*LI*SI*Nt_inv, envGetGrid(FermionField));
  if( bMakeSink )
    envCreate(std::vector<FermionField>,   SinkName, 1, nnoise*LI*SI*Nt_inv, envGetGrid(FermionField));

  grid4d = env().getGrid();
  Coordinate latt_size   = GridDefaultLatt();
  Coordinate simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();
  Coordinate simd_layout_3 = GridDefaultSimd(Nd-1, vComplex::Nsimd());
  latt_size[Nd-1] = 1;
  simd_layout_3.push_back( 1 );
  mpi_layout[Nd-1] = 1;
  grid3d = MakeLowerDimGrid(grid4d);
  
  
  envTmp(LatticeSpinColourVector, "tmp2",1,LatticeSpinColourVector(grid4d));
  envTmp(LatticeSpinColourVector, "tmp3d",1,LatticeSpinColourVector(grid3d));
  envTmp(LatticeColourVector, "tmp3d_nospin",1,LatticeColourVector(grid3d));
  envTmp(LatticeSpinColourVector, "sink_tslice",1,LatticeSpinColourVector(grid3d));
  envTmp(LatticeColourVector, "evec3d",1,LatticeColourVector(grid3d));
}

// clean up any temporaries created by setup (that aren't stored in the environment)
template <typename FImpl>
void TDistilVectors<FImpl>::Cleanup(void)
{
  if( grid3d != nullptr ) {
    delete grid3d;
    grid3d = nullptr;
  }
  grid4d = nullptr;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilVectors<FImpl>::execute(void)
{
  auto &noise        = envGet(NoiseTensor, NoiseVectorName);
  auto &perambulator = envGet(PerambTensor, PerambulatorName);
  auto &epack        = envGet(Grid::Hadrons::EigenPack<LatticeColourVector>, LapEvecName);
  
  envGetTmp(LatticeSpinColourVector, tmp2);
  envGetTmp(LatticeSpinColourVector, tmp3d);
  envGetTmp(LatticeColourVector,     tmp3d_nospin);
  envGetTmp(LatticeSpinColourVector, sink_tslice);
  envGetTmp(LatticeColourVector,     evec3d);

  const int Ntlocal{ grid4d->LocalDimensions()[3] };
  const int Ntfirst{ grid4d->LocalStarts()[3] };

  const int Nt{ env().getDim(Tdir) };
  const int TI{ Hadrons::MDistil::DistilParameters::ParameterDefault( par().TI, Nt, false ) };
  const int LI{ static_cast<int>( perambulator.tensor.dimension(2) ) };
  const int SI{ static_cast<int>( perambulator.tensor.dimension(5) ) };
  const int Nt_inv{ static_cast<int>( perambulator.tensor.dimension(4) ) };
  const int nnoise{ static_cast<int>( perambulator.tensor.dimension(3) ) };
  // Nvec defaults to what's in the perambulator unless overriden
  const int nvec{Hadrons::MDistil::DistilParameters::ParameterDefault(par().nvec, static_cast<int>( perambulator.tensor.dimension(1) ), false)};
  const int tsrc{ par().tsrc };
  const bool full_tdil{ TI==Nt };

  int vecindex;
  int t_inv;
  if( bMakeSource ) {
    auto &rho = envGet(std::vector<FermionField>, SourceName);
    for( int inoise = 0; inoise < nnoise; inoise++ ) {
      for( int dk = 0; dk < LI; dk++ ) {
        for( int dt = 0; dt < Nt_inv; dt++ ) {
          for( int ds = 0; ds < SI; ds++ ) {
            vecindex = inoise + nnoise * dk + nnoise * LI * ds + nnoise *LI * SI*dt;
            rho[vecindex] = 0;
            tmp3d_nospin = 0;
            for (int it = dt; it < Nt; it += TI){
              if (full_tdil) t_inv = tsrc; else t_inv = it;
              if( t_inv >= Ntfirst && t_inv < Ntfirst + Ntlocal ) {
                for (int ik = dk; ik < nvec; ik += LI){
                  for (int is = ds; is < Ns; is += SI){
                    ExtractSliceLocal(evec3d,epack.evec[ik],0,t_inv-Ntfirst,Tdir);
                    tmp3d_nospin = evec3d * noise(inoise, t_inv, ik, is);
                    tmp3d=0;
                    pokeSpin(tmp3d,tmp3d_nospin,is);
                    tmp2=0;
                    InsertSliceLocal(tmp3d,tmp2,0,t_inv-Ntfirst,Tdir);
                    rho[vecindex] += tmp2;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if( bMakeSink ) {
    auto &phi = envGet(std::vector<FermionField>, SinkName);
    for( int inoise = 0; inoise < nnoise; inoise++ ) {
      for( int dk = 0; dk < LI; dk++ ) {
        for( int dt = 0; dt < Nt_inv; dt++ ) {
          for( int ds = 0; ds < SI; ds++ ) {
            vecindex = inoise + nnoise * dk + nnoise * LI * ds + nnoise *LI * SI*dt;
            phi[vecindex] = 0;
            for (int t = Ntfirst; t < Ntfirst + Ntlocal; t++) {
              sink_tslice=0;
              for (int ivec = 0; ivec < nvec; ivec++) {
                ExtractSliceLocal(evec3d,epack.evec[ivec],0,t-Ntfirst,Tdir);
                sink_tslice += evec3d * perambulator.tensor(t, ivec, dk, inoise,dt,ds);
              }
              InsertSliceLocal(sink_tslice,phi[vecindex],0,t-Ntfirst,Tdir);
            }
          }
        }
      }
    }
  }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif // Hadrons_MDistil_DistilVectors_hpp_
