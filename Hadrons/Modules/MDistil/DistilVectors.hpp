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
                                    std::string, rho,
                                    std::string, phi,
                                    std::string, DistilPar);
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
    std::string DParName;
    bool bMakeRho;
    bool bMakePhi;
    std::string RhoName;
    std::string PhiName;
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
    NoiseVectorName = par().noise;
    LapEvecName = par().lapevec;
    DParName = par().DistilPar;
    return { PerambulatorName, NoiseVectorName, LapEvecName, DParName };
}

template <typename FImpl>
std::vector<std::string> TDistilVectors<FImpl>::getOutput(void)
{
    RhoName  = par().rho;
    PhiName  = par().phi;
    bMakeRho = ( RhoName.size() > 0 );
    bMakePhi = ( PhiName.size() > 0 );
    if (!bMakeRho && !bMakePhi)
    {
        RhoName = getName();
        PhiName = RhoName;
        RhoName.append("_rho");
        PhiName.append("_phi");
        bMakeRho = true;
        bMakePhi = true;
    }
    std::vector<std::string> out;
    if (bMakeRho)
        out.push_back(RhoName);
    if (bMakePhi)
        out.push_back(PhiName);
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilVectors<FImpl>::setup(void)
{
    Cleanup();
    auto &noise        = envGet(NoiseTensor,  NoiseVectorName);
    auto &perambulator = envGet(PerambTensor, PerambulatorName);
    auto &DPar         = envGet(MDistil::DistilParameters,  DParName);
    
    // We expect the perambulator to have been created with these indices
    assert( perambulator.ValidateIndexNames() && "Perambulator index names bad" );
    
    const int Nt{env().getDim(Tdir)}; 
    const int nvec{DPar.nvec}; 
    const int nnoise{DPar.nnoise}; 
    const int TI{DPar.TI}; 
    const int LI{DPar.LI}; 
    const int SI{DPar.SI}; 
    const bool full_tdil{ TI == Nt }; 
    const int Nt_inv{ full_tdil ? 1 : TI };
    
    if (bMakeRho)
        envCreate(std::vector<FermionField>, RhoName, 1, nnoise*LI*SI*Nt_inv, envGetGrid(FermionField));
    if (bMakePhi)
        envCreate(std::vector<FermionField>,   PhiName, 1, nnoise*LI*SI*Nt_inv, envGetGrid(FermionField));
    
    grid4d = env().getGrid();
    Coordinate latt_size   = GridDefaultLatt();
    Coordinate mpi_layout  = GridDefaultMpi();
    Coordinate simd_layout_3 = GridDefaultSimd(Nd-1, vComplex::Nsimd());
    latt_size[Nd-1] = 1;
    simd_layout_3.push_back( 1 );
    mpi_layout[Nd-1] = 1;
    grid3d = MakeLowerDimGrid(grid4d);
    
    envTmp(LatticeSpinColourVector, "source4d",1,LatticeSpinColourVector(grid4d));
    envTmp(LatticeSpinColourVector, "source3d",1,LatticeSpinColourVector(grid3d));
    envTmp(LatticeColourVector, "source3d_nospin",1,LatticeColourVector(grid3d));
    envTmp(LatticeSpinColourVector, "sink3d",1,LatticeSpinColourVector(grid3d));
    envTmp(LatticeColourVector, "evec3d",1,LatticeColourVector(grid3d));
}

// clean up any temporaries created by setup (that aren't stored in the environment)
template <typename FImpl>
void TDistilVectors<FImpl>::Cleanup(void)
{
    if ( grid3d != nullptr)
    {
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
    auto &DPar         = envGet(MDistil::DistilParameters,  DParName);
    
    envGetTmp(LatticeSpinColourVector, source4d);
    envGetTmp(LatticeSpinColourVector, source3d);
    envGetTmp(LatticeColourVector,     source3d_nospin);
    envGetTmp(LatticeSpinColourVector, sink3d);
    envGetTmp(LatticeColourVector,     evec3d);
    
    const int Ntlocal{ grid4d->LocalDimensions()[3] };
    const int Ntfirst{ grid4d->LocalStarts()[3] };
    
    const int Nt{env().getDim(Tdir)}; 
    const int nvec{DPar.nvec}; 
    const int nnoise{DPar.nnoise}; 
    const int tsrc{DPar.tsrc}; 
    const int TI{DPar.TI}; 
    const int LI{DPar.LI}; 
    const int SI{DPar.SI}; 
    const bool full_tdil{ TI == Nt }; 
    const int Nt_inv{ full_tdil ? 1 : TI };
    
    int vecindex;
    int t_inv;
    if (bMakeRho)
    {
        auto &rho = envGet(std::vector<FermionField>, RhoName);
        for (int inoise = 0; inoise < nnoise; inoise++) {
            for (int dk = 0; dk < LI; dk++) {
                for (int dt = 0; dt < Nt_inv; dt++) {
                    for (int ds = 0; ds < SI; ds++) {
                        vecindex = inoise + nnoise * dk + nnoise * LI * ds + nnoise *LI * SI*dt;
                        rho[vecindex] = 0;
                        source3d_nospin = 0;
                        for (int it = dt; it < Nt; it += TI){
                            if (full_tdil) t_inv = tsrc; else t_inv = it;
                            if (t_inv >= Ntfirst && t_inv < Ntfirst + Ntlocal) {
                                for (int ik = dk; ik < nvec; ik += LI){
                                    for (int is = ds; is < Ns; is += SI){
                                        ExtractSliceLocal(evec3d,epack.evec[ik],0,t_inv-Ntfirst,Tdir);
                                        source3d_nospin = evec3d * noise.tensor(inoise, t_inv, ik, is);
                                        source3d=0;
                                        pokeSpin(source3d,source3d_nospin,is);
                                        source4d=0;
                                        InsertSliceLocal(source3d,source4d,0,t_inv-Ntfirst,Tdir);
                                        rho[vecindex] += source4d;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if (bMakePhi) {
        auto &phi = envGet(std::vector<FermionField>, PhiName);
        for (int inoise = 0; inoise < nnoise; inoise++) {
            for (int dk = 0; dk < LI; dk++) {
                for (int dt = 0; dt < Nt_inv; dt++) {
                    for (int ds = 0; ds < SI; ds++) {
                        vecindex = inoise + nnoise * dk + nnoise * LI * ds + nnoise *LI * SI*dt;
                        phi[vecindex] = 0;
                        for (int t = Ntfirst; t < Ntfirst + Ntlocal; t++) {
                            sink3d=0;
                            for (int ivec = 0; ivec < nvec; ivec++) {
                                ExtractSliceLocal(evec3d,epack.evec[ivec],0,t-Ntfirst,Tdir);
                                sink3d += evec3d * perambulator.tensor(t, ivec, dk, inoise,dt,ds);
                            }
                            InsertSliceLocal(sink3d,phi[vecindex],0,t-Ntfirst,Tdir);
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
