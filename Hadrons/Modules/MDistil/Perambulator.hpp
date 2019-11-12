/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: Hadrons/Modules/MDistil/Perambulator.hpp
 
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

#ifndef Hadrons_MDistil_Perambulator_hpp_
#define Hadrons_MDistil_Perambulator_hpp_

#include <Hadrons/Modules/MDistil/DistilCommon.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MDistil)

/******************************************************************************
 *                             Perambulator                                    *
 ******************************************************************************/

class PerambulatorPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PerambulatorPar,
                                    std::string, lapevec,
                                    std::string, solver,
                                    std::string, noise,
                                    std::string, PerambFileName,
                                    std::string, UnsmearedSinkFileName,
                                    std::string, UnsmearedSinkMultiFile,
                                    std::string, DistilPar);
};

template <typename FImpl>
class TPerambulator: public Module<PerambulatorPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    // constructor
    TPerambulator(const std::string name);
    // destructor
    virtual ~TPerambulator(void);
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
protected:
    virtual void Cleanup(void);
protected:
    // These variables are created in setup() and freed in Cleanup()
    GridCartesian * grid3d; // Owned by me, so I must delete it
    GridCartesian * grid4d; // Owned by environment (so I won't delete it)
    // Other members
    unsigned int Ls_;
    std::string sLapEvecName;
    std::string sNoiseName;
    std::string DParName;
};

MODULE_REGISTER_TMP(Perambulator, TPerambulator<FIMPL>, MDistil);

/******************************************************************************
 *                 TPerambulator implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TPerambulator<FImpl>::TPerambulator(const std::string name)
: grid3d{nullptr}, grid4d{nullptr}, Module<PerambulatorPar>(name)
{}

// destructor
template <typename FImpl>
TPerambulator<FImpl>::~TPerambulator(void)
{
    Cleanup();
};

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TPerambulator<FImpl>::getInput(void)
{
    sLapEvecName = par().lapevec;
    sNoiseName = par().noise;
    DParName = par().DistilPar;
    return {sLapEvecName, par().solver, sNoiseName, DParName };
}

template <typename FImpl>
std::vector<std::string> TPerambulator<FImpl>::getOutput(void)
{
    return {getName(), getName() + "_unsmeared_sink"};
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPerambulator<FImpl>::setup(void)
{
    Cleanup();
    grid4d = env().getGrid();
    grid3d = MakeLowerDimGrid(grid4d);
    auto &DPar         = envGet(MDistil::DistilParameters,  DParName);
    const int Nt{env().getDim(Tdir)}; 
    const int nvec{DPar.nvec}; 
    const int nnoise{DPar.nnoise}; 
    const int tsrc{DPar.tsrc}; 
    const int TI{DPar.TI}; 
    const int LI{DPar.LI}; 
    const int SI{DPar.SI}; 
    const bool full_tdil{ TI == Nt }; 
    const int Nt_inv{ full_tdil ? 1 : TI };
    const std::string UnsmearedSinkFileName{ par().UnsmearedSinkFileName };
    if (!UnsmearedSinkFileName.empty())
        //bool bMulti = ( Hadrons::MDistil::DistilParameters::ParameterDefault( par().UnsmearedSinkMultiFile, 1, true ) != 0 );
        bool bMulti = 0;
    
    envCreate(PerambTensor, getName(), 1, Nt,nvec,LI,nnoise,Nt_inv,SI);
    envCreate(std::vector<FermionField>, getName() + "_unsmeared_sink", 1,
              nnoise*LI*Ns*Nt_inv, envGetGrid(FermionField));
    
    envTmpLat(LatticeSpinColourVector,   "dist_source");
    envTmpLat(LatticeSpinColourVector,   "source4d");
    envTmp(LatticeSpinColourVector,      "source3d",1,LatticeSpinColourVector(grid3d));
    envTmp(LatticeColourVector,          "source3d_nospin",1,LatticeColourVector(grid3d));
    envTmpLat(LatticeSpinColourVector,   "result4d");
    envTmpLat(LatticeColourVector,       "result4d_nospin");
    envTmp(LatticeColourVector,          "result3d_nospin",1,LatticeColourVector(grid3d));
    envTmp(LatticeColourVector,          "evec3d",1,LatticeColourVector(grid3d));
    
    Ls_ = env().getObjectLs(par().solver);
    envTmpLat(FermionField, "v4dtmp");
    envTmpLat(FermionField, "v5dtmp", Ls_);
    envTmpLat(FermionField, "v5dtmp_sol", Ls_);
}

// clean up any temporaries created by setup (that aren't stored in the environment)
template <typename FImpl>
void TPerambulator<FImpl>::Cleanup(void)
{
    if (grid3d != nullptr)
    {
        delete grid3d;
        grid3d = nullptr;
    }
    grid4d = nullptr;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPerambulator<FImpl>::execute(void)
{
    auto &DPar         = envGet(MDistil::DistilParameters,  DParName);
    const int Nt{env().getDim(Tdir)}; 
    const int nvec{DPar.nvec}; 
    const int nnoise{DPar.nnoise}; 
    const int tsrc{DPar.tsrc}; 
    const int TI{DPar.TI}; 
    const int LI{DPar.LI}; 
    const int SI{DPar.SI}; 
    const bool full_tdil{ TI == Nt }; 
    const int Nt_inv{ full_tdil ? 1 : TI };

    auto &solver=envGet(Solver, par().solver);
    auto &mat = solver.getFMat();
    envGetTmp(FermionField, v4dtmp);
    envGetTmp(FermionField, v5dtmp);
    envGetTmp(FermionField, v5dtmp_sol);
    auto &noise = envGet(NoiseTensor, sNoiseName);
    auto &perambulator = envGet(PerambTensor, getName());
    auto &epack = envGet(LapEvecs, sLapEvecName);
    auto &unsmeared_sink = envGet(std::vector<FermionField>, getName() + "_unsmeared_sink");
    envGetTmp(LatticeSpinColourVector, dist_source);
    envGetTmp(LatticeSpinColourVector, source4d);
    envGetTmp(LatticeSpinColourVector, source3d);
    envGetTmp(LatticeColourVector, source3d_nospin);
    envGetTmp(LatticeSpinColourVector, result4d);
    envGetTmp(LatticeColourVector, result4d_nospin);
    envGetTmp(LatticeColourVector, result3d_nospin);
    envGetTmp(LatticeColourVector, evec3d);
    const int Ntlocal{grid4d->LocalDimensions()[3]};
    const int Ntfirst{grid4d->LocalStarts()[3]};
    const std::string UnsmearedSinkFileName{ par().UnsmearedSinkFileName };

    for (int inoise = 0; inoise < nnoise; inoise++)
    {
        for (int dk = 0; dk < LI; dk++)
        {
            for (int dt = 0; dt < Nt_inv; dt++)
            {
                for (int ds = 0; ds < SI; ds++)
                {
                    LOG(Message) <<  "LapH source vector from noise " << inoise << " and dilution component (d_k,d_t,d_alpha) : (" << dk << ","<< dt << "," << ds << ")" << std::endl;
                    dist_source = 0;
                    source3d_nospin = 0;
                    evec3d = 0;
                    for (int it = dt; it < Nt; it += TI)
                    {
                        const int t_inv{full_tdil ? tsrc : it};
                        if( t_inv >= Ntfirst && t_inv < Ntfirst + Ntlocal )
                        {
                            for (int ik = dk; ik < nvec; ik += LI)
                            {
                                for (int is = ds; is < Ns; is += SI)
                                {
                                    ExtractSliceLocal(evec3d,epack.evec[ik],0,t_inv-Ntfirst,Tdir);
                                    source3d_nospin = evec3d * noise.tensor(inoise, t_inv, ik, is);
                                    source3d=0;
                                    pokeSpin(source3d,source3d_nospin,is);
                                    source4d=0;
                                    InsertSliceLocal(source3d,source4d,0,t_inv-Ntfirst,Tdir);
                                    dist_source += source4d;
                                }
                            }
                        }
                    }
                    result4d=0;
                    v4dtmp = dist_source;
                    if (Ls_ == 1)
                        solver(result4d, v4dtmp);
                    else
                    {
                        mat.ImportPhysicalFermionSource(v4dtmp, v5dtmp);
                        solver(v5dtmp_sol, v5dtmp);
                        mat.ExportPhysicalFermionSolution(v5dtmp_sol, v4dtmp);
                        result4d = v4dtmp;
                    }
                    if (!UnsmearedSinkFileName.empty())
                        unsmeared_sink[inoise+nnoise*(dk+LI*(dt+Nt_inv*ds))] = result4d;
                    for (int is = 0; is < Ns; is++)
                    {
                        result4d_nospin = peekSpin(result4d,is);
                        for (int t = Ntfirst; t < Ntfirst + Ntlocal; t++)
                        {
                            ExtractSliceLocal(result3d_nospin,result4d_nospin,0,t-Ntfirst,Tdir); 
			    for (int ivec = 0; ivec < nvec; ivec++)
                            {
                                ExtractSliceLocal(evec3d,epack.evec[ivec],0,t-Ntfirst,Tdir);
                                pokeSpin(perambulator.tensor(t, ivec, dk, inoise,dt,ds),static_cast<Complex>(innerProduct(evec3d, result3d_nospin)),is);
                            }
                        }
                    }
                }
            }
        }
    }
    // Now share my timeslice data with other members of the grid
    const int NumSlices{grid4d->_processors[Tdir] / grid3d->_processors[Tdir]};
    if (NumSlices > 1)
    {
        LOG(Debug) <<  "Sharing perambulator data with other nodes" << std::endl;
        const int MySlice {grid4d->_processor_coor[Tdir]};
        const int SliceCount {static_cast<int>(perambulator.tensor.size()/NumSlices)};
        PerambTensor::Scalar * const MyData {perambulator.tensor.data()+MySlice*SliceCount};
        Coordinate coor(Nd);
        for (int i = 0 ; i < Tdir ; i++) coor[i] = grid4d->_processor_coor[i];
        std::vector<CommsRequest_t> reqs(0);
        for (int i = 1; i < NumSlices ; i++)
        {
            coor[Tdir] = (MySlice+i)%NumSlices;
            const int SendRank { grid4d->RankFromProcessorCoor(coor) };
            const int RecvSlice { ( MySlice - i + NumSlices ) % NumSlices };
            coor[Tdir] = RecvSlice;
            const auto RecvRank = grid4d->RankFromProcessorCoor(coor);
            grid4d->SendToRecvFromBegin(reqs,MyData,SendRank, perambulator.tensor.data()
                                        + RecvSlice*SliceCount,RecvRank,SliceCount*sizeof(PerambTensor::Scalar));
        }
        grid4d->SendToRecvFromComplete(reqs);
    }
    
    // Save the perambulator to disk from the boss node
    if (grid4d->IsBoss())
    {
        std::string sPerambName {par().PerambFileName};
        sPerambName.append(".");
        sPerambName.append(std::to_string(vm().getTrajectory()));
        perambulator.write(sPerambName.c_str());
    }
    
    //Save the unsmeared sinks if filename specified
    if (!UnsmearedSinkFileName.empty())
    {
       // bool bMulti = ( Hadrons::MDistil::DistilParameters::ParameterDefault( par().UnsmearedSinkMultiFile, 1, false ) != 0 );
        bool bMulti =0;
        LOG(Message) << "Writing unsmeared sink to " << UnsmearedSinkFileName << std::endl;
        A2AVectorsIo::write(UnsmearedSinkFileName, unsmeared_sink, bMulti, vm().getTrajectory());
    }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_Perambulator_hpp_
